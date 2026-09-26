import json
import os
import xml.etree.ElementTree as ET

from literature_search import EMAIL, TOOL
from taxonomy_aliases import ncbi_get, normalize_name


CACHE_FILE = os.environ.get(
    "VIRUS_TAXONOMY_CACHE",
    "virus_taxonomy_cache_correctness.json",
)


# ============================================================
# CACHE
# ============================================================


def load_cache():
    if not os.path.exists(CACHE_FILE):
        return {}

    try:
        with open(
            CACHE_FILE,
            "r",
            encoding="utf-8",
        ) as handle:
            return json.load(handle)
    except Exception:
        return {}



def save_cache(cache):
    with open(
        CACHE_FILE,
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            cache,
            handle,
            indent=2,
            ensure_ascii=False,
            sort_keys=True,
        )


# ============================================================
# NCBI TAXONOMY PARSING
# ============================================================


def _all_name_texts(taxon):
    names = set()

    scientific_name = (
        taxon.findtext("ScientificName")
        or ""
    ).strip()

    if scientific_name:
        names.add(scientific_name)

    other_names = taxon.find("OtherNames")

    if other_names is not None:
        accepted_tags = {
            "Synonym",
            "EquivalentName",
            "GenbankCommonName",
            "CommonName",
            "Acronym",
            "GenbankSynonym",
        }

        for child in other_names:
            tag = child.tag.split("}")[-1]
            text = (child.text or "").strip()

            if tag in accepted_tags and text:
                names.add(text)

    return names



def _parse_taxon(taxon, query):
    if taxon is None:
        return {
            "query": query,
            "tax_id": None,
            "scientific_name": None,
            "rank": None,
            "aliases": [query] if query else [],
            "lineage": [],
            "by_rank": {},
            "resolved": False,
        }

    aliases = _all_name_texts(taxon)

    if query:
        aliases.add(query)

    lineage = []
    by_rank = {}

    lineage_ex = taxon.find("LineageEx")

    if lineage_ex is not None:
        for node in lineage_ex.findall("Taxon"):
            item = {
                "tax_id": node.findtext("TaxId"),
                "scientific_name": node.findtext("ScientificName"),
                "rank": node.findtext("Rank") or "no rank",
            }
            lineage.append(item)

            rank = item["rank"]
            name = item["scientific_name"]

            if rank and rank != "no rank" and name:
                by_rank[rank] = name

    scientific_name = taxon.findtext("ScientificName") or ""
    rank = taxon.findtext("Rank") or "no rank"

    lineage.append({
        "tax_id": taxon.findtext("TaxId"),
        "scientific_name": scientific_name,
        "rank": rank,
    })

    if rank and rank != "no rank" and scientific_name:
        by_rank[rank] = scientific_name

    aliases = sorted(
        aliases,
        key=lambda value: (
            0 if normalize_name(value) == normalize_name(query) else 1,
            -len(value),
            value.lower(),
        ),
    )

    return {
        "query": query,
        "tax_id": taxon.findtext("TaxId"),
        "scientific_name": scientific_name or None,
        "rank": rank,
        "aliases": aliases,
        "lineage": lineage,
        "by_rank": by_rank,
        "resolved": True,
    }


# ============================================================
# LOOKUP
# ============================================================


def _fetch_by_tax_id(tax_id, query):
    response = ncbi_get(
        "efetch.fcgi",
        {
            "db": "taxonomy",
            "id": str(tax_id),
            "retmode": "xml",
            "email": EMAIL,
            "tool": TOOL,
        },
    )

    root = ET.fromstring(response.text)
    taxon = root.find(".//Taxon")

    if taxon is not None and normalize_name(query) not in {normalize_name(n) for n in _all_name_texts(taxon)}:
        return _parse_taxon(None, query)

    return _parse_taxon(
        taxon,
        query,
    )



def _fetch_by_name(query):
    response = ncbi_get(
        "esearch.fcgi",
        {
            "db": "taxonomy",
            "term": f'"{query}"',
            "retmode": "json",
            "retmax": 10,
            "email": EMAIL,
            "tool": TOOL,
        },
    )

    ids = (
        response
        .json()
        .get("esearchresult", {})
        .get("idlist", [])
    )

    if not ids:
        return _parse_taxon(None, query)

    response = ncbi_get(
        "efetch.fcgi",
        {
            "db": "taxonomy",
            "id": ",".join(ids),
            "retmode": "xml",
            "email": EMAIL,
            "tool": TOOL,
        },
    )

    root = ET.fromstring(response.text)
    target_norm = normalize_name(query)
    chosen = None

    # Prefer a scientific-name OR synonym match to the query.
    for taxon in root.findall("./Taxon"):
        candidate_names = _all_name_texts(taxon)

        if any(
            normalize_name(name) == target_norm
            for name in candidate_names
        ):
            chosen = taxon
            break

    if chosen is None:
        return _parse_taxon(None, query)

    return _parse_taxon(chosen, query)



def get_virus_taxonomy_context(
    virus,
    tax_id=None,
    force_refresh=False,
):
    """
    Resolve a target virus to NCBI Taxonomy and collect accepted names,
    synonyms, historical/equivalent names, lineage and taxonomic rank.

    A supplied tax_id is preferred because virus names change over time.
    """

    cache = load_cache()

    cache_key = (
        f"taxid:{tax_id}"
        if tax_id
        else f"name:{normalize_name(virus)}"
    )

    if (
        not force_refresh
        and cache_key in cache
    ):
        return cache[cache_key]

    context = None

    if tax_id:
        try:
            context = _fetch_by_tax_id(
                tax_id,
                virus,
            )
        except Exception as error:
            print(
                "WARNING: Virus taxonomy lookup by tax ID "
                f"failed for {virus}: {error}"
            )

    if (
        context is None
        or not context.get("resolved")
    ):
        try:
            context = _fetch_by_name(virus)
        except Exception as error:
            print(
                "WARNING: Virus taxonomy lookup by name "
                f"failed for {virus}: {error}"
            )
            context = _parse_taxon(None, virus)

    cache[cache_key] = context
    save_cache(cache)

    return context



def get_virus_aliases(
    virus,
    tax_id=None,
    force_refresh=False,
):
    context = get_virus_taxonomy_context(
        virus,
        tax_id=tax_id,
        force_refresh=force_refresh,
    )

    aliases = list(
        context.get("aliases", [])
        or []
    )

    if virus and all(
        normalize_name(alias)
        != normalize_name(virus)
        for alias in aliases
    ):
        aliases.insert(0, virus)

    return aliases
