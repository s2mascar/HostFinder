import csv
import json
import os
import re
import xml.etree.ElementTree as ET

from collections import Counter

from literature_search import EMAIL, TOOL

from taxonomy_aliases import (
    ncbi_get,
    normalize_name,
)


VIRUS_HOST_DB = os.environ.get(
    "VIRUS_HOST_DB",
    "data/virushostdb.daily.tsv",
)

CACHE_FILE = os.environ.get(
    "BIOLOGICAL_CONTEXT_CACHE",
    "biological_context_cache_v05.json",
)


_VH_ROWS = None
_VH_NAME_INDEX = {}
_VH_TAXID_INDEX = {}
_SOURCE_TAXONOMY_CACHE = {}


# ============================================================
# NORMALIZATION
# ============================================================


def normalize_header(text):
    text = str(text or "").lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_")


# ============================================================
# BROAD BIOLOGICAL GROUP
# ============================================================


def broad_group_from_lineage(lineage):
    text = str(lineage or "").lower()

    if "viridiplantae" in text:
        return "PLANT"

    if "metazoa" in text:
        return "ANIMAL"

    if "fungi" in text:
        return "FUNGUS"

    if "bacteria" in text:
        return "BACTERIA"

    if "archaea" in text:
        return "ARCHAEA"

    if "viruses" in text:
        return "VIRUS"

    if "eukaryota" in text:
        return "OTHER_EUKARYOTE"

    return "UNKNOWN"


# ============================================================
# NCBI TAXONOMY
# ============================================================


def empty_taxonomy_context(name):
    return {
        "query": name,
        "tax_id": None,
        "scientific_name": None,
        "rank": None,
        "lineage": [],
        "by_rank": {},
        "broad_group": "UNKNOWN",
        "exact_name_match": False,
    }


def get_taxonomy_context(name):
    """Retrieve one NCBI taxonomy record and its ranked lineage."""

    if not name:
        return empty_taxonomy_context(name)

    search_params = {
        "db": "taxonomy",
        "term": f'"{name}"',
        "retmode": "json",
        "retmax": 5,
        "email": EMAIL,
        "tool": TOOL,
    }

    response = ncbi_get(
        "esearch.fcgi",
        search_params,
    )

    ids = (
        response
        .json()
        .get("esearchresult", {})
        .get("idlist", [])
    )

    if not ids:
        return empty_taxonomy_context(name)

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
    target_norm = normalize_name(name)

    chosen = None
    exact_match = False

    for taxon in root.findall(".//Taxon"):
        scientific_name = (
            taxon.findtext("ScientificName")
            or ""
        )

        if normalize_name(scientific_name) == target_norm:
            chosen = taxon
            exact_match = True
            break

    if chosen is None:
        chosen = root.find(".//Taxon")

    if chosen is None:
        return empty_taxonomy_context(name)

    lineage = []
    lineage_ex = chosen.find("LineageEx")

    if lineage_ex is not None:
        for node in lineage_ex.findall("Taxon"):
            lineage.append({
                "tax_id": node.findtext("TaxId"),
                "scientific_name": node.findtext("ScientificName"),
                "rank": node.findtext("Rank") or "no rank",
            })

    lineage.append({
        "tax_id": chosen.findtext("TaxId"),
        "scientific_name": chosen.findtext("ScientificName"),
        "rank": chosen.findtext("Rank") or "no rank",
    })

    by_rank = {}

    for node in lineage:
        rank = node.get("rank")
        scientific_name = node.get("scientific_name")

        if (
            rank
            and rank != "no rank"
            and scientific_name
        ):
            by_rank[rank] = scientific_name

    lineage_string = "; ".join(
        node["scientific_name"]
        for node in lineage
        if node.get("scientific_name")
    )

    return {
        "query": name,
        "tax_id": chosen.findtext("TaxId"),
        "scientific_name": chosen.findtext("ScientificName"),
        "rank": chosen.findtext("Rank") or "no rank",
        "lineage": lineage,
        "by_rank": by_rank,
        "broad_group": broad_group_from_lineage(lineage_string),
        "exact_name_match": exact_match,
    }


# ============================================================
# TAXONOMY LOOKUP BY TAXONOMY ID
# ============================================================


def get_taxonomy_context_by_taxid(tax_id):
    """Resolve a source-organism taxonomy ID from Virus-Host DB."""

    global _SOURCE_TAXONOMY_CACHE

    tax_id = str(tax_id or "").strip()

    if (
        not tax_id
        or not tax_id.isdigit()
        or tax_id == "1"
    ):
        return empty_taxonomy_context(tax_id)

    if tax_id in _SOURCE_TAXONOMY_CACHE:
        return _SOURCE_TAXONOMY_CACHE[tax_id]

    try:
        response = ncbi_get(
            "efetch.fcgi",
            {
                "db": "taxonomy",
                "id": tax_id,
                "retmode": "xml",
                "email": EMAIL,
                "tool": TOOL,
            },
        )

        root = ET.fromstring(response.text)
        chosen = root.find(".//Taxon")

        if chosen is None:
            result = empty_taxonomy_context(tax_id)
        else:
            lineage = []
            lineage_ex = chosen.find("LineageEx")

            if lineage_ex is not None:
                for node in lineage_ex.findall("Taxon"):
                    lineage.append({
                        "tax_id": node.findtext("TaxId"),
                        "scientific_name": node.findtext("ScientificName"),
                        "rank": node.findtext("Rank") or "no rank",
                    })

            lineage.append({
                "tax_id": chosen.findtext("TaxId"),
                "scientific_name": chosen.findtext("ScientificName"),
                "rank": chosen.findtext("Rank") or "no rank",
            })

            by_rank = {}

            for node in lineage:
                rank = node.get("rank")
                name = node.get("scientific_name")

                if rank and rank != "no rank" and name:
                    by_rank[rank] = name

            lineage_string = "; ".join(
                node["scientific_name"]
                for node in lineage
                if node.get("scientific_name")
            )

            result = {
                "query": tax_id,
                "tax_id": chosen.findtext("TaxId"),
                "scientific_name": chosen.findtext("ScientificName"),
                "rank": chosen.findtext("Rank") or "no rank",
                "lineage": lineage,
                "by_rank": by_rank,
                "broad_group": broad_group_from_lineage(lineage_string),
                "exact_name_match": True,
            }

    except Exception as error:
        print(
            "WARNING: Source-organism taxonomy lookup failed:",
            tax_id,
            error,
        )
        result = empty_taxonomy_context(tax_id)

    _SOURCE_TAXONOMY_CACHE[tax_id] = result
    return result


# ============================================================
# VIRUS-HOST DB
# ============================================================


def get_field(row, *possible_names):
    normalized = {
        normalize_header(key): value
        for key, value in row.items()
    }

    for name in possible_names:
        key = normalize_header(name)

        if key in normalized:
            return str(normalized[key] or "").strip()

    return ""


def canonicalize_vh_row(row):
    """
    Virus-Host DB columns documented in its README include:
    virus tax id, virus name, virus lineage, host tax id,
    host name, host lineage, pmid, evidence, sample type,
    and source organism.
    """

    return {
        "virus_tax_id": get_field(
            row,
            "virus tax id",
            "virus_tax_id",
        ),
        "virus_name": get_field(
            row,
            "virus name",
            "virus_name",
        ),
        "virus_lineage": get_field(
            row,
            "virus lineage",
            "virus_lineage",
        ),
        "host_tax_id": get_field(
            row,
            "host tax id",
            "host_tax_id",
        ),
        "host_name": get_field(
            row,
            "host name",
            "host_name",
        ),
        "host_lineage": get_field(
            row,
            "host lineage",
            "host_lineage",
        ),
        "pmid": get_field(
            row,
            "pmid",
        ),
        "evidence": get_field(
            row,
            "evidence",
        ),
        "sample_type": get_field(
            row,
            "sample type",
            "sample_type",
        ),
        "source_organism": get_field(
            row,
            "source organism",
            "source_organism",
        ),
    }


def is_environmental_host_record(row):
    """
    Virus-Host DB can use taxonomy ID 1 / name "root" for viral
    sequences whose source is environmental rather than a resolved
    biological host. Such records must not be treated as known hosts or
    used for taxonomic-distance calculations.
    """

    host_tax_id = str(
        row.get("host_tax_id")
        or ""
    ).strip()

    host_name = normalize_name(
        row.get("host_name")
        or ""
    )

    placeholder_names = {
        "root",
        "environmental sample",
        "environmental samples",
        "environmental samples metagenomes",
        "uncultured organism",
        "unclassified sequences",
    }

    return (
        host_tax_id == "1"
        or host_name in placeholder_names
    )


def load_virus_host_db():
    global _VH_ROWS
    global _VH_NAME_INDEX
    global _VH_TAXID_INDEX

    if _VH_ROWS is not None:
        return

    _VH_ROWS = []
    _VH_NAME_INDEX = {}
    _VH_TAXID_INDEX = {}

    if not os.path.exists(VIRUS_HOST_DB):
        print(
            "WARNING: Virus-Host DB not found at "
            f"{VIRUS_HOST_DB}. Biological context will use "
            "NCBI taxonomy only."
        )
        return

    with open(
        VIRUS_HOST_DB,
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for raw_row in reader:
            row = canonicalize_vh_row(raw_row)

            if not row["virus_name"]:
                continue

            _VH_ROWS.append(row)

            name_key = normalize_name(
                row["virus_name"]
            )

            _VH_NAME_INDEX.setdefault(
                name_key,
                [],
            ).append(row)

            tax_id = row["virus_tax_id"]

            if tax_id:
                _VH_TAXID_INDEX.setdefault(
                    str(tax_id),
                    [],
                ).append(row)

    print(
        "Loaded Virus-Host DB records:",
        len(_VH_ROWS),
    )


# ============================================================
# TAXONOMIC DISTANCE
# ============================================================


RANK_DISTANCE = {
    "species": 0,
    "subspecies": 0,
    "genus": 1,
    "family": 2,
    "order": 3,
    "class": 4,
    "phylum": 5,
    "kingdom": 6,
    "superkingdom": 7,
}


def host_taxonomic_distance(
    target_taxonomy,
    other_host_name,
    other_lineage,
):
    """
    Return a simple rank-based distance between the target host and a
    known host. The lowest common ancestor is inferred by matching the
    target's NCBI lineage names against the Virus-Host DB lineage.
    """

    target_name = normalize_name(
        target_taxonomy.get("scientific_name")
    )

    if (
        target_name
        and target_name == normalize_name(other_host_name)
    ):
        return {
            "distance": 0,
            "lca_rank": "species",
            "lca_name": target_taxonomy.get("scientific_name"),
        }

    other_names = {
        normalize_name(part)
        for part in str(other_lineage or "").split(";")
        if part.strip()
    }

    if other_host_name:
        other_names.add(
            normalize_name(other_host_name)
        )

    for node in reversed(
        target_taxonomy.get("lineage", [])
    ):
        name = node.get("scientific_name")
        rank = node.get("rank")

        if normalize_name(name) in other_names:
            return {
                "distance": RANK_DISTANCE.get(rank, 8),
                "lca_rank": rank,
                "lca_name": name,
            }

    return {
        "distance": 9,
        "lca_rank": None,
        "lca_name": None,
    }


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
        )


# ============================================================
# CONTEXT AGENT
# ============================================================


def run_biological_context_agent(
    host,
    virus,
    force_refresh=False,
):
    """
    Build a biological prior for a host-virus pair.

    IMPORTANT:
    The exact target-host record is removed from the context shown to
    downstream agents so that benchmark labels cannot leak through a
    curated interaction database.
    """

    cache = load_cache()
    cache_key = (
        normalize_name(host)
        + "||"
        + normalize_name(virus)
    )

    if (
        not force_refresh
        and cache_key in cache
    ):
        print(
            "Biological context loaded from cache: ",
            f"{host} / {virus}",
        )
        return cache[cache_key]

    print(
        "Building biological context for:",
        host,
        "/",
        virus,
    )

    try:
        host_taxonomy = get_taxonomy_context(host)
    except Exception as error:
        print(
            "WARNING: Host taxonomy lookup failed:",
            error,
        )
        host_taxonomy = empty_taxonomy_context(host)

    try:
        virus_taxonomy = get_taxonomy_context(virus)
    except Exception as error:
        print(
            "WARNING: Virus taxonomy lookup failed:",
            error,
        )
        virus_taxonomy = empty_taxonomy_context(virus)

    load_virus_host_db()

    records = list(
        _VH_NAME_INDEX.get(
            normalize_name(virus),
            [],
        )
    )

    match_level = (
        "EXACT_VIRUS"
        if records
        else None
    )

    # Try exact virus taxonomy ID if name lookup did not match.
    if (
        not records
        and virus_taxonomy.get("tax_id")
    ):
        records = list(
            _VH_TAXID_INDEX.get(
                str(virus_taxonomy["tax_id"]),
                [],
            )
        )

        if records:
            match_level = "VIRUS_TAXID"

    # Fall back to the viral family as a weaker host-range prior.
    if not records and _VH_ROWS:
        family = (
            virus_taxonomy
            .get("by_rank", {})
            .get("family")
        )

        if family:
            family_norm = normalize_name(family)

            records = [
                row
                for row in _VH_ROWS
                if family_norm in normalize_name(
                    row["virus_lineage"]
                )
            ]

            if records:
                match_level = "VIRUS_FAMILY"

    # ========================================================
    # SEPARATE REAL HOSTS FROM ENVIRONMENTAL PLACEHOLDERS
    # ========================================================

    unique_hosts = {}
    environmental_records = []
    environmental_seen = set()

    for row in records:

        if is_environmental_host_record(row):

            env_key = (
                row.get("sample_type", ""),
                row.get("source_organism", ""),
                row.get("pmid", ""),
                row.get("evidence", ""),
            )

            if env_key not in environmental_seen:
                environmental_seen.add(env_key)

                source_organism_raw = row.get(
                    "source_organism",
                    ""
                )

                source_taxonomy = (
                    get_taxonomy_context_by_taxid(
                        source_organism_raw
                    )
                    if str(source_organism_raw).strip().isdigit()
                    else empty_taxonomy_context(
                        source_organism_raw
                    )
                )

                environmental_records.append({
                    "sample_type": row.get(
                        "sample_type",
                        ""
                    ),
                    "source_organism": source_organism_raw,
                    "source_organism_tax_id": source_taxonomy.get(
                        "tax_id"
                    ),
                    "source_organism_name": source_taxonomy.get(
                        "scientific_name"
                    ),
                    "source_organism_group": source_taxonomy.get(
                        "broad_group",
                        "UNKNOWN"
                    ),
                    "source_organism_family": source_taxonomy.get(
                        "by_rank",
                        {}
                    ).get(
                        "family"
                    ),
                    "source_organism_order": source_taxonomy.get(
                        "by_rank",
                        {}
                    ).get(
                        "order"
                    ),
                    "pmid": row.get(
                        "pmid",
                        ""
                    ),
                    "evidence": row.get(
                        "evidence",
                        ""
                    ),
                })

            continue

        key = (
            row["host_tax_id"]
            or normalize_name(
                row["host_name"]
            )
        )

        if not key:
            continue

        unique_hosts[key] = row

    target_tax_id = str(
        host_taxonomy.get("tax_id")
        or ""
    )

    target_name_norm = normalize_name(host)

    public_hosts = []
    direct_target_match = False

    for row in unique_hosts.values():
        is_target = (
            (
                target_tax_id
                and row["host_tax_id"] == target_tax_id
            )
            or normalize_name(row["host_name"]) == target_name_norm
        )

        if is_target:
            # Hidden diagnostic only. Never expose this to the policy.
            direct_target_match = True
            continue

        distance = host_taxonomic_distance(
            host_taxonomy,
            row["host_name"],
            row["host_lineage"],
        )

        public_hosts.append({
            "host_tax_id": row["host_tax_id"],
            "host_name": row["host_name"],
            "broad_group": broad_group_from_lineage(
                row["host_lineage"]
            ),
            "distance": distance["distance"],
            "lca_rank": distance["lca_rank"],
            "lca_name": distance["lca_name"],
            "evidence": row["evidence"],
        })

    public_hosts.sort(
        key=lambda item: (
            item["distance"],
            item["host_name"],
        )
    )

    group_counts = Counter(
        item["broad_group"]
        for item in public_hosts
        if item["broad_group"] != "UNKNOWN"
    )

    target_group = host_taxonomy.get(
        "broad_group",
        "UNKNOWN",
    )

    if not public_hosts:
        prior = "UNKNOWN"

        if environmental_records:
            prior_reason = (
                "Only environmental or unresolved Virus-Host DB "
                "records were available. These records are not treated "
                "as biological hosts and cannot establish host-range "
                "compatibility."
            )
        else:
            prior_reason = (
                "No non-target host-range records were available for "
                "the exact virus or its taxonomic fallback."
            )

    elif target_group == "UNKNOWN":
        prior = "UNKNOWN"
        prior_reason = (
            "The target host broad biological group could not be "
            "resolved reliably."
        )

    elif target_group in group_counts:
        prior = "BIOLOGICALLY_CONSISTENT"
        prior_reason = (
            "The target host belongs to a broad biological group "
            "already represented among known non-target hosts."
        )

    else:
        prior = "BIOLOGICALLY_DISCORDANT"
        prior_reason = (
            "The target host belongs to a different broad biological "
            "group from the currently observed non-target host range. "
            "This is only a prior and must not override explicit "
            "literature evidence."
        )

    context = {
        "target_host": {
            "scientific_name": host_taxonomy.get("scientific_name"),
            "tax_id": host_taxonomy.get("tax_id"),
            "broad_group": target_group,
            "family": host_taxonomy.get("by_rank", {}).get("family"),
            "order": host_taxonomy.get("by_rank", {}).get("order"),
            "class": host_taxonomy.get("by_rank", {}).get("class"),
        },
        "target_virus": {
            "scientific_name": virus_taxonomy.get("scientific_name"),
            "tax_id": virus_taxonomy.get("tax_id"),
            "family": virus_taxonomy.get("by_rank", {}).get("family"),
            "genus": virus_taxonomy.get("by_rank", {}).get("genus"),
        },
        "host_range": {
            "match_level": match_level,
            "known_non_target_hosts": len(public_hosts),
            "broad_group_counts": dict(group_counts),
            "closest_known_hosts": public_hosts[:8],
            "environmental_record_count": len(
                environmental_records
            ),
            "environmental_contexts": environmental_records[:8],
        },
        "biological_prior": {
            "status": prior,
            "reason": prior_reason,
        },
        "_diagnostics": {
            "target_was_present_in_database": direct_target_match,
            "raw_host_records": len(unique_hosts),
            "environmental_records": len(environmental_records),
            "virus_host_db_path": VIRUS_HOST_DB,
        },
    }

    cache[cache_key] = context
    save_cache(cache)

    return context


# ============================================================
# SAFE PROMPT VIEW
# ============================================================


def context_for_prompt(context):
    """Remove hidden diagnostic fields before passing context to an LLM."""

    safe_context = {
        key: value
        for key, value in (context or {}).items()
        if not key.startswith("_")
    }

    return json.dumps(
        safe_context,
        indent=2,
        ensure_ascii=False,
    )
