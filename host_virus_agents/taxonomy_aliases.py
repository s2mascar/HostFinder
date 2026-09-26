import json
import os
import re
import time
import xml.etree.ElementTree as ET

import requests

from literature_search import (
    NCBI_BASE,
    EMAIL,
    TOOL,
    TIMEOUT
)


CACHE_FILE = os.environ.get("HOST_ALIAS_CACHE", "host_alias_cache_correctness.json")
MAX_RETRIES = 4


# ============================================================
# NORMALIZATION
# ============================================================

def normalize_name(text):

    if not text:
        return ""

    text = str(text).lower()

    text = re.sub(
        r"[^a-z0-9]+",
        " ",
        text
    )

    return " ".join(
        text.split()
    )


# ============================================================
# CACHE
# ============================================================

def load_cache():

    if not os.path.exists(
        CACHE_FILE
    ):
        return {}

    try:

        with open(
            CACHE_FILE,
            "r",
            encoding="utf-8"
        ) as f:

            return json.load(f)

    except Exception:

        return {}


def save_cache(cache):

    with open(
        CACHE_FILE,
        "w",
        encoding="utf-8"
    ) as f:

        json.dump(
            cache,
            f,
            indent=2,
            ensure_ascii=False
        )


# ============================================================
# REQUEST WITH RETRIES
# ============================================================

def ncbi_get(
    endpoint,
    params
):

    last_error = None

    for attempt in range(
        MAX_RETRIES
    ):

        try:

            response = requests.get(
                f"{NCBI_BASE}/{endpoint}",
                params=params,
                timeout=TIMEOUT
            )

            if response.status_code == 429:

                wait = 2 ** attempt

                print(
                    "NCBI Taxonomy rate limited. "
                    f"Retrying in {wait}s..."
                )

                time.sleep(wait)

                continue

            if response.status_code in {
                500,
                502,
                503,
                504
            }:

                wait = 2 ** attempt

                time.sleep(wait)

                continue

            response.raise_for_status()

            # Stay comfortably below unauthenticated
            # NCBI request limits.
            time.sleep(0.4)

            return response

        except requests.RequestException as error:

            last_error = error

            if attempt == (
                MAX_RETRIES - 1
            ):
                break

            time.sleep(
                2 ** attempt
            )

    raise RuntimeError(
        f"NCBI Taxonomy request failed: "
        f"{last_error}"
    )


# ============================================================
# GENUS ABBREVIATION
# ============================================================

def scientific_abbreviation(
    scientific_name
):

    parts = scientific_name.split()

    if len(parts) < 2:
        return None

    genus = parts[0]
    rest = " ".join(
        parts[1:]
    )

    if not genus:
        return None

    return (
        f"{genus[0]}. {rest}"
    )


# ============================================================
# TAXONOMY LOOKUP
# ============================================================

def fetch_taxonomy_aliases(
    host
):

    # --------------------------------------------------------
    # Search taxonomy
    # --------------------------------------------------------

    search_params = {
        "db": "taxonomy",
        "term": f'"{host}"',
        "retmode": "json",
        "retmax": 5,
        "email": EMAIL,
        "tool": TOOL
    }

    response = ncbi_get(
        "esearch.fcgi",
        search_params
    )

    ids = (
        response
        .json()
        .get(
            "esearchresult",
            {}
        )
        .get(
            "idlist",
            []
        )
    )

    if not ids:
        return [
            host
        ]

    # --------------------------------------------------------
    # Fetch possible taxonomy records
    # --------------------------------------------------------

    fetch_params = {
        "db": "taxonomy",
        "id": ",".join(ids),
        "retmode": "xml",
        "email": EMAIL,
        "tool": TOOL
    }

    response = ncbi_get(
        "efetch.fcgi",
        fetch_params
    )

    root = ET.fromstring(
        response.text
    )

    target_norm = normalize_name(
        host
    )

    chosen_taxon = None

    # Prefer an exact scientific-name match.
    for taxon in root.findall(
        "./Taxon"
    ):

        scientific_name = (
            taxon.findtext(
                "ScientificName"
            )
            or ""
        )

        if (
            target_norm in {normalize_name(scientific_name), *(
                normalize_name(node.text or "")
                for tag in ("Synonym", "EquivalentName", "GenbankSynonym", "CommonName", "GenbankCommonName")
                for node in taxon.findall("OtherNames/" + tag)
            )}
        ):

            chosen_taxon = taxon
            break

    # An unmatched search hit is not an alias resolution.
    if chosen_taxon is None:
        return [host]

    if chosen_taxon is None:

        return [
            host
        ]

    aliases = set()

    # --------------------------------------------------------
    # Scientific name
    # --------------------------------------------------------

    scientific_name = (
        chosen_taxon.findtext(
            "ScientificName"
        )
        or host
    )

    aliases.add(
        scientific_name
    )

    aliases.add(
        host
    )

    abbreviation = (
        scientific_abbreviation(
            scientific_name
        )
    )

    if abbreviation:
        aliases.add(
            abbreviation
        )

    # --------------------------------------------------------
    # NCBI synonyms / common names
    # --------------------------------------------------------

    other_names = chosen_taxon.find(
        "OtherNames"
    )

    if other_names is not None:

        accepted_tags = {
            "Synonym",
            "EquivalentName",
            "GenbankCommonName",
            "CommonName",
            "GenbankSynonym"
        }

        for child in other_names:

            tag = child.tag.split(
                "}"
            )[-1]

            if (
                tag in accepted_tags
                and child.text
            ):

                name = (
                    child.text
                    .strip()
                )

                if name:
                    aliases.add(
                        name
                    )

    # --------------------------------------------------------
    # Clean
    # --------------------------------------------------------

    aliases = [
        alias
        for alias in aliases
        if alias
        and len(alias) >= 2
    ]

    aliases = sorted(
        aliases,
        key=lambda x: (
            -len(x),
            x.lower()
        )
    )

    return aliases


# ============================================================
# PUBLIC FUNCTION
# ============================================================

def get_host_aliases(
    host
):

    cache = load_cache()

    cache_key = normalize_name(
        host
    )

    if cache_key in cache:

        aliases = cache[
            cache_key
        ]

        print(
            f"Host aliases loaded from cache: "
            f"{aliases}"
        )

        return aliases

    try:

        aliases = (
            fetch_taxonomy_aliases(
                host
            )
        )

    except Exception as error:

        print(
            "WARNING: Could not retrieve "
            f"NCBI taxonomy aliases: {error}"
        )

        aliases = [
            host
        ]

    cache[
        cache_key
    ] = aliases

    save_cache(
        cache
    )

    print(
        f"Host aliases: {aliases}"
    )

    return aliases


# ============================================================
# ALIAS MATCHING
# ============================================================

def text_contains_alias(
    text,
    aliases
):

    if not text:
        return False

    text_norm = normalize_name(
        text
    )

    for alias in aliases:

        alias_norm = normalize_name(
            alias
        )

        if not alias_norm:
            continue

        pattern = (
            r"(?<!\w)"
            + re.escape(
                alias_norm
            )
            + r"(?!\w)"
        )

        if re.search(
            pattern,
            text_norm
        ):

            return True

    return False


def name_matches_host(
    reported_name,
    aliases
):

    if not reported_name:
        return False

    reported_norm = (
        normalize_name(
            reported_name
        )
    )

    for alias in aliases:

        if (
            reported_norm
            == normalize_name(
                alias
            )
        ):
            return True

    return False
