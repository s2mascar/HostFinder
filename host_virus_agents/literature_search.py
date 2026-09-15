import requests
import xml.etree.ElementTree as ET
import sys

# ============================================================
# SETTINGS
# ============================================================

NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

EUROPE_PMC_BASE = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
)

EMAIL = "s2mascar@uwaterloo.ca"
TOOL = "hostfinder_literature_agent"

TIMEOUT = 30


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def clean_text(text):
    """Remove extra whitespace."""

    if not text:
        return ""

    return " ".join(text.split())


def element_text(element):
    """Extract all text from an XML element."""

    if element is None:
        return ""

    return clean_text(
        "".join(element.itertext())
    )


def normalize_pmcid(pmcid):
    """Convert 1234567 -> PMC1234567."""

    if not pmcid:
        return None

    pmcid = str(pmcid).strip()

    if pmcid.upper().startswith("PMC"):
        return pmcid.upper()

    return f"PMC{pmcid}"


# ============================================================
# PUBMED SEARCH
# ============================================================

def search_pubmed(host, virus, max_results=20):

    query = f'"{host}" AND "{virus}"'

    print("\n" + "=" * 80)
    print("SEARCHING PUBMED")
    print("=" * 80)
    print(query)

    # --------------------------------------------------------
    # Find matching PubMed IDs
    # --------------------------------------------------------

    search_params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": max_results,
        "email": EMAIL,
        "tool": TOOL,
    }

    response = requests.get(
        f"{NCBI_BASE}/esearch.fcgi",
        params=search_params,
        timeout=TIMEOUT
    )

    response.raise_for_status()

    data = response.json()

    pmids = data["esearchresult"]["idlist"]

    print(f"Found {len(pmids)} PubMed records.")

    if not pmids:
        return []

    # --------------------------------------------------------
    # Download PubMed records
    # --------------------------------------------------------

    fetch_params = {
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
        "email": EMAIL,
        "tool": TOOL,
    }

    response = requests.get(
        f"{NCBI_BASE}/efetch.fcgi",
        params=fetch_params,
        timeout=TIMEOUT
    )

    response.raise_for_status()

    root = ET.fromstring(response.text)

    papers = []

    for article in root.findall(".//PubmedArticle"):

        pmid = article.findtext(".//PMID")

        # Title
        title_element = article.find(".//ArticleTitle")
        title = element_text(title_element)

        # Abstract
        abstract_parts = []

        for abstract_element in article.findall(
            ".//AbstractText"
        ):
            text = element_text(abstract_element)

            if text:
                abstract_parts.append(text)

        abstract = " ".join(abstract_parts)

        # IDs
        pmcid = None
        doi = None

        for article_id in article.findall(
            ".//ArticleId"
        ):

            id_type = article_id.attrib.get(
                "IdType"
            )

            if id_type == "pmc":
                pmcid = normalize_pmcid(
                    article_id.text
                )

            elif id_type == "doi":
                doi = article_id.text

        papers.append({
            "title": title,
            "abstract": abstract,
            "full_text": "",
            "pmid": pmid,
            "pmcid": pmcid,
            "doi": doi,
            "sources": ["PubMed"],
            "search_queries": []
        })

    return papers


# ============================================================
# PMC FULL-TEXT SEARCH
# ============================================================

def search_pmc(host, virus, max_results=20):

    query = f'"{host}" AND "{virus}"'

    print("\n" + "=" * 80)
    print("SEARCHING PMC FULL TEXT")
    print("=" * 80)
    print(query)

    # --------------------------------------------------------
    # Search PMC
    # --------------------------------------------------------

    search_params = {
        "db": "pmc",
        "term": query,
        "retmode": "json",
        "retmax": max_results,
        "email": EMAIL,
        "tool": TOOL,
    }

    response = requests.get(
        f"{NCBI_BASE}/esearch.fcgi",
        params=search_params,
        timeout=TIMEOUT
    )

    response.raise_for_status()

    data = response.json()

    pmc_ids = data["esearchresult"]["idlist"]

    print(f"Found {len(pmc_ids)} PMC records.")

    if not pmc_ids:
        return []

    # --------------------------------------------------------
    # Download PMC full text
    # --------------------------------------------------------

    fetch_params = {
        "db": "pmc",
        "id": ",".join(pmc_ids),
        "retmode": "xml",
        "email": EMAIL,
        "tool": TOOL,
    }

    response = requests.get(
        f"{NCBI_BASE}/efetch.fcgi",
        params=fetch_params,
        timeout=TIMEOUT
    )

    response.raise_for_status()

    root = ET.fromstring(response.text)

    papers = []

    for article in root.findall(".//article"):

        pmid = None
        pmcid = None
        doi = None

        # IDs
        for article_id in article.findall(
            ".//article-id"
        ):

            id_type = article_id.attrib.get(
                "pub-id-type"
            )

            if id_type == "pmid":
                pmid = article_id.text

            elif id_type == "pmc":
                pmcid = normalize_pmcid(
                    article_id.text
                )

            elif id_type == "doi":
                doi = article_id.text

        # Title
        title = element_text(
            article.find(".//article-title")
        )

        # Abstract
        abstract = element_text(
            article.find(".//abstract")
        )

        # Full text
        body = article.find(".//body")

        full_text = element_text(body)

        papers.append({
            "title": title,
            "abstract": abstract,
            "full_text": full_text,
            "pmid": pmid,
            "pmcid": pmcid,
            "doi": doi,
            "sources": ["PMC"],
            "search_queries": []
        })

    return papers


# ============================================================
# EUROPE PMC SEARCH
# ============================================================

def search_europe_pmc(
    host,
    virus,
    max_results=20
):

    query = f'"{host}" AND "{virus}"'

    print("\n" + "=" * 80)
    print("SEARCHING EUROPE PMC")
    print("=" * 80)
    print(query)

    params = {
        "query": query,
        "format": "json",
        "resultType": "core",
        "pageSize": max_results,
    }

    response = requests.get(
        EUROPE_PMC_BASE,
        params=params,
        timeout=TIMEOUT
    )

    response.raise_for_status()

    data = response.json()

    results = (
        data
        .get("resultList", {})
        .get("result", [])
    )

    print(
        f"Found {len(results)} Europe PMC records."
    )

    papers = []

    for result in results:

        pmid = result.get("pmid")

        pmcid = normalize_pmcid(
            result.get("pmcid")
        )

        doi = result.get("doi")

        title = clean_text(
            result.get("title", "")
        )

        abstract = clean_text(
            result.get("abstractText", "")
        )

        papers.append({
            "title": title,
            "abstract": abstract,
            "full_text": "",
            "pmid": pmid,
            "pmcid": pmcid,
            "doi": doi,
            "sources": ["Europe PMC"],
            "search_queries": []
        })

    return papers


# ============================================================
# DEDUPLICATION
# ============================================================

def paper_key(paper):
    """
    Generate an identifier used to determine whether
    two search results are the same paper.
    """

    if paper.get("pmid"):
        return f"pmid:{paper['pmid']}"

    if paper.get("pmcid"):
        return f"pmcid:{paper['pmcid']}"

    if paper.get("doi"):
        return f"doi:{paper['doi'].lower()}"

    title = paper.get(
        "title",
        ""
    ).strip().lower()

    return f"title:{title}"


def merge_papers(all_papers):

    merged = {}

    for paper in all_papers:

        key = paper_key(paper)

        if key not in merged:

            merged[key] = paper.copy()

            merged[key]["sources"] = list(
                paper.get(
                    "sources",
                    []
                )
            )

            merged[key]["search_queries"] = list(
                paper.get(
                    "search_queries",
                    []
                )
            )

            continue

        existing = merged[key]

        # ----------------------------------------------------
        # Merge sources
        # ----------------------------------------------------

        existing["sources"] = sorted(
            set(
                existing.get(
                    "sources",
                    []
                )
                +
                paper.get(
                    "sources",
                    []
                )
            )
        )

        # ----------------------------------------------------
        # Merge search queries
        # ----------------------------------------------------

        for search_query in paper.get(
            "search_queries",
            []
        ):

            if search_query not in existing[
                "search_queries"
            ]:

                existing[
                    "search_queries"
                ].append(
                    search_query
                )

        # ----------------------------------------------------
        # Keep longest abstract
        # ----------------------------------------------------

        if (
            len(
                paper.get(
                    "abstract",
                    ""
                )
            )
            >
            len(
                existing.get(
                    "abstract",
                    ""
                )
            )
        ):

            existing["abstract"] = (
                paper["abstract"]
            )

        # ----------------------------------------------------
        # Keep longest full text
        # ----------------------------------------------------

        if (
            len(
                paper.get(
                    "full_text",
                    ""
                )
            )
            >
            len(
                existing.get(
                    "full_text",
                    ""
                )
            )
        ):

            existing["full_text"] = (
                paper["full_text"]
            )

        # ----------------------------------------------------
        # Fill missing IDs
        # ----------------------------------------------------

        for field in [
            "pmid",
            "pmcid",
            "doi"
        ]:

            if not existing.get(field):

                existing[field] = (
                    paper.get(field)
                )

    return list(
        merged.values()
    )


# ============================================================
# SEARCH ALL SOURCES
# ============================================================

def search_literature(
    host,
    virus,
    max_results=20
):

    pubmed_results = search_pubmed(
        host,
        virus,
        max_results
    )

    pmc_results = search_pmc(
        host,
        virus,
        max_results
    )

    europe_pmc_results = search_europe_pmc(
        host,
        virus,
        max_results
    )

    all_results = (
        pubmed_results
        +
        pmc_results
        +
        europe_pmc_results
    )

    merged_results = merge_papers(
        all_results
    )

    print("\n" + "=" * 80)
    print("SEARCH SUMMARY")
    print("=" * 80)

    print(
        f"PubMed:                     "
        f"{len(pubmed_results)}"
    )

    print(
        f"PMC full text:               "
        f"{len(pmc_results)}"
    )

    print(
        f"Europe PMC:                  "
        f"{len(europe_pmc_results)}"
    )

    print(
        f"Unique papers after merging: "
        f"{len(merged_results)}"
    )

    return merged_results


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    if len(sys.argv) != 3:

        print(
            'Usage:\n'
            'python literature_search.py '
            '"HOST" "VIRUS"'
        )

        sys.exit(1)

    host = sys.argv[1]
    virus = sys.argv[2]

    papers = search_literature(
        host,
        virus
    )

    for i, paper in enumerate(
        papers,
        start=1
    ):

        print("\n" + "=" * 80)
        print(f"PAPER {i}")
        print("=" * 80)

        print(
            "Sources:",
            ", ".join(
                paper.get(
                    "sources",
                    []
                )
            )
        )

        print(
            "PMID:",
            paper.get("pmid")
        )

        print(
            "PMCID:",
            paper.get("pmcid")
        )

        print(
            "DOI:",
            paper.get("doi")
        )

        print("\nTITLE:")

        print(
            paper.get(
                "title",
                ""
            )
        )

        print("\nABSTRACT:")

        abstract = paper.get(
            "abstract",
            ""
        )

        if abstract:
            print(
                abstract[:1500]
            )
        else:
            print(
                "No abstract available."
            )

        if paper.get(
            "full_text"
        ):

            print("\nFULL TEXT PREVIEW:")

            print(
                paper[
                    "full_text"
                ][:2000]
            )