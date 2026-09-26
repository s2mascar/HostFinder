import json
import re
import sys

from evidence_semantics import decide, retrieval_sufficient


def run_evidence_agent(*args, **kwargs):
    from evidence_agent import run_evidence_agent as run
    return run(*args, **kwargs)


# ============================================================
# RETRIEVAL SUFFICIENCY
# ============================================================

def determine_final_status(evidence_results, search_metadata):
    return decide(evidence_results, search_metadata)["literature_status"]


# ============================================================
# DETERMINISTIC REASON
# ============================================================

def build_reason(
    host,
    virus,
    evidence_results,
    final_status,
    search_metadata
):

    if final_status == "KNOWN":

        exact = [
            result
            for result
            in evidence_results
            if result.get(
                "classification"
            )
            == "EXACT_SUPPORT"
        ]

        return (
            f"At least one retrieved paper "
            f"explicitly associates {virus} "
            f"with {host}. "
            f"{len(exact)} paper(s) provided "
            f"verified exact support."
        )

    if (
        final_status
        == "POSSIBLY_KNOWN"
    ):

        related = [
            result
            for result
            in evidence_results
            if result.get(
                "classification"
            )
            == "TARGET_HOST_RELATED"
        ]

        strongest = (
            related[0]
            if related
            else {}
        )

        reported_virus = (
            strongest.get(
                "reported_virus_name"
            )
            or "a related virus"
        )

        return (
            f"No retrieved paper provides "
            f"exact support for {host} and "
            f"{virus}. However, literature "
            f"from the target host reports "
            f"{reported_virus} in a directly "
            f"related viral context."
        )

    if (
        final_status
        == "NO_EVIDENCE_FOUND"
    ):

        other_host_count = sum(
            1
            for result
            in evidence_results
            if result.get(
                "classification"
            )
            == "VIRUS_OTHER_HOST"
        )

        if other_host_count:

            return (
                f"No retrieved paper provides "
                f"evidence associating {virus} "
                f"with {host}. "
                f"{other_host_count} retrieved "
                f"paper(s) discussed the target "
                f"virus in other hosts or samples, "
                f"which does not support the target "
                f"interaction."
            )

        return (
            f"No retrieved paper provided "
            f"exact or target-host-related "
            f"evidence for {host} and {virus}. "
            f"This means no evidence was found "
            f"by the current literature pipeline; "
            f"it does not prove the biological "
            f"interaction is absent."
        )

    return (
        "A verified passage contained the "
        "target host and virus, but the "
        "relationship could not be resolved "
        "with sufficient confidence."
    )


# ============================================================
# STRONGEST EVIDENCE SUMMARY
# ============================================================

def build_evidence_summary(
    evidence_results
):

    sections = []

    for i, result in enumerate(
        evidence_results,
        start=1
    ):

        sections.append(
            f"""
PAPER {i}

Title:
{result.get("title")}

Evidence class:
{result.get("classification")}

Relationship type:
{result.get("relationship_type")}

Reported virus:
{result.get("reported_virus_name")}

Reported host:
{result.get("reported_host_or_sample")}

Verified passage:
{result.get("relationship_passage")}
"""
        )

    return "\n".join(
        sections
    )


# ============================================================
# CONFIDENCE
# ============================================================

def default_confidence(
    final_status,
    search_metadata
):

    if final_status == "UNCLEAR":
        return "LOW"

    if not retrieval_sufficient(
        search_metadata
    ):
        return "LOW"

    if final_status == "KNOWN":
        return "HIGH"

    if (
        final_status
        == "POSSIBLY_KNOWN"
    ):
        return "MEDIUM"

    return "MEDIUM"


# ============================================================
# JUDGE
# ============================================================

def judge_interaction(
    host,
    virus,
    evidence_results,
    search_metadata
):

    decision = decide(evidence_results, search_metadata)
    final_status = decision["literature_status"]
    reason = decision["reason"]

    exact_supporting_papers = []

    related_papers = []

    other_host_papers = []

    for result in evidence_results:

        paper = {
            "title":
                result.get(
                    "title"
                ),

            "pmid":
                result.get(
                    "pmid"
                ),

            "pmcid":
                result.get(
                    "pmcid"
                ),

            "doi":
                result.get(
                    "doi"
                )
        }

        classification = (
            result.get(
                "classification"
            )
        )

        if (
            classification
            == "EXACT_SUPPORT"
        ):

            exact_supporting_papers.append(
                paper
            )

        elif (
            classification
            == "TARGET_HOST_RELATED"
        ):

            related_papers.append(
                paper
            )

        elif (
            classification
            == "VIRUS_OTHER_HOST"
        ):

            other_host_papers.append(
                paper
            )

    confidence = decision["confidence"]
    strongest_evidence = next((r["structured_evidence"]["supporting_text"]
                               for r in evidence_results
                               if r.get("structured_evidence", {}).get("verified")), "")

    return {
        "classification": decision["classification"],
        "evidence_results": evidence_results,
        "search_metadata": search_metadata,
        "host":
            host,

        "virus":
            virus,

        "literature_status":
            final_status,

        "papers_examined":
            len(
                evidence_results
            ),

        "exact_supporting_papers":
            exact_supporting_papers,

        "related_papers":
            related_papers,

        "other_host_papers":
            other_host_papers,

        "confidence":
            confidence,

        "reason":
            reason,

        "strongest_evidence":
            strongest_evidence,

        "retrieval_complete":
            search_metadata.get(
                "retrieval_complete",
                False
            ),

        "search_source_failures":
            search_metadata.get(
                "source_failures",
                []
            )
    }


# ============================================================
# COMPLETE PIPELINE
# ============================================================

def run_judge_agent(
    host,
    virus
):

    output = run_evidence_agent(
        host,
        virus
    )

    return judge_interaction(
        host,
        virus,
        output[
            "evidence_results"
        ],
        output[
            "search_metadata"
        ]
    )


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    if len(sys.argv) != 3:

        print(
            'Usage:\n'
            'python judge_agent.py '
            '"HOST" "VIRUS"'
        )

        sys.exit(1)

    result = run_judge_agent(
        sys.argv[1],
        sys.argv[2]
    )

    print("\n" + "=" * 80)
    print("FINAL LITERATURE VERDICT")
    print("=" * 80)

    print(
        "Host:",
        result["host"]
    )

    print(
        "Virus:",
        result["virus"]
    )

    print(
        "\nLiterature status:",
        result[
            "literature_status"
        ]
    )

    print(
        "Papers examined:",
        result[
            "papers_examined"
        ]
    )

    print(
        "Confidence:",
        result[
            "confidence"
        ]
    )

    print("\nReason:")

    print(
        result[
            "reason"
        ]
    )

    print(
        "\nExact supporting papers:",
        len(
            result[
                "exact_supporting_papers"
            ]
        )
    )

    print(
        "Target-host-related papers:",
        len(
            result[
                "related_papers"
            ]
        )
    )

    print(
        "Target-virus other-host papers:",
        len(
            result[
                "other_host_papers"
            ]
        )
    )