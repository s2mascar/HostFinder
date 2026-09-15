import json
import re
import sys

from search_agent import (
    run_search_agent,
    generate_text,
    cleanup_gpu
)

from taxonomy_aliases import (
    get_host_aliases,
    text_contains_alias,
    name_matches_host
)


# ============================================================
# SETTINGS
# ============================================================

SNIPPET_WINDOW = 900

MAX_HOST_SNIPPETS = 4
MAX_VIRUS_SNIPPETS = 6

MAX_OCCURRENCES_PER_TERM = 25

MAX_ABSTRACT_CHARS = 7000

MAX_INPUT_TOKENS = 5500


# ============================================================
# NORMALIZATION
# ============================================================

def normalize_entity(
    text
):

    if not text:
        return ""

    text = str(
        text
    ).lower()

    text = re.sub(
        r"[^a-z0-9]+",
        " ",
        text
    )

    return " ".join(
        text.split()
    )


def normalize_whitespace(
    text
):

    if not text:
        return ""

    return " ".join(
        str(
            text
        ).split()
    )


def same_virus(
    reported,
    target
):

    reported_norm = (
        normalize_entity(
            reported
        )
    )

    target_norm = (
        normalize_entity(
            target
        )
    )

    if (
        not reported_norm
        or not target_norm
    ):
        return False

    if (
        reported_norm
        == target_norm
    ):
        return True

    # Allow extended names such as:
    #
    # SARS-CoV-2
    # SARS-CoV-2 (severe acute ...)
    if (
        len(target_norm) >= 8
        and target_norm
        in reported_norm
    ):
        return True

    return False


def text_contains_name(
    text,
    name
):

    if (
        not text
        or not name
    ):
        return False

    return (
        normalize_entity(
            name
        )
        in normalize_entity(
            text
        )
    )


# ============================================================
# PASSAGE VERIFICATION
# ============================================================

def verify_passage(
    passage,
    source_text
):

    if not passage:
        return False

    passage_norm = (
        normalize_whitespace(
            passage
        )
    )

    source_norm = (
        normalize_whitespace(
            source_text
        )
    )

    return (
        passage_norm
        in source_norm
    )


# ============================================================
# SNIPPET RANKING
# ============================================================

def rank_snippets(
    candidates,
    max_snippets
):

    unique = {}

    for score, snippet in candidates:

        cleaned = (
            normalize_whitespace(
                snippet
            )
        )

        if not cleaned:
            continue

        if (
            cleaned not in unique
            or score
            > unique[cleaned]
        ):

            unique[
                cleaned
            ] = score

    ranked = sorted(
        unique.items(),
        key=lambda item:
            item[1],
        reverse=True
    )

    selected = [
        snippet
        for snippet, score
        in ranked[
            :max_snippets
        ]
    ]

    return (
        "\n\n---\n\n".join(
            selected
        )
    )


# ============================================================
# HOST-CONTEXT SNIPPETS
# ============================================================

def get_host_context_snippets(
    text,
    host_aliases,
    window=SNIPPET_WINDOW,
    max_snippets=MAX_HOST_SNIPPETS
):

    if not text:
        return ""

    text_lower = (
        text.lower()
    )

    candidates = []

    for alias in host_aliases:

        if not alias:
            continue

        alias_lower = (
            alias.lower()
        )

        start = 0
        count = 0

        while (
            count
            < MAX_OCCURRENCES_PER_TERM
        ):

            position = (
                text_lower.find(
                    alias_lower,
                    start
                )
            )

            if position == -1:
                break

            left = max(
                0,
                position - window
            )

            right = min(
                len(text),
                position
                + len(alias)
                + window
            )

            snippet = (
                text[
                    left:right
                ]
            )

            score = 5

            lower = (
                snippet.lower()
            )

            for word in [
                "virus",
                "viral",
                "sequenc",
                "detected",
                "identified",
                "sample",
                "collected",
                "isolated",
                "transcriptome",
                "rna"
            ]:

                if word in lower:
                    score += 1

            candidates.append(
                (
                    score,
                    snippet
                )
            )

            start = (
                position
                + len(alias)
            )

            count += 1

    return rank_snippets(
        candidates,
        max_snippets
    )


# ============================================================
# TARGET-VIRUS CONTEXT
# ============================================================

def get_virus_context_snippets(
    text,
    host_aliases,
    virus,
    window=SNIPPET_WINDOW,
    max_snippets=MAX_VIRUS_SNIPPETS
):

    if not text:
        return ""

    text_lower = (
        text.lower()
    )

    virus_lower = (
        virus.lower()
    )

    candidates = []

    start = 0
    count = 0

    while (
        count
        < MAX_OCCURRENCES_PER_TERM
    ):

        position = (
            text_lower.find(
                virus_lower,
                start
            )
        )

        if position == -1:
            break

        left = max(
            0,
            position - window
        )

        right = min(
            len(text),
            position
            + len(virus)
            + window
        )

        snippet = (
            text[
                left:right
            ]
        )

        score = 7

        if text_contains_alias(
            snippet,
            host_aliases
        ):
            score += 12

        lower = (
            snippet.lower()
        )

        for word in [
            "similar",
            "identity",
            "phylogen",
            "closest",
            "related",
            "detected",
            "isolated",
            "identified",
            "sequenced",
            "infection",
            "infected"
        ]:

            if word in lower:
                score += 1

        candidates.append(
            (
                score,
                snippet
            )
        )

        start = (
            position
            + len(virus)
        )

        count += 1

    return rank_snippets(
        candidates,
        max_snippets
    )


# ============================================================
# PYTHON EVIDENCE CLASSIFICATION
# ============================================================

def classify_relationship(
    extraction,
    host_aliases,
    target_virus,
    evidence_source
):

    study_host = (
        extraction.get(
            "study_host",
            ""
        )
        or ""
    )

    study_host_passage = (
        extraction.get(
            "study_host_passage",
            ""
        )
        or ""
    )

    host_virus_name = (
        extraction.get(
            "host_virus_name",
            ""
        )
        or ""
    )

    host_virus_passage = (
        extraction.get(
            "host_virus_passage",
            ""
        )
        or ""
    )

    comparison_virus_name = (
        extraction.get(
            "comparison_virus_name",
            ""
        )
        or ""
    )

    comparison_host_or_sample = (
        extraction.get(
            "comparison_host_or_sample",
            ""
        )
        or ""
    )

    relationship_type = (
        extraction.get(
            "relationship_type",
            "UNCLEAR"
        )
        .upper()
        .strip()
    )

    relationship_passage = (
        extraction.get(
            "relationship_passage",
            ""
        )
        or ""
    )

    # ========================================================
    # VERIFY TEXT
    # ========================================================

    study_host_passage_verified = (
        verify_passage(
            study_host_passage,
            evidence_source
        )
    )

    host_virus_passage_verified = (
        verify_passage(
            host_virus_passage,
            evidence_source
        )
    )

    relationship_passage_verified = (
        verify_passage(
            relationship_passage,
            evidence_source
        )
    )

    # ========================================================
    # RESOLVE HOST
    # ========================================================

    study_host_matches_target = (
        name_matches_host(
            study_host,
            host_aliases
        )
        or (
            study_host_passage_verified
            and text_contains_alias(
                study_host_passage,
                host_aliases
            )
        )
    )

    # This replaces the old strong_study_host_context rule.
    #
    # If Qwen extracted the correct study host AND supplied
    # a genuine passage from the paper establishing study
    # context, that is enough.
    verified_target_host_context = (
        study_host_matches_target
        and study_host_passage_verified
    )

    # ========================================================
    # RESOLVE VIRUSES
    # ========================================================

    host_virus_matches_target = (
        same_virus(
            host_virus_name,
            target_virus
        )
    )

    comparison_virus_matches_target = (
        same_virus(
            comparison_virus_name,
            target_virus
        )
    )

    host_virus_passage_mentions_host_virus = (
        text_contains_name(
            host_virus_passage,
            host_virus_name
        )
        if host_virus_name
        else False
    )

    relationship_mentions_target = (
        text_contains_name(
            relationship_passage,
            target_virus
        )
    )

    # ========================================================
    # SAVE DIAGNOSTICS
    # ========================================================

    extraction[
        "study_host_passage_verified"
    ] = study_host_passage_verified

    extraction[
        "study_host_matches_target"
    ] = study_host_matches_target

    extraction[
        "verified_target_host_context"
    ] = verified_target_host_context

    # Backwards-compatible diagnostic field.
    extraction[
        "strong_study_host_context"
    ] = verified_target_host_context

    extraction[
        "host_virus_passage_verified"
    ] = host_virus_passage_verified

    extraction[
        "host_virus_matches_target"
    ] = host_virus_matches_target

    extraction[
        "comparison_virus_matches_target"
    ] = comparison_virus_matches_target

    extraction[
        "relationship_passage_verified"
    ] = relationship_passage_verified

    extraction[
        "relationship_mentions_target_virus"
    ] = relationship_mentions_target

    # Backwards compatibility with the existing Judge
    # and your diagnostic script.
    extraction[
        "reported_virus_name"
    ] = host_virus_name

    extraction[
        "reported_host_or_sample"
    ] = study_host

    extraction[
        "reported_virus_matches_target"
    ] = host_virus_matches_target

    extraction[
        "reported_host_matches_target"
    ] = study_host_matches_target

    extraction[
        "reported_virus_passage_verified"
    ] = host_virus_passage_verified

    # ========================================================
    # DIRECT / RELATED TYPES
    # ========================================================

    direct_types = {
        "DETECTED_IN",
        "ISOLATED_FROM",
        "INFECTION_OF",
        "EXPERIMENTAL_INFECTION",
        "SEQUENCED_FROM",
        "IDENTIFIED_IN"
    }

    related_types = {
        "SEQUENCE_SIMILARITY",
        "PHYLOGENETIC_COMPARISON",
        "RELATED_VIRUS"
    }

    # ========================================================
    # EXACT SUPPORT
    #
    # The virus actually recovered from the study host
    # is the target virus.
    # ========================================================

    if (
        verified_target_host_context
        and host_virus_matches_target
        and host_virus_passage_verified
        and host_virus_passage_mentions_host_virus
    ):

        return "EXACT_SUPPORT"

    # ========================================================
    # TARGET HOST RELATED
    #
    # The target host contains virus A.
    # Virus A is explicitly related to target virus B.
    # ========================================================

    if (
        verified_target_host_context
        and host_virus_name
        and not host_virus_matches_target
        and host_virus_passage_verified
        and host_virus_passage_mentions_host_virus
        and comparison_virus_matches_target
        and relationship_type
            in related_types
        and relationship_passage_verified
        and relationship_mentions_target
    ):

        return "TARGET_HOST_RELATED"

    # ========================================================
    # TARGET VIRUS BELONGS TO ANOTHER HOST
    # ========================================================

    if (
        comparison_virus_matches_target
        and comparison_host_or_sample
        and not (
            verified_target_host_context
            and host_virus_name
            and relationship_type
                in related_types
            and relationship_passage_verified
        )
    ):

        return "VIRUS_OTHER_HOST"

    # If the paper's study host is not our target,
    # but the virus found there is the exact target virus.
    if (
        not verified_target_host_context
        and host_virus_matches_target
        and host_virus_passage_verified
    ):

        return "VIRUS_OTHER_HOST"

    # ========================================================
    # MODEL SAYS RELATED BUT REQUIRED LINK IS INCOMPLETE
    # ========================================================

    if (
        relationship_type
        in related_types
    ):

        # If the target virus is actually present in the
        # relationship passage but we cannot establish the
        # target host -> host virus -> target virus chain,
        # treat as no support rather than a false positive.
        return "NO_SUPPORT"

    # ========================================================
    # EXPLICIT NON-SUPPORT
    # ========================================================

    if relationship_type == "MENTION_ONLY":
        return "MENTION_ONLY"

    if relationship_type in {
        "NO_RELATION",
        "DIFFERENT_HOST"
    }:

        if (
            host_virus_matches_target
            or comparison_virus_matches_target
        ):
            return "VIRUS_OTHER_HOST"

        return "NO_SUPPORT"

    # ========================================================
    # GENUINE UNCERTAINTY
    # ========================================================

    if relationship_type == "UNCLEAR":

        if (
            verified_target_host_context
            and (
                host_virus_matches_target
                or comparison_virus_matches_target
            )
        ):

            return "UNCLEAR"

        return "NO_SUPPORT"

    return "NO_SUPPORT"


# ============================================================
# ANALYZE ONE PAPER
# ============================================================

def analyze_paper(
    host,
    host_aliases,
    virus,
    paper
):

    full_text = (
        paper.get(
            "full_text",
            ""
        )
        or ""
    )

    abstract = (
        paper.get(
            "abstract",
            ""
        )
        or ""
    )

    abstract = (
        abstract[
            :MAX_ABSTRACT_CHARS
        ]
    )

    host_context = (
        get_host_context_snippets(
            full_text,
            host_aliases
        )
    )

    virus_context = (
        get_virus_context_snippets(
            full_text,
            host_aliases,
            virus
        )
    )

    evidence_source = f"""
TITLE:
{paper.get("title", "")}

ABSTRACT:
{abstract}

HOST/STUDY CONTEXT:
{host_context}

TARGET-VIRUS CONTEXT:
{virus_context}
"""

    aliases_text = "\n".join(
        f"- {alias}"
        for alias
        in host_aliases
    )

    # ========================================================
    # EXTRACTION PROMPT
    # ========================================================

    system_prompt = """
You are a scientific evidence extraction agent.

You are NOT responsible for the final classification.

Your job is to identify:

1. What organism/sample the paper is studying.
2. What virus was actually identified, detected, sequenced,
   isolated, or otherwise associated with that study host.
3. Whether that HOST-ASSOCIATED VIRUS is compared with another
   virus.
4. What comparison virus is used.

IMPORTANT DISTINCTION:

HOST VIRUS
=
the virus actually recovered from, identified in, detected in,
sequenced from, or associated with the study host/sample.

COMPARISON VIRUS
=
a different virus used for sequence similarity, phylogenetic
comparison, closest-match comparison, taxonomy, or reference.

NEVER put a comparison virus into host_virus_name merely because
it is the TARGET VIRUS.

Example:

TARGET HOST:
Lampyris noctiluca

TARGET VIRUS:
Hubei partiti-like virus 31

Paper says:

"Lampyris noctiluca partiti-like virus 1 (LnoPLV1) ..."

and later:

"the protein is similar to RdRP of Hubei partiti-like virus 31."

Correct extraction:

study_host:
Lampyris noctiluca

host_virus_name:
Lampyris noctiluca partiti-like virus 1

comparison_virus_name:
Hubei partiti-like virus 31

relationship_type:
SEQUENCE_SIMILARITY


Second example:

TARGET HOST:
Lampyris noctiluca

TARGET VIRUS:
Lampyris noctiluca partiti-like virus 1

The same paper should still produce:

host_virus_name:
Lampyris noctiluca partiti-like virus 1

comparison_virus_name:
Hubei partiti-like virus 31

The TARGET does not determine which field the virus belongs in.
The paper does.


Choose one relationship_type:

DETECTED_IN
IDENTIFIED_IN
ISOLATED_FROM
SEQUENCED_FROM
INFECTION_OF
EXPERIMENTAL_INFECTION
SEQUENCE_SIMILARITY
PHYLOGENETIC_COMPARISON
RELATED_VIRUS
DIFFERENT_HOST
MENTION_ONLY
NO_RELATION
UNCLEAR


PASSAGE RULES:

study_host_passage:
A verbatim passage establishing the study organism/sample.

host_virus_passage:
A verbatim passage naming the virus associated with that
study organism/sample.

relationship_passage:
A verbatim passage showing the relationship between the
host-associated virus and the comparison virus.

If a passage does not exist, return "".

Do NOT invent text.

Return ONLY valid JSON:

{
    "study_host": "...",

    "study_host_passage": "...",

    "host_virus_name": "...",

    "host_virus_passage": "...",

    "comparison_virus_name": "...",

    "comparison_host_or_sample": "...",

    "relationship_type": "...",

    "relationship_passage": "...",

    "reason": "..."
}
"""

    user_prompt = f"""
TARGET HOST:
{host}

VALIDATED HOST ALIASES:
{aliases_text}

TARGET VIRUS:
{virus}

PAPER:
{evidence_source}
"""

    messages = [
        {
            "role":
                "system",

            "content":
                system_prompt
        },
        {
            "role":
                "user",

            "content":
                user_prompt
        }
    ]

    print(
        "\nEvidence Agent analyzing:"
    )

    print(
        paper.get(
            "title",
            ""
        )
    )

    try:

        response = generate_text(
            messages,
            max_new_tokens=400,
            max_input_tokens=
                MAX_INPUT_TOKENS
        )

        # More robust than greedy regex.
        decoder = json.JSONDecoder()

        extraction = None

        for i, char in enumerate(
            response
        ):

            if char != "{":
                continue

            try:

                obj, _ = (
                    decoder.raw_decode(
                        response[i:]
                    )
                )

                if isinstance(
                    obj,
                    dict
                ):

                    extraction = obj
                    break

            except json.JSONDecodeError:
                continue

        if extraction is None:

            raise ValueError(
                "No valid JSON object returned."
            )

    except Exception as error:

        print(
            "WARNING: Evidence extraction "
            f"failed: {error}"
        )

        cleanup_gpu()

        extraction = {
            "study_host":
                "",

            "study_host_passage":
                "",

            "host_virus_name":
                "",

            "host_virus_passage":
                "",

            "comparison_virus_name":
                "",

            "comparison_host_or_sample":
                "",

            "relationship_type":
                "UNCLEAR",

            "relationship_passage":
                "",

            "reason":
                (
                    "Extraction failed: "
                    f"{error}"
                )
        }

    classification = (
        classify_relationship(
            extraction,
            host_aliases,
            virus,
            evidence_source
        )
    )

    extraction[
        "classification"
    ] = classification

    extraction[
        "host_aliases"
    ] = host_aliases

    extraction[
        "title"
    ] = paper.get(
        "title"
    )

    extraction[
        "pmid"
    ] = paper.get(
        "pmid"
    )

    extraction[
        "pmcid"
    ] = paper.get(
        "pmcid"
    )

    extraction[
        "doi"
    ] = paper.get(
        "doi"
    )

    extraction[
        "sources"
    ] = paper.get(
        "sources",
        []
    )

    extraction[
        "supporting_quote"
    ] = extraction.get(
        "relationship_passage",
        ""
    )

    cleanup_gpu()

    return extraction


# ============================================================
# RUN EVIDENCE AGENT
# ============================================================

def run_evidence_agent(
    host,
    virus
):

    # Search agent now resolves aliases too.
    search_result = (
        run_search_agent(
            host,
            virus
        )
    )

    papers = (
        search_result[
            "papers"
        ]
    )

    search_metadata = (
        search_result[
            "search_metadata"
        ]
    )

    host_aliases = (
        search_result.get(
            "host_aliases"
        )
        or get_host_aliases(
            host
        )
    )

    print(
        "\n"
        + "=" * 80
    )

    print(
        "STARTING HOST-VIRUS/"
        "COMPARISON-VIRUS ANALYSIS"
    )

    print(
        "=" * 80
    )

    print(
        "Validated host aliases:",
        host_aliases
    )

    results = []

    for paper in papers:

        result = analyze_paper(
            host,
            host_aliases,
            virus,
            paper
        )

        results.append(
            result
        )

    return {
        "evidence_results":
            results,

        "search_metadata":
            search_metadata,

        "host_aliases":
            host_aliases
    }


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":

    if len(sys.argv) != 3:

        print(
            'Usage: python evidence_agent.py '
            '"HOST" "VIRUS"'
        )

        sys.exit(1)

    host = sys.argv[1]
    virus = sys.argv[2]

    output = (
        run_evidence_agent(
            host,
            virus
        )
    )

    for result in output[
        "evidence_results"
    ]:

        print(
            "\n"
            + "=" * 80
        )

        print(
            "Paper:",
            result.get(
                "title"
            )
        )

        print(
            "Evidence class:",
            result.get(
                "classification"
            )
        )

        print(
            "\nStudy host:",
            result.get(
                "study_host"
            )
        )

        print(
            "Study host matches target:",
            result.get(
                "study_host_matches_target"
            )
        )

        print(
            "Verified target-host context:",
            result.get(
                "verified_target_host_context"
            )
        )

        print(
            "\nHost-associated virus:",
            result.get(
                "host_virus_name"
            )
        )

        print(
            "Host virus matches target:",
            result.get(
                "host_virus_matches_target"
            )
        )

        print(
            "\nComparison virus:",
            result.get(
                "comparison_virus_name"
            )
        )

        print(
            "Comparison virus matches target:",
            result.get(
                "comparison_virus_matches_target"
            )
        )

        print(
            "\nRelationship type:",
            result.get(
                "relationship_type"
            )
        )

        print(
            "\nStudy-host passage:"
        )

        print(
            result.get(
                "study_host_passage"
            )
        )

        print(
            "\nHost-virus passage:"
        )

        print(
            result.get(
                "host_virus_passage"
            )
        )

        print(
            "\nRelationship passage:"
        )

        print(
            result.get(
                "relationship_passage"
            )
        )

        print(
            "\nReason:"
        )

        print(
            result.get(
                "reason"
            )
        )