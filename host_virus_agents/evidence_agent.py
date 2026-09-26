import json
import re
import sys

from evidence_semantics import classify as classify_relationship, entity_equal, role_conflict


# Keep deterministic verification importable without loading a model or network
# packages. The live pipeline still uses the same implementations.
def run_search_agent(*args, **kwargs):
    from search_agent import run_search_agent as run
    return run(*args, **kwargs)


def generate_text(*args, **kwargs):
    from search_agent import generate_text as run
    return run(*args, **kwargs)


def cleanup_gpu():
    from search_agent import cleanup_gpu as run
    return run()


def get_host_aliases(*args, **kwargs):
    from taxonomy_aliases import get_host_aliases as run
    return run(*args, **kwargs)


def text_contains_alias(*args, **kwargs):
    from taxonomy_aliases import text_contains_alias as run
    return run(*args, **kwargs)


def name_matches_host(*args, **kwargs):
    from taxonomy_aliases import name_matches_host as run
    return run(*args, **kwargs)


def run_biological_context_agent(*args, **kwargs):
    from bioresearch_env.biological_context_agent import run_biological_context_agent as run
    return run(*args, **kwargs)


def context_for_prompt(*args, **kwargs):
    from bioresearch_env.biological_context_agent import context_for_prompt as run
    return run(*args, **kwargs)


def get_virus_taxonomy_context(*args, **kwargs):
    from bioresearch_env.virus_taxonomy_agent import get_virus_taxonomy_context as run
    return run(*args, **kwargs)


def get_pair_taxonomy(host, virus):
    """Resolve identity separately from the historical biological prior cache."""
    from bioresearch_env.biological_context_agent import get_taxonomy_context
    errors = []
    try:
        host_record = get_taxonomy_context(host)
    except Exception as error:
        host_record = {}
        errors.append("Host taxonomy: " + str(error))
    try:
        virus_record = get_virus_taxonomy_context(virus)
    except Exception as error:
        virus_record = {}
        errors.append("Virus taxonomy: " + str(error))
    return {
        "host": host_record, "virus": virus_record, "errors": errors,
        "resolved": bool(host_record.get("tax_id") and host_record.get("exact_name_match")
                         and host_record.get("rank") == "species"
                         and virus_record.get("resolved") and virus_record.get("rank") not in {"genus", "family", "order"}),
    }

from bioresearch_env.relationship_language import (
    has_direct_interaction_language,
    has_host_specific_study_context,
    has_human_clinical_context,
    relationship_type_supported_by_text,
    related_type_supported_by_text
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
    target,
    target_aliases=None,
):
    """
    Compare virus names using the original target name plus validated
    NCBI Taxonomy aliases/synonyms when available.

    This deliberately does not infer equivalence from family/genus alone.
    """

    return entity_equal(reported, target, target_aliases)


def text_contains_any_name(
    text,
    names,
):
    if not text:
        return False

    return any(
        text_contains_name(text, name)
        for name in (names or [])
        if name
    )


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
    virus_aliases=None,
    window=SNIPPET_WINDOW,
    max_snippets=MAX_VIRUS_SNIPPETS
):
    if not text:
        return ""

    text_lower = text.lower()

    search_names = [virus]
    search_names.extend(virus_aliases or [])

    # Deduplicate normalized aliases while keeping the requested virus first.
    unique_names = []
    seen = set()

    for name in search_names:
        if not name:
            continue

        key = normalize_entity(name)

        if not key or key in seen:
            continue

        seen.add(key)
        unique_names.append(name)

    candidates = []

    for search_name in unique_names:
        name_lower = search_name.lower()
        start = 0
        count = 0

        while count < MAX_OCCURRENCES_PER_TERM:
            position = text_lower.find(
                name_lower,
                start,
            )

            if position == -1:
                break

            left = max(0, position - window)
            right = min(
                len(text),
                position + len(search_name) + window,
            )

            snippet = text[left:right]

            # Exact requested target name is strongest; validated aliases are
            # still useful but receive a slightly smaller base score.
            score = (
                8
                if normalize_entity(search_name) == normalize_entity(virus)
                else 6
            )

            if text_contains_alias(snippet, host_aliases):
                score += 12

            lower = snippet.lower()

            for word in [
                "similar",
                "identity",
                "phylogen",
                "closest",
                "related",
                "detected",
                "isolated",
                "identified",
                "discovered",
                "sequenced",
                "infection",
                "infected",
            ]:
                if word in lower:
                    score += 1

            candidates.append((score, snippet))

            start = position + len(search_name)
            count += 1

    return rank_snippets(
        candidates,
        max_snippets,
    )


# ============================================================
# PYTHON EVIDENCE CLASSIFICATION
# ============================================================


HOST_ASSOCIATION_TYPES = {
    "DETECTED_IN",
    "IDENTIFIED_IN",
    "DISCOVERED_IN",
    "ISOLATED_FROM",
    "RECOVERED_FROM",
    "SEQUENCED_FROM",
    "FOUND_IN",
    "PRESENT_IN",
    "OBTAINED_FROM",
    "AMPLIFIED_FROM",
    "INFECTION_OF",
    "NATURAL_INFECTION",
    "EXPERIMENTAL_INFECTION",
    "REPLICATES_IN",
    "TRANSMITTED_TO",
    "ASSOCIATED_WITH",
    "STUDY_ASSOCIATION",
}

COMPARISON_RELATIONSHIP_TYPES = {
    "SEQUENCE_SIMILARITY",
    "PHYLOGENETIC_COMPARISON",
    "RELATED_VIRUS",
    "SAME_FAMILY",
    "SAME_GENUS",
}


def _normalize_host_relationship_type(value):
    value = str(value or "UNCLEAR").upper().strip()

    aliases = {
        "DISCOVERED": "DISCOVERED_IN",
        "IDENTIFIED": "IDENTIFIED_IN",
        "DETECTED": "DETECTED_IN",
        "ISOLATED": "ISOLATED_FROM",
        "RECOVERED": "RECOVERED_FROM",
        "SEQUENCED": "SEQUENCED_FROM",
    }

    return aliases.get(value, value)


def _normalize_comparison_relationship_type(value, passage=""):
    value = str(value or "UNCLEAR").upper().strip()

    aliases = {
        "PHYLOGENETIC_SIMILARITY": "PHYLOGENETIC_COMPARISON",
        "PHYLOGENETIC_RELATEDNESS": "PHYLOGENETIC_COMPARISON",
        "SEQUENCE_COMPARISON": "SEQUENCE_SIMILARITY",
        "SEQUENCE_IDENTITY": "SEQUENCE_SIMILARITY",
        "SIMILARITY": "SEQUENCE_SIMILARITY",
    }

    value = aliases.get(value, value)

    if value in COMPARISON_RELATIONSHIP_TYPES or value in {
        "REFERENCE_ONLY",
        "DIFFERENT_HOST",
        "NO_RELATION",
        "UNCLEAR",
        "MENTION_ONLY",
    }:
        return value

    lower = str(passage or "").lower()

    if "phylogen" in lower or "clustered with" in lower:
        return "PHYLOGENETIC_COMPARISON"

    if (
        "similar" in lower
        or "identity" in lower
        or "blast" in lower
    ):
        return "SEQUENCE_SIMILARITY"

    if "related" in lower or "closest" in lower:
        return "RELATED_VIRUS"

    return value

def _legacy_classify_relationship_v06(
    extraction,
    host_aliases,
    target_virus,
    evidence_source,
    biological_context=None,
    virus_aliases=None,
):
    """
    Deterministically classify one paper using two strictly separated
    evidence edges:

        host edge:
            target study host -> host-associated virus

        comparison edge:
            host-associated virus -> comparison/target virus

    v0.6 safety rule:
    virus-virus comparison language (sequence similarity, phylogeny,
    same family/genus) can NEVER establish a host edge.
    """

    biological_context = biological_context or {}

    # ========================================================
    # EXTRACTED FIELDS
    # ========================================================

    study_host = str(extraction.get("study_host", "") or "").strip()
    study_host_passage = str(
        extraction.get("study_host_passage", "") or ""
    ).strip()

    host_virus_name = str(
        extraction.get("host_virus_name", "") or ""
    ).strip()
    host_virus_passage = str(
        extraction.get("host_virus_passage", "") or ""
    ).strip()

    host_virus_relationship_type = _normalize_host_relationship_type(
        extraction.get("host_virus_relationship_type")
        or extraction.get("relationship_type", "UNCLEAR")
        or "UNCLEAR"
    )

    comparison_source_virus_name = str(
        extraction.get("comparison_source_virus_name", "") or ""
    ).strip()
    comparison_virus_name = str(
        extraction.get("comparison_virus_name", "") or ""
    ).strip()
    comparison_host_or_sample = str(
        extraction.get("comparison_host_or_sample", "") or ""
    ).strip()

    comparison_relationship_passage = str(
        extraction.get("comparison_relationship_passage")
        or extraction.get("relationship_passage", "")
        or ""
    ).strip()

    comparison_relationship_type = _normalize_comparison_relationship_type(
        extraction.get("comparison_relationship_type")
        or extraction.get("relationship_type", "UNCLEAR")
        or "UNCLEAR",
        comparison_relationship_passage,
    )

    # If the model omitted the comparison source, the host-associated
    # virus is the only safe source candidate. Record this inference.
    comparison_source_inferred = False
    if (
        not comparison_source_virus_name
        and host_virus_name
        and comparison_virus_name
    ):
        comparison_source_virus_name = host_virus_name
        comparison_source_inferred = True

    # ========================================================
    # VERIFY GROUNDED SOURCE TEXT
    # ========================================================

    study_host_passage_verified = verify_passage(
        study_host_passage,
        evidence_source,
    )
    host_virus_passage_verified = verify_passage(
        host_virus_passage,
        evidence_source,
    )
    comparison_relationship_passage_verified = verify_passage(
        comparison_relationship_passage,
        evidence_source,
    )

    # ========================================================
    # TARGET HOST RESOLUTION
    # ========================================================

    target_host_name = (
        biological_context
        .get("target_host", {})
        .get("scientific_name")
        or ""
    )

    human_target = (
        normalize_entity(target_host_name) == "homo sapiens"
        or any(
            normalize_entity(alias) == "homo sapiens"
            for alias in (host_aliases or [])
        )
    )

    model_host_name_matches_target = name_matches_host(
        study_host,
        host_aliases,
    )

    passage_contains_target_host = (
        study_host_passage_verified
        and text_contains_alias(
            study_host_passage,
            host_aliases,
        )
    )

    human_clinical_context = (
        human_target
        and study_host_passage_verified
        and (
            has_human_clinical_context(study_host_passage)
            or has_human_clinical_context(study_host)
        )
    )

    # CRITICAL v0.6 rule: an LLM-provided study_host string is not enough.
    # The grounded passage itself must establish the target host, except for
    # the explicitly handled Homo sapiens clinical normalization.
    study_host_matches_target = (
        passage_contains_target_host
        or human_clinical_context
    )

    verified_target_host_context = (
        study_host_passage_verified
        and study_host_matches_target
    )

    # ========================================================
    # VIRUS ENTITY RESOLUTION
    # ========================================================

    host_virus_matches_target = same_virus(
        host_virus_name,
        target_virus,
        target_aliases=virus_aliases,
    )

    comparison_virus_matches_target = same_virus(
        comparison_virus_name,
        target_virus,
        target_aliases=virus_aliases,
    )

    comparison_source_matches_host_virus = (
        bool(comparison_source_virus_name)
        and bool(host_virus_name)
        and same_virus(
            comparison_source_virus_name,
            host_virus_name,
        )
    )

    host_virus_passage_mentions_host_virus = (
        text_contains_name(host_virus_passage, host_virus_name)
        if host_virus_name
        else False
    )

    host_virus_passage_contains_target_host = (
        host_virus_passage_verified
        and text_contains_alias(
            host_virus_passage,
            host_aliases,
        )
    )

    target_virus_names = [target_virus]
    target_virus_names.extend(virus_aliases or [])

    comparison_relationship_mentions_target = text_contains_any_name(
        comparison_relationship_passage,
        target_virus_names,
    )

    comparison_relationship_mentions_source = (
        text_contains_name(
            comparison_relationship_passage,
            comparison_source_virus_name,
        )
        if comparison_source_virus_name
        else False
    )

    # ========================================================
    # STRICT EDGE TYPE SETS
    # ========================================================

    host_association_types = {
        "DETECTED_IN",
        "IDENTIFIED_IN",
        "DISCOVERED_IN",
        "ISOLATED_FROM",
        "RECOVERED_FROM",
        "SEQUENCED_FROM",
        "FOUND_IN",
        "PRESENT_IN",
        "OBTAINED_FROM",
        "AMPLIFIED_FROM",
        "INFECTION_OF",
        "NATURAL_INFECTION",
        "EXPERIMENTAL_INFECTION",
        "REPLICATES_IN",
        "TRANSMITTED_TO",
        "ASSOCIATED_WITH",
        "STUDY_ASSOCIATION",
    }

    # Only these host-edge types may use two passages from the same
    # host-specific discovery study.
    contextual_host_types = {
        "DETECTED_IN",
        "IDENTIFIED_IN",
        "DISCOVERED_IN",
        "SEQUENCED_FROM",
        "FOUND_IN",
        "PRESENT_IN",
        "STUDY_ASSOCIATION",
    }

    related_types = {
        "SEQUENCE_SIMILARITY",
        "PHYLOGENETIC_COMPARISON",
        "RELATED_VIRUS",
        "SAME_FAMILY",
        "SAME_GENUS",
    }

    # ========================================================
    # HOST EDGE SUPPORT
    # ========================================================

    combined_host_context = normalize_whitespace(
        study_host_passage + " " + host_virus_passage
    )

    host_specific_study_context = has_host_specific_study_context(
        combined_host_context
    )
    paper_has_host_specific_study_context = has_host_specific_study_context(
        evidence_source
    )

    # For ordinary direct evidence, the same grounded virus passage must
    # contain the target host (or clear human clinical context). This avoids
    # stitching together unrelated sections of a long paper.
    direct_passage_has_target_host = (
        host_virus_passage_contains_target_host
        or (
            human_target
            and has_human_clinical_context(host_virus_passage)
        )
    )

    direct_host_relationship_supported = (
        host_virus_relationship_type in host_association_types
        and host_virus_passage_verified
        and host_virus_passage_mentions_host_virus
        and direct_passage_has_target_host
        and relationship_type_supported_by_text(
            host_virus_relationship_type,
            host_virus_passage,
        )
    )

    # Cross-passage host evidence is allowed only for genuine host-specific
    # discovery designs and only for host-association relationship types.
    contextual_host_association_supported = (
        host_virus_relationship_type in contextual_host_types
        and verified_target_host_context
        and host_virus_passage_verified
        and host_virus_passage_mentions_host_virus
        and (
            host_specific_study_context
            or paper_has_host_specific_study_context
        )
    )

    human_clinical_association_supported = (
        human_target
        and host_virus_relationship_type in host_association_types
        and verified_target_host_context
        and host_virus_passage_verified
        and host_virus_passage_mentions_host_virus
        and (
            has_human_clinical_context(host_virus_passage)
            or has_human_clinical_context(combined_host_context)
            or has_human_clinical_context(evidence_source)
        )
        and (
            relationship_type_supported_by_text(
                host_virus_relationship_type,
                host_virus_passage,
            )
            or host_virus_relationship_type in contextual_host_types
        )
    )

    if direct_host_relationship_supported:
        host_edge_support_mode = "DIRECT"
    elif contextual_host_association_supported:
        host_edge_support_mode = "STUDY_CONTEXT"
    elif human_clinical_association_supported:
        host_edge_support_mode = "HUMAN_CLINICAL"
    else:
        host_edge_support_mode = "NONE"

    host_edge_verified = (
        host_virus_relationship_type in host_association_types
        and verified_target_host_context
        and host_virus_passage_verified
        and host_virus_passage_mentions_host_virus
        and host_edge_support_mode != "NONE"
    )

    # ========================================================
    # COMPARISON EDGE SUPPORT
    # ========================================================

    related_relationship_supported = (
        comparison_relationship_type in related_types
        and comparison_relationship_passage_verified
        and comparison_relationship_mentions_target
        and related_type_supported_by_text(
            comparison_relationship_type,
            comparison_relationship_passage,
        )
    )

    comparison_source_link_supported = (
        comparison_source_matches_host_virus
        and (
            comparison_relationship_mentions_source
            or not comparison_source_inferred
            or host_virus_passage_verified
        )
    )

    comparison_edge_verified = (
        bool(host_virus_name)
        and comparison_virus_matches_target
        and comparison_source_link_supported
        and related_relationship_supported
    )

    # ========================================================
    # EXPLICIT EDGE OBJECTS
    # ========================================================

    host_edge = {
        "host": study_host,
        "virus": host_virus_name,
        "relationship_type": host_virus_relationship_type,
        "study_host_passage": study_host_passage,
        "virus_passage": host_virus_passage,
        "study_host_passage_verified": study_host_passage_verified,
        "virus_passage_verified": host_virus_passage_verified,
        "model_host_name_matches_target": model_host_name_matches_target,
        "passage_contains_target_host": passage_contains_target_host,
        "target_host_match": study_host_matches_target,
        "target_virus_match": host_virus_matches_target,
        "support_mode": host_edge_support_mode,
        "verified": host_edge_verified,
    }

    comparison_edge = {
        "source_virus": comparison_source_virus_name,
        "comparison_virus": comparison_virus_name,
        "comparison_host_or_sample": comparison_host_or_sample,
        "relationship_type": comparison_relationship_type,
        "relationship_passage": comparison_relationship_passage,
        "passage_verified": comparison_relationship_passage_verified,
        "target_virus_match": comparison_virus_matches_target,
        "source_matches_host_virus": comparison_source_matches_host_virus,
        "source_inferred": comparison_source_inferred,
        "passage_mentions_source": comparison_relationship_mentions_source,
        "passage_mentions_target": comparison_relationship_mentions_target,
        "relationship_supported": related_relationship_supported,
        "verified": comparison_edge_verified,
    }

    extraction["host_edge"] = host_edge
    extraction["comparison_edge"] = comparison_edge

    # ========================================================
    # BACKWARDS-COMPATIBLE DIAGNOSTICS
    # ========================================================

    extraction["study_host_passage_verified"] = study_host_passage_verified
    extraction["study_host_matches_target"] = study_host_matches_target
    extraction["model_host_name_matches_target"] = model_host_name_matches_target
    extraction["passage_contains_target_host"] = passage_contains_target_host
    extraction["human_clinical_context"] = human_clinical_context
    extraction["verified_target_host_context"] = verified_target_host_context
    extraction["strong_study_host_context"] = verified_target_host_context
    extraction["host_virus_passage_verified"] = host_virus_passage_verified
    extraction["host_virus_matches_target"] = host_virus_matches_target
    extraction["comparison_virus_matches_target"] = comparison_virus_matches_target
    extraction["comparison_source_matches_host_virus"] = comparison_source_matches_host_virus
    extraction["comparison_relationship_passage_verified"] = comparison_relationship_passage_verified
    extraction["relationship_passage_verified"] = comparison_relationship_passage_verified
    extraction["comparison_relationship_mentions_target_virus"] = comparison_relationship_mentions_target
    extraction["relationship_mentions_target_virus"] = comparison_relationship_mentions_target
    extraction["comparison_relationship_mentions_source_virus"] = comparison_relationship_mentions_source
    extraction["host_specific_study_context"] = host_specific_study_context
    extraction["paper_has_host_specific_study_context"] = paper_has_host_specific_study_context
    extraction["direct_host_relationship_supported"] = direct_host_relationship_supported
    extraction["contextual_host_association_supported"] = contextual_host_association_supported
    extraction["human_clinical_association_supported"] = human_clinical_association_supported
    extraction["host_association_supported"] = host_edge_verified
    extraction["related_relationship_supported"] = related_relationship_supported
    extraction["comparison_edge_verified"] = comparison_edge_verified
    extraction["reported_virus_name"] = host_virus_name
    extraction["reported_host_or_sample"] = study_host
    extraction["reported_virus_matches_target"] = host_virus_matches_target
    extraction["reported_host_matches_target"] = study_host_matches_target
    extraction["reported_virus_passage_verified"] = host_virus_passage_verified

    # ========================================================
    # EXACT SUPPORT
    # ========================================================

    if host_edge_verified and host_virus_matches_target:
        extraction["relationship_type"] = host_virus_relationship_type
        extraction["relationship_passage"] = host_virus_passage
        extraction["classification_basis"] = (
            "VERIFIED_HOST_EDGE_TO_TARGET_VIRUS"
        )
        return "EXACT_SUPPORT"

    # ========================================================
    # TARGET HOST RELATED
    # ========================================================

    if (
        host_edge_verified
        and not host_virus_matches_target
        and comparison_edge_verified
    ):
        extraction["relationship_type"] = comparison_relationship_type
        extraction["relationship_passage"] = comparison_relationship_passage
        extraction["classification_basis"] = (
            "VERIFIED_HOST_EDGE_PLUS_COMPARISON_EDGE"
        )
        return "TARGET_HOST_RELATED"

    # ========================================================
    # TARGET VIRUS EXPLICITLY ASSOCIATED WITH ANOTHER HOST
    # ========================================================

    non_target_host_direct_edge = (
        not verified_target_host_context
        and host_virus_matches_target
        and host_virus_relationship_type in host_association_types
        and host_virus_passage_verified
        and host_virus_passage_mentions_host_virus
        and relationship_type_supported_by_text(
            host_virus_relationship_type,
            host_virus_passage,
        )
    )

    explicit_comparison_other_host = (
        comparison_virus_matches_target
        and comparison_relationship_type == "DIFFERENT_HOST"
        and comparison_relationship_passage_verified
        and comparison_relationship_mentions_target
        and bool(comparison_host_or_sample)
    )

    if non_target_host_direct_edge or explicit_comparison_other_host:
        extraction["classification_basis"] = (
            "TARGET_VIRUS_VERIFIED_IN_NON_TARGET_HOST"
        )
        return "VIRUS_OTHER_HOST"

    # ========================================================
    # INCOMPLETE RELATED CLAIM
    # ========================================================

    if comparison_relationship_type in related_types:
        extraction["relationship_type"] = comparison_relationship_type
        extraction["relationship_passage"] = comparison_relationship_passage
        extraction["classification_basis"] = (
            "RELATED_LANGUAGE_BUT_INCOMPLETE_EDGE_CHAIN"
        )
        return "NO_SUPPORT"

    # A comparison-type label is never valid as a host edge. If the model
    # produced one, fail closed rather than turning it into exact support.
    if host_virus_relationship_type in related_types:
        extraction["classification_basis"] = (
            "INVALID_COMPARISON_TYPE_ON_HOST_EDGE"
        )
        return "NO_SUPPORT"

    # ========================================================
    # EXPLICIT NON-SUPPORT
    # ========================================================

    if host_virus_relationship_type == "MENTION_ONLY":
        extraction["relationship_type"] = "MENTION_ONLY"
        extraction["relationship_passage"] = host_virus_passage
        extraction["classification_basis"] = "MENTION_ONLY"
        return "MENTION_ONLY"

    if comparison_relationship_type == "REFERENCE_ONLY":
        extraction["relationship_type"] = "MENTION_ONLY"
        extraction["relationship_passage"] = comparison_relationship_passage
        extraction["classification_basis"] = "REFERENCE_ONLY"
        return "MENTION_ONLY"

    if (
        host_virus_relationship_type in {"DIFFERENT_HOST", "NO_RELATION"}
        or comparison_relationship_type in {"DIFFERENT_HOST", "NO_RELATION"}
    ):
        if non_target_host_direct_edge or explicit_comparison_other_host:
            extraction["classification_basis"] = "EXPLICIT_NON_TARGET_HOST"
            return "VIRUS_OTHER_HOST"

        extraction["classification_basis"] = "EXPLICIT_NO_RELATION"
        return "NO_SUPPORT"

    # ========================================================
    # GENUINE UNCERTAINTY
    # ========================================================

    if (
        host_virus_relationship_type == "UNCLEAR"
        or comparison_relationship_type == "UNCLEAR"
    ):
        if (
            verified_target_host_context
            and (
                host_virus_matches_target
                or comparison_virus_matches_target
            )
        ):
            extraction["classification_basis"] = (
                "TARGET_ENTITY_PRESENT_BUT_EDGE_UNRESOLVED"
            )
            return "UNCLEAR"

    extraction["classification_basis"] = "NO_VERIFIED_SUPPORT"
    return "NO_SUPPORT"


# ============================================================
# ANALYZE ONE PAPER
# ============================================================

def _extract_json_object(response):
    decoder = json.JSONDecoder()

    for i, char in enumerate(response or ""):
        if char != "{":
            continue

        try:
            obj, _ = decoder.raw_decode(response[i:])

            if isinstance(obj, dict):
                return obj
        except json.JSONDecodeError:
            continue

    raise ValueError("No valid JSON object returned.")


def _empty_extraction(reason=""):
    return {
        "study_host": "",
        "study_host_passage": "",
        "host_virus_name": "",
        "host_virus_passage": "",
        "host_virus_relationship_type": "UNCLEAR",
        "comparison_source_virus_name": "",
        "comparison_virus_name": "",
        "comparison_host_or_sample": "",
        "comparison_relationship_type": "UNCLEAR",
        "comparison_relationship_passage": "",
        "reason": reason,
    }


def _extraction_rescue_reason(
    extraction,
    evidence_source,
    target_virus,
    target_virus_names,
):
    """Return a short rescue reason, or an empty string if no rescue is needed."""

    extraction = extraction or {}
    conflict = role_conflict(extraction)
    if conflict:
        return conflict

    core_fields = [
        extraction.get("study_host"),
        extraction.get("host_virus_name"),
        extraction.get("comparison_virus_name"),
    ]

    if not any(str(value or "").strip() for value in core_fields):
        return "EMPTY_EXTRACTION"

    host_type = _normalize_host_relationship_type(
        extraction.get("host_virus_relationship_type")
        or extraction.get("relationship_type", "UNCLEAR")
    )

    # A virus-virus comparison label on the host edge is the exact failure
    # mode that caused v0.5 to promote comparison viruses to KNOWN.
    if host_type in COMPARISON_RELATIONSHIP_TYPES:
        return "COMPARISON_TYPE_ON_HOST_EDGE"

    target_visible = text_contains_any_name(
        evidence_source,
        target_virus_names,
    )

    extracted_virus = (
        str(extraction.get("host_virus_name") or "").strip()
        or str(extraction.get("comparison_virus_name") or "").strip()
    )

    if target_visible and not extracted_virus:
        return "TARGET_VISIBLE_BUT_NOT_EXTRACTED"

    # If the requested target has been put on the host edge with a comparison
    # relationship, force a role-correction pass.
    host_virus_name = str(extraction.get("host_virus_name") or "").strip()
    if (
        host_virus_name
        and same_virus(
            host_virus_name,
            target_virus,
            target_aliases=target_virus_names,
        )
        and host_type in COMPARISON_RELATIONSHIP_TYPES
    ):
        return "TARGET_PROMOTED_FROM_COMPARISON"

    return ""


def _run_extraction(messages, max_new_tokens=500):
    response = generate_text(
        messages,
        max_new_tokens=max_new_tokens,
        max_input_tokens=MAX_INPUT_TOKENS,
    )
    try:
        result = _extract_json_object(response)
    except Exception as error:
        error.raw_response = response
        raise
    result["raw_model_response"] = response
    return result


def _flexible_name_matches(text, name):
    """Yield raw-text regex matches tolerant of punctuation/hyphen changes."""
    if not text or not name:
        return []

    tokens = normalize_entity(name).split()
    if not tokens:
        return []

    pattern = r"(?<!\w)" + r"[\W_]+".join(
        re.escape(token) for token in tokens
    ) + r"(?!\w)"

    return list(
        re.finditer(
            pattern,
            text,
            flags=re.IGNORECASE,
        )
    )


def _source_window(source_text, start, end, window=650):
    """Return an exact substring from source_text around a match."""
    left = max(0, start - window)
    right = min(len(source_text), end + window)
    return source_text[left:right].strip()


def _best_grounded_entity_passage(
    source_text,
    names,
    host_aliases=None,
    prefer_direct=False,
    prefer_study_context=False,
    prefer_related=False,
):
    """Find the highest-scoring exact source window containing an entity."""

    candidates = []
    seen_names = set()

    for index, name in enumerate(names or []):
        key = normalize_entity(name)
        if not key or key in seen_names:
            continue
        seen_names.add(key)

        for match in _flexible_name_matches(source_text, name):
            snippet = _source_window(
                source_text,
                match.start(),
                match.end(),
            )

            score = 20 if index == 0 else 15

            if host_aliases and text_contains_alias(snippet, host_aliases):
                score += 15

            if prefer_direct and has_direct_interaction_language(snippet):
                score += 10

            if prefer_study_context and has_host_specific_study_context(snippet):
                score += 10

            if prefer_related and related_type_supported_by_text(
                "RELATED_VIRUS",
                snippet,
            ):
                score += 8

            # Shorter windows are easier to verify and less likely to combine
            # unrelated parts of the paper when scores tie.
            candidates.append((score, -len(snippet), snippet))

    if not candidates:
        return ""

    candidates.sort(reverse=True)
    return candidates[0][2]


def _best_grounded_comparison_passage(
    source_text,
    source_virus,
    comparison_virus,
    comparison_aliases=None,
):
    target_names = [comparison_virus]
    target_names.extend(comparison_aliases or [])

    candidates = []

    for name in target_names:
        for match in _flexible_name_matches(source_text, name):
            snippet = _source_window(
                source_text,
                match.start(),
                match.end(),
                window=800,
            )

            score = 10

            if source_virus and text_contains_name(snippet, source_virus):
                score += 20

            if (
                "phylogen" in snippet.lower()
                or "similar" in snippet.lower()
                or "identity" in snippet.lower()
                or "closest" in snippet.lower()
                or "related" in snippet.lower()
                or "blast" in snippet.lower()
            ):
                score += 20

            candidates.append((score, -len(snippet), snippet))

    if not candidates:
        return ""

    candidates.sort(reverse=True)
    return candidates[0][2]


def _ground_extraction_passages(
    extraction,
    evidence_source,
    host_aliases,
    target_virus,
    virus_aliases,
):
    """
    Replace model-generated/reconstructed quotations with exact spans copied
    from evidence_source. Entity roles and relationship labels are NOT changed.
    """

    extraction = dict(extraction or {})
    replacements = []

    # Study-host passage: prefer the scientific name first, then validated
    # aliases. This independently verifies that the target host is actually
    # present in the supplied paper evidence.
    host_names = list(host_aliases or [])

    existing = str(extraction.get("study_host_passage") or "")
    existing_ok = (
        verify_passage(existing, evidence_source)
        and text_contains_alias(existing, host_aliases)
    )

    if not existing_ok:
        grounded = _best_grounded_entity_passage(
            evidence_source,
            host_names,
            host_aliases=host_aliases,
            prefer_direct=True,
            prefer_study_context=True,
        )
        if grounded:
            extraction["study_host_passage"] = grounded
            replacements.append("study_host_passage")

    # Host-associated-virus passage.
    host_virus_name = str(extraction.get("host_virus_name") or "").strip()
    if host_virus_name:
        existing = str(extraction.get("host_virus_passage") or "")
        existing_ok = (
            verify_passage(existing, evidence_source)
            and text_contains_name(existing, host_virus_name)
        )

        if not existing_ok:
            grounded = _best_grounded_entity_passage(
                evidence_source,
                [host_virus_name],
                host_aliases=host_aliases,
                prefer_direct=True,
                prefer_study_context=True,
            )
            if grounded:
                extraction["host_virus_passage"] = grounded
                replacements.append("host_virus_passage")

    # Virus-virus comparison passage.
    comparison_virus = str(
        extraction.get("comparison_virus_name") or ""
    ).strip()
    comparison_source = str(
        extraction.get("comparison_source_virus_name")
        or extraction.get("host_virus_name")
        or ""
    ).strip()

    if comparison_virus:
        existing = str(
            extraction.get("comparison_relationship_passage")
            or extraction.get("relationship_passage")
            or ""
        )
        existing_ok = (
            verify_passage(existing, evidence_source)
            and text_contains_name(existing, comparison_virus)
        )

        comparison_aliases = (
            virus_aliases
            if same_virus(
                comparison_virus,
                target_virus,
                target_aliases=virus_aliases,
            )
            else []
        )

        if not existing_ok:
            grounded = _best_grounded_comparison_passage(
                evidence_source,
                comparison_source,
                comparison_virus,
                comparison_aliases=comparison_aliases,
            )
            if grounded:
                extraction["comparison_relationship_passage"] = grounded
                replacements.append("comparison_relationship_passage")

    extraction["grounding_applied"] = bool(replacements)
    extraction["grounding_replacements"] = replacements
    return extraction


def _run_rescue_extraction(
    host,
    host_aliases,
    virus,
    virus_aliases,
    paper,
    abstract,
    host_context,
    virus_context,
    rescue_reason,
):
    """Focused second pass used only for empty or role-conflicted extraction."""

    alias_text = "\n".join(
        f"- {alias}" for alias in (host_aliases or [])
    )
    virus_alias_text = "\n".join(
        f"- {alias}" for alias in (virus_aliases or [])[:20]
    )

    rescue_system_prompt = """
You are correcting a failed scientific evidence extraction.

The paper can contain TWO DIFFERENT relations:

1. HOST EDGE: study host -> host-associated virus
2. COMPARISON EDGE: host-associated virus -> comparison virus

CRITICAL ROLE CONSTRAINTS:
- host_virus_relationship_type MUST describe a HOST-to-VIRUS association.
- NEVER use SEQUENCE_SIMILARITY, PHYLOGENETIC_COMPARISON, RELATED_VIRUS,
  SAME_FAMILY, or SAME_GENUS as host_virus_relationship_type.
- If the TARGET VIRUS appears only in text such as "similar to", "closest
  to", "phylogenetically related to", or a taxonomy comparison, it belongs
  in comparison_virus_name, NOT host_virus_name.
- comparison_source_virus_name is the virus being compared TO the target.
- Do not promote the requested target virus to host_virus_name simply
  because it is the requested target.

Valid HOST relationship types:
DETECTED_IN
IDENTIFIED_IN
DISCOVERED_IN
ISOLATED_FROM
RECOVERED_FROM
SEQUENCED_FROM
FOUND_IN
PRESENT_IN
OBTAINED_FROM
AMPLIFIED_FROM
INFECTION_OF
NATURAL_INFECTION
EXPERIMENTAL_INFECTION
REPLICATES_IN
TRANSMITTED_TO
ASSOCIATED_WITH
STUDY_ASSOCIATION
MENTION_ONLY
DIFFERENT_HOST
NO_RELATION
UNCLEAR

Valid COMPARISON relationship types:
SEQUENCE_SIMILARITY
PHYLOGENETIC_COMPARISON
RELATED_VIRUS
SAME_FAMILY
SAME_GENUS
REFERENCE_ONLY
DIFFERENT_HOST
NO_RELATION
UNCLEAR

Do not invent passages. The Python verifier will ground exact source spans
separately, so prioritize correct ENTITY ROLES and RELATIONSHIP TYPES.

Return ONLY JSON:
{
  "study_host": "...",
  "study_host_passage": "...",
  "host_virus_name": "...",
  "host_virus_passage": "...",
  "host_virus_relationship_type": "...",
  "comparison_source_virus_name": "...",
  "comparison_virus_name": "...",
  "comparison_host_or_sample": "...",
  "comparison_relationship_type": "...",
  "comparison_relationship_passage": "...",
  "reason": "..."
}
"""

    rescue_user_prompt = f"""
RESCUE REASON:
{rescue_reason}

TARGET HOST:
{host}

VALIDATED HOST ALIASES:
{alias_text}

TARGET VIRUS:
{virus}

VALIDATED TARGET-VIRUS ALIASES:
{virus_alias_text}

TITLE:
{paper.get('title', '')}

ABSTRACT:
{abstract}

HOST CONTEXT:
{host_context}

TARGET VIRUS CONTEXT:
{virus_context}
"""

    return _run_extraction(
        [
            {"role": "system", "content": rescue_system_prompt},
            {"role": "user", "content": rescue_user_prompt},
        ],
        max_new_tokens=500,
    )


def analyze_paper(
    host,
    host_aliases,
    virus,
    paper,
    biological_context=None,
):
    full_text = paper.get("full_text", "") or ""
    abstract = (paper.get("abstract", "") or "")[:MAX_ABSTRACT_CHARS]

    context_tax_id = (
        (biological_context or {})
        .get("target_virus", {})
        .get("tax_id")
    )

    virus_taxonomy = get_virus_taxonomy_context(
        virus,
        tax_id=context_tax_id,
    )

    virus_aliases = list(virus_taxonomy.get("aliases", []) or [])
    if not any(
        normalize_entity(alias) == normalize_entity(virus)
        for alias in virus_aliases
    ):
        virus_aliases.insert(0, virus)

    host_context = get_host_context_snippets(
        full_text,
        host_aliases,
    )
    virus_context = get_virus_context_snippets(
        full_text,
        host_aliases,
        virus,
        virus_aliases=virus_aliases,
    )

    evidence_source = f"""
TITLE:
{paper.get('title', '')}

ABSTRACT:
{abstract}

HOST/STUDY CONTEXT:
{host_context}

TARGET-VIRUS CONTEXT:
{virus_context}
"""

    aliases_text = "\n".join(f"- {alias}" for alias in host_aliases)
    virus_aliases_text = "\n".join(
        f"- {alias}" for alias in virus_aliases[:20]
    )
    biological_context_text = (
        context_for_prompt(biological_context)
        if biological_context
        else "{}"
    )

    system_prompt = """
You are a scientific evidence extraction agent.
Extract a complete literal assertion binding subject host, predicate and object
virus. A name fragment is not an assertion. Preserve other hosts and other
viruses; never swap comparison endpoints to make them fit the target.
Report host_scope and virus_scope as SPECIES, GENUS, FAMILY or UNRESOLVED.
For cross-passage discovery report shared_specimen only if BOTH quotations
explicitly identify the same named specimen collected from the stated host.
General virome study context, environmental source, exposure to viral proteins,
sequence similarity, and references to earlier work are not exact support.
Do not guess an omitted full scientific name from an ambiguous abbreviation.

You extract evidence; you do NOT make the final classification.

Separate TWO relations:
A. study host -> host-associated virus
B. host-associated virus -> comparison virus

CRITICAL ROLE CONSTRAINTS:
- host_virus_relationship_type MUST be a HOST-to-VIRUS relationship.
- NEVER assign SEQUENCE_SIMILARITY, PHYLOGENETIC_COMPARISON,
  RELATED_VIRUS, SAME_FAMILY, or SAME_GENUS to host_virus_relationship_type.
- If the requested TARGET VIRUS appears only as a sequence/phylogenetic/
  taxonomic comparison, it MUST be comparison_virus_name, not
  host_virus_name.
- comparison_source_virus_name is the host-associated virus on the left
  side of that virus-virus comparison.
- A target name or validated synonym does not by itself prove a host link.

Valid HOST relationship types:
DETECTED_IN
IDENTIFIED_IN
DISCOVERED_IN
ISOLATED_FROM
RECOVERED_FROM
SEQUENCED_FROM
FOUND_IN
PRESENT_IN
OBTAINED_FROM
AMPLIFIED_FROM
INFECTION_OF
NATURAL_INFECTION
EXPERIMENTAL_INFECTION
REPLICATES_IN
TRANSMITTED_TO
ASSOCIATED_WITH
STUDY_ASSOCIATION
MENTION_ONLY
DIFFERENT_HOST
NO_RELATION
UNCLEAR

Valid COMPARISON relationship types:
SEQUENCE_SIMILARITY
PHYLOGENETIC_COMPARISON
RELATED_VIRUS
SAME_FAMILY
SAME_GENUS
REFERENCE_ONLY
DIFFERENT_HOST
NO_RELATION
UNCLEAR

A host-specific virome/transcriptome/HTS discovery paper can establish a
host-associated virus across separate passages, but a general paper that
mentions the host and virus in unrelated sections cannot.

BIOLOGICAL CONTEXT is a PRIOR only, never evidence.
VALIDATED TARGET-VIRUS ALIASES are entity-resolution hints only.

For Homo sapiens, clearly human patients/clinical specimens may normalize
to Homo sapiens. A human cell line by itself is not natural-host evidence.

Return ONLY JSON:
{
    "study_host": "...",
    "host_scope": "SPECIES|GENUS|FAMILY|UNRESOLVED",
    "virus_scope": "SPECIES|GENUS|FAMILY|UNRESOLVED",
    "shared_specimen": "literal specimen identifier or empty string",
    "study_host_passage": "...",
    "host_virus_name": "...",
    "host_virus_passage": "...",
    "host_virus_relationship_type": "...",
    "comparison_source_virus_name": "...",
    "comparison_virus_name": "...",
    "comparison_host_or_sample": "...",
    "comparison_relationship_type": "...",
    "comparison_relationship_passage": "...",
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

VALIDATED TARGET-VIRUS ALIASES:
{virus_aliases_text}

BIOLOGICAL CONTEXT PRIOR:
{biological_context_text}

PAPER:
{evidence_source}
"""

    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": user_prompt},
    ]

    print("\nEvidence Agent analyzing:")
    print(paper.get("title", ""))

    extraction_status = "SUCCESS"
    extraction_attempts = 1
    extraction_error = ""
    rescue_reason = ""
    extraction_trace = []

    try:
        extraction = _run_extraction(messages)
        extraction_trace.append({"attempt": "primary", "output": dict(extraction)})
    except Exception as error:
        extraction = _empty_extraction(
            f"Primary extraction failed: {error}"
        )
        extraction_error = str(error)
        extraction_trace.append({"attempt": "primary", "error": str(error), "raw_response": getattr(error, "raw_response", None)})

    target_virus_names = [virus]
    target_virus_names.extend(virus_aliases)

    rescue_reason = _extraction_rescue_reason(
        extraction,
        evidence_source,
        virus,
        target_virus_names,
    )
    if not rescue_reason:
        preview = dict(extraction)
        if classify_relationship(preview, host_aliases, virus, evidence_source,
                                 biological_context=biological_context,
                                 virus_aliases=virus_aliases) == "UNCLEAR":
            rescue_reason = "UNBOUND_DIRECTED_ASSERTION: " + preview["classification_basis"]

    if rescue_reason:
        extraction_attempts = 2
        print(
            "Evidence extraction needs rescue: "
            f"{rescue_reason}"
        )

        try:
            rescue = _run_rescue_extraction(
                host,
                host_aliases,
                virus,
                virus_aliases,
                paper,
                abstract,
                host_context,
                virus_context,
                rescue_reason,
            )

            second_reason = _extraction_rescue_reason(
                rescue,
                evidence_source,
                virus,
                target_virus_names,
            )
            extraction_trace.append({"attempt": "rescue", "output": dict(rescue), "residual_warning": second_reason})

            # A non-empty rescue is still useful even if it has a residual
            # role warning; the deterministic edge verifier will fail closed.
            if not any(
                str(rescue.get(field) or "").strip()
                for field in (
                    "study_host",
                    "host_virus_name",
                    "comparison_virus_name",
                )
            ):
                extraction_status = "FAILED"
                extraction_error = (
                    extraction_error
                    or "Rescue extraction returned no usable entity fields."
                )
                extraction = rescue
            else:
                extraction_status = "RESCUED"
                extraction_error = ""
                extraction = rescue
                if second_reason:
                    extraction["residual_rescue_warning"] = second_reason

        except Exception as error:
            extraction_status = "FAILED"
            extraction_error = str(error)
            extraction_trace.append({"attempt": "rescue", "error": str(error), "raw_response": getattr(error, "raw_response", None)})
            cleanup_gpu()

    # Preserve proposed quotations and roles. Invalid quotations are uncertainty;
    # the verifier expands only literal source spans within assertion boundaries.
    extraction["extraction_status"] = extraction_status
    extraction["extraction_error"] = extraction_error
    extraction["virus_taxonomy_resolved"] = virus_taxonomy.get("resolved", False)
    extraction["evidence_source"] = evidence_source
    extraction["source_paper"] = dict(paper)
    extraction["virus_taxonomy"] = virus_taxonomy
    extraction["extraction_trace"] = extraction_trace

    classification = classify_relationship(
        extraction,
        host_aliases,
        virus,
        evidence_source,
        biological_context=biological_context,
        virus_aliases=virus_aliases,
    )

    extraction["classification"] = classification
    extraction["host_aliases"] = host_aliases
    extraction["target_virus_aliases"] = virus_aliases
    extraction["target_virus_tax_id"] = virus_taxonomy.get("tax_id")
    extraction["target_virus_scientific_name"] = virus_taxonomy.get(
        "scientific_name"
    )
    extraction["virus_taxonomy_resolved"] = virus_taxonomy.get(
        "resolved",
        False,
    )
    extraction["extraction_status"] = extraction_status
    extraction["extraction_attempts"] = extraction_attempts
    extraction["extraction_error"] = extraction_error
    extraction["rescue_reason"] = rescue_reason

    extraction["title"] = paper.get("title")
    extraction["pmid"] = paper.get("pmid")
    extraction["pmcid"] = paper.get("pmcid")
    extraction["doi"] = paper.get("doi")
    extraction["sources"] = paper.get("sources", [])

    extraction["biological_prior"] = (
        (biological_context or {})
        .get("biological_prior", {})
        .get("status")
    )

    if classification == "EXACT_SUPPORT":
        supporting_quote = extraction.get("host_virus_passage", "")
    else:
        supporting_quote = extraction.get(
            "comparison_relationship_passage",
            "",
        )

    extraction["supporting_quote"] = supporting_quote

    cleanup_gpu()
    return extraction


# ============================================================
# RUN EVIDENCE AGENT
# ============================================================

def run_evidence_agent(
    host,
    virus
):

    biological_context = (
        run_biological_context_agent(
            host,
            virus
        )
    )

    # Search agent resolves aliases and retrieves/ranks papers.
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
    search_metadata["taxonomy_resolution"] = get_pair_taxonomy(host, virus)
    search_metadata["taxonomy_resolved"] = search_metadata["taxonomy_resolution"]["resolved"]

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

    print(
        "Biological prior:",
        biological_context
        .get(
            "biological_prior",
            {}
        )
        .get(
            "status",
            "UNKNOWN"
        )
    )

    results = []

    for paper in papers:

        result = analyze_paper(
            host,
            host_aliases,
            virus,
            paper,
            biological_context=
                biological_context
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
            host_aliases,

        "biological_context":
            biological_context
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
            "\nHost-virus relationship type:",
            result.get(
                "host_virus_relationship_type"
            )
        )

        print(
            "Comparison relationship type:",
            result.get(
                "comparison_relationship_type"
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
            "\nComparison relationship passage:"
        )

        print(
            result.get(
                "comparison_relationship_passage"
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
