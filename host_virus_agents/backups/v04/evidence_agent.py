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

from bioresearch_env.biological_context_agent import (
    run_biological_context_agent,
    context_for_prompt
)

from bioresearch_env.virus_taxonomy_agent import (
    get_virus_taxonomy_context,
)

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

    reported_norm = normalize_entity(reported)

    if not reported_norm:
        return False

    candidate_names = [target]
    candidate_names.extend(target_aliases or [])

    seen = set()

    for candidate in candidate_names:
        candidate_norm = normalize_entity(candidate)

        if (
            not candidate_norm
            or candidate_norm in seen
        ):
            continue

        seen.add(candidate_norm)

        if reported_norm == candidate_norm:
            return True

        # Allow a canonical name followed by descriptive text. Keep the
        # threshold reasonably long so generic names do not over-match.
        if (
            len(candidate_norm) >= 8
            and candidate_norm in reported_norm
        ):
            return True

        if (
            len(reported_norm) >= 8
            and reported_norm in candidate_norm
        ):
            return True

    return False


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

def classify_relationship(
    extraction,
    host_aliases,
    target_virus,
    evidence_source,
    biological_context=None,
    virus_aliases=None,
):
    """
    Deterministically classify one paper using two explicit evidence
    edges:

        host edge:
            study host -> host-associated virus

        comparison edge:
            host-associated virus -> comparison/target virus

    This separation prevents a comparison virus from being treated as a
    virus of the target host merely because it is mentioned in the same
    paper, while still allowing valid related-virus evidence chains.
    """

    biological_context = biological_context or {}

    # ========================================================
    # EXTRACTED FIELDS
    # ========================================================

    study_host = (
        extraction.get("study_host", "")
        or ""
    )

    study_host_passage = (
        extraction.get("study_host_passage", "")
        or ""
    )

    host_virus_name = (
        extraction.get("host_virus_name", "")
        or ""
    )

    host_virus_passage = (
        extraction.get("host_virus_passage", "")
        or ""
    )

    host_virus_relationship_type = (
        extraction.get("host_virus_relationship_type")
        or extraction.get("relationship_type", "UNCLEAR")
        or "UNCLEAR"
    ).upper().strip()

    comparison_source_virus_name = (
        extraction.get("comparison_source_virus_name", "")
        or ""
    )

    comparison_virus_name = (
        extraction.get("comparison_virus_name", "")
        or ""
    )

    comparison_host_or_sample = (
        extraction.get("comparison_host_or_sample", "")
        or ""
    )

    comparison_relationship_type = (
        extraction.get("comparison_relationship_type")
        or extraction.get("relationship_type", "UNCLEAR")
        or "UNCLEAR"
    ).upper().strip()

    comparison_relationship_passage = (
        extraction.get("comparison_relationship_passage")
        or extraction.get("relationship_passage", "")
        or ""
    )

    # If the model omitted the comparison source, the only legitimate
    # source candidate in this extraction is the host-associated virus.
    # Record the inference so diagnostics distinguish it from an explicit
    # model extraction.
    comparison_source_inferred = False

    if (
        not comparison_source_virus_name
        and host_virus_name
        and comparison_virus_name
    ):
        comparison_source_virus_name = host_virus_name
        comparison_source_inferred = True

    # ========================================================
    # VERIFY QUOTED TEXT
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
    # RESOLVE TARGET HOST
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
            for alias in host_aliases
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

    study_host_matches_target = (
        name_matches_host(
            study_host,
            host_aliases,
        )
        or (
            study_host_passage_verified
            and text_contains_alias(
                study_host_passage,
                host_aliases,
            )
        )
        or human_clinical_context
    )

    verified_target_host_context = (
        study_host_matches_target
        and study_host_passage_verified
    )

    # ========================================================
    # RESOLVE VIRUSES
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
        text_contains_name(
            host_virus_passage,
            host_virus_name,
        )
        if host_virus_name
        else False
    )

    target_virus_names = [target_virus]
    target_virus_names.extend(virus_aliases or [])

    comparison_relationship_mentions_target = (
        text_contains_any_name(
            comparison_relationship_passage,
            target_virus_names,
        )
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
    # RELATIONSHIP TYPES
    # ========================================================

    direct_types = {
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

    related_types = {
        "SEQUENCE_SIMILARITY",
        "PHYLOGENETIC_COMPARISON",
        "RELATED_VIRUS",
        "SAME_FAMILY",
        "SAME_GENUS",
    }

    negative_host_types = {
        "MENTION_ONLY",
        "DIFFERENT_HOST",
        "NO_RELATION",
    }

    # ========================================================
    # HOST EDGE SUPPORT
    # ========================================================

    combined_host_context = normalize_whitespace(
        study_host_passage
        + " "
        + host_virus_passage
    )

    host_virus_passage_has_direct_language = (
        has_direct_interaction_language(
            host_virus_passage
        )
    )

    combined_context_has_direct_language = (
        has_direct_interaction_language(
            combined_host_context
        )
    )

    host_specific_study_context = (
        has_host_specific_study_context(
            combined_host_context
        )
    )

    # v0.4: also inspect the supplied paper evidence as a whole. In
    # discovery papers, the virome/transcriptome framing may occur in the
    # title or abstract while the verified host and virus passages occur
    # elsewhere.
    paper_has_host_specific_study_context = (
        has_host_specific_study_context(
            evidence_source
        )
    )

    direct_host_relationship_supported = (
        host_virus_relationship_type in direct_types
        and (
            relationship_type_supported_by_text(
                host_virus_relationship_type,
                host_virus_passage,
            )
            or relationship_type_supported_by_text(
                host_virus_relationship_type,
                combined_host_context,
            )
        )
    )

    contextual_host_association_supported = (
        verified_target_host_context
        and host_virus_passage_verified
        and host_virus_passage_mentions_host_virus
        and host_virus_relationship_type
            not in negative_host_types
        and (
            host_specific_study_context
            or paper_has_host_specific_study_context
        )
    )

    human_clinical_association_supported = (
        human_clinical_context
        and host_virus_relationship_type
            not in negative_host_types
        and host_virus_passage_verified
        and host_virus_passage_mentions_host_virus
        and (
            combined_context_has_direct_language
            or has_human_clinical_context(
                combined_host_context
            )
            or has_human_clinical_context(
                evidence_source
            )
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
        verified_target_host_context
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

    # The relationship passage should identify the source virus when
    # possible. If it uses a pronoun/abbreviation, an explicitly extracted
    # comparison_source_virus_name matching host_virus_name can still link
    # the edge, provided the passage itself is verified and contains the
    # target comparison virus plus relationship language.
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
    # SAVE BACKWARDS-COMPATIBLE DIAGNOSTICS
    # ========================================================

    extraction["study_host_passage_verified"] = (
        study_host_passage_verified
    )
    extraction["study_host_matches_target"] = (
        study_host_matches_target
    )
    extraction["human_clinical_context"] = (
        human_clinical_context
    )
    extraction["verified_target_host_context"] = (
        verified_target_host_context
    )
    extraction["strong_study_host_context"] = (
        verified_target_host_context
    )
    extraction["host_virus_passage_verified"] = (
        host_virus_passage_verified
    )
    extraction["host_virus_matches_target"] = (
        host_virus_matches_target
    )
    extraction["comparison_virus_matches_target"] = (
        comparison_virus_matches_target
    )
    extraction["comparison_source_matches_host_virus"] = (
        comparison_source_matches_host_virus
    )
    extraction["comparison_relationship_passage_verified"] = (
        comparison_relationship_passage_verified
    )
    extraction["relationship_passage_verified"] = (
        comparison_relationship_passage_verified
    )
    extraction["comparison_relationship_mentions_target_virus"] = (
        comparison_relationship_mentions_target
    )
    extraction["relationship_mentions_target_virus"] = (
        comparison_relationship_mentions_target
    )
    extraction["comparison_relationship_mentions_source_virus"] = (
        comparison_relationship_mentions_source
    )
    extraction["host_virus_passage_has_direct_language"] = (
        host_virus_passage_has_direct_language
    )
    extraction["combined_context_has_direct_language"] = (
        combined_context_has_direct_language
    )
    extraction["host_specific_study_context"] = (
        host_specific_study_context
    )
    extraction["paper_has_host_specific_study_context"] = (
        paper_has_host_specific_study_context
    )
    extraction["direct_host_relationship_supported"] = (
        direct_host_relationship_supported
    )
    extraction["contextual_host_association_supported"] = (
        contextual_host_association_supported
    )
    extraction["human_clinical_association_supported"] = (
        human_clinical_association_supported
    )
    extraction["host_association_supported"] = (
        host_edge_verified
    )
    extraction["related_relationship_supported"] = (
        related_relationship_supported
    )
    extraction["comparison_edge_verified"] = (
        comparison_edge_verified
    )

    extraction["reported_virus_name"] = host_virus_name
    extraction["reported_host_or_sample"] = study_host
    extraction["reported_virus_matches_target"] = (
        host_virus_matches_target
    )
    extraction["reported_host_matches_target"] = (
        study_host_matches_target
    )
    extraction["reported_virus_passage_verified"] = (
        host_virus_passage_verified
    )

    # ========================================================
    # EXACT SUPPORT
    # ========================================================

    if (
        host_edge_verified
        and host_virus_matches_target
    ):
        extraction["relationship_type"] = (
            host_virus_relationship_type
        )
        extraction["relationship_passage"] = (
            host_virus_passage
        )
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
        extraction["relationship_type"] = (
            comparison_relationship_type
        )
        extraction["relationship_passage"] = (
            comparison_relationship_passage
        )
        extraction["classification_basis"] = (
            "VERIFIED_HOST_EDGE_PLUS_COMPARISON_EDGE"
        )
        return "TARGET_HOST_RELATED"

    # ========================================================
    # TARGET VIRUS EXPLICITLY ASSOCIATED WITH ANOTHER HOST
    # ========================================================
    # v0.4 deliberately does NOT classify a target virus as
    # VIRUS_OTHER_HOST merely because it appears as a comparison virus
    # known from another host. That is normal for phylogenetic/reference
    # comparisons. The target virus must itself be the host-associated
    # virus in a non-target study, or be explicitly labelled DIFFERENT_HOST
    # in a verified comparison passage.
    # ========================================================

    non_target_host_direct_edge = (
        not verified_target_host_context
        and host_virus_matches_target
        and host_virus_passage_verified
        and host_virus_passage_mentions_host_virus
        and direct_host_relationship_supported
    )

    explicit_comparison_other_host = (
        comparison_virus_matches_target
        and comparison_relationship_type == "DIFFERENT_HOST"
        and comparison_relationship_passage_verified
        and comparison_relationship_mentions_target
        and bool(comparison_host_or_sample)
    )

    if (
        non_target_host_direct_edge
        or explicit_comparison_other_host
    ):
        extraction["classification_basis"] = (
            "TARGET_VIRUS_VERIFIED_IN_NON_TARGET_HOST"
        )
        return "VIRUS_OTHER_HOST"

    # ========================================================
    # RELATED CLAIM WITHOUT A COMPLETE TWO-EDGE CHAIN
    # ========================================================

    if comparison_relationship_type in related_types:
        extraction["relationship_type"] = (
            comparison_relationship_type
        )
        extraction["relationship_passage"] = (
            comparison_relationship_passage
        )
        extraction["classification_basis"] = (
            "RELATED_LANGUAGE_BUT_INCOMPLETE_EDGE_CHAIN"
        )
        return "NO_SUPPORT"

    # ========================================================
    # EXPLICIT NON-SUPPORT / MENTION ONLY
    # ========================================================

    if host_virus_relationship_type == "MENTION_ONLY":
        extraction["relationship_type"] = "MENTION_ONLY"
        extraction["relationship_passage"] = host_virus_passage
        extraction["classification_basis"] = "MENTION_ONLY"
        return "MENTION_ONLY"

    if comparison_relationship_type == "REFERENCE_ONLY":
        extraction["relationship_type"] = "MENTION_ONLY"
        extraction["relationship_passage"] = (
            comparison_relationship_passage
        )
        extraction["classification_basis"] = "REFERENCE_ONLY"
        return "MENTION_ONLY"

    if (
        host_virus_relationship_type
        in {"DIFFERENT_HOST", "NO_RELATION"}
        or comparison_relationship_type
        in {"DIFFERENT_HOST", "NO_RELATION"}
    ):
        if (
            non_target_host_direct_edge
            or explicit_comparison_other_host
        ):
            extraction["classification_basis"] = (
                "EXPLICIT_NON_TARGET_HOST"
            )
            return "VIRUS_OTHER_HOST"

        extraction["classification_basis"] = (
            "EXPLICIT_NO_RELATION"
        )
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


def _extraction_needs_rescue(
    extraction,
    evidence_source,
    target_virus_names,
):
    """
    Rescue only genuine extraction failures.

    Trigger when the model produced no meaningful entity fields, or when
    the supplied evidence visibly contains a target-virus name/alias but
    neither host_virus_name nor comparison_virus_name was extracted.
    """

    extraction = extraction or {}

    core_fields = [
        extraction.get("study_host"),
        extraction.get("host_virus_name"),
        extraction.get("comparison_virus_name"),
    ]

    if not any(str(value or "").strip() for value in core_fields):
        return True

    target_visible = text_contains_any_name(
        evidence_source,
        target_virus_names,
    )

    extracted_virus = (
        str(extraction.get("host_virus_name") or "").strip()
        or str(extraction.get("comparison_virus_name") or "").strip()
    )

    return bool(target_visible and not extracted_virus)


def _run_extraction(messages, max_new_tokens=500):
    response = generate_text(
        messages,
        max_new_tokens=max_new_tokens,
        max_input_tokens=MAX_INPUT_TOKENS,
    )

    return _extract_json_object(response)


def _run_rescue_extraction(
    host,
    host_aliases,
    virus,
    virus_aliases,
    paper,
    abstract,
    host_context,
    virus_context,
):
    """
    Second-pass extraction for cases where the normal prompt returns an
    empty structure. It intentionally uses only focused evidence snippets.
    """

    alias_text = "\n".join(
        f"- {alias}"
        for alias in (host_aliases or [])
    )

    virus_alias_text = "\n".join(
        f"- {alias}"
        for alias in (virus_aliases or [])[:20]
    )

    rescue_system_prompt = """
You are rescuing a failed scientific evidence extraction.

Use ONLY the supplied title, abstract, HOST CONTEXT, and TARGET VIRUS
CONTEXT. Do not invent passages.

Your goal is to recover two explicit edges when supported:

1. study host -> host-associated virus
2. host-associated virus -> comparison virus

A host-specific virome/transcriptome/HTS study can establish the first
edge across two passages: one establishes the target host/study material
and another names a virus discovered/identified/characterized in that
study.

If a virus is described as discovered in the target host, use
DISCOVERED_IN. If the paper only establishes a host-specific discovery
study across passages, STUDY_ASSOCIATION is allowed.

For the second edge, preserve the virus on the LEFT side of a sequence or
phylogenetic comparison as comparison_source_virus_name and the virus on
the RIGHT side as comparison_virus_name.

The target virus may be written using a validated taxonomic synonym.
Do not treat a comparison virus as a virus of the target host unless the
paper actually supports that host association.

Return ONLY JSON with this schema:
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

    messages = [
        {
            "role": "system",
            "content": rescue_system_prompt,
        },
        {
            "role": "user",
            "content": rescue_user_prompt,
        },
    ]

    return _run_extraction(
        messages,
        max_new_tokens=500,
    )


def analyze_paper(
    host,
    host_aliases,
    virus,
    paper,
    biological_context=None
):
    full_text = paper.get("full_text", "") or ""
    abstract = paper.get("abstract", "") or ""
    abstract = abstract[:MAX_ABSTRACT_CHARS]

    # --------------------------------------------------------
    # Resolve target-virus names once per paper. The underlying resolver
    # is disk-cached, so repeated papers/tasks do not repeatedly hit NCBI.
    # Prefer the tax ID already supplied by the Biological Context Agent.
    # --------------------------------------------------------

    context_tax_id = (
        (biological_context or {})
        .get("target_virus", {})
        .get("tax_id")
    )

    virus_taxonomy = get_virus_taxonomy_context(
        virus,
        tax_id=context_tax_id,
    )

    virus_aliases = list(
        virus_taxonomy.get("aliases", [])
        or []
    )

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
        for alias in host_aliases
    )

    virus_aliases_text = "\n".join(
        f"- {alias}"
        for alias in virus_aliases[:20]
    )

    biological_context_text = (
        context_for_prompt(biological_context)
        if biological_context
        else "{}"
    )

    system_prompt = """
You are a scientific evidence extraction agent.

You are NOT responsible for the final classification.
Your job is to extract only what the paper explicitly supports.

A host and virus merely appearing in the same paper is NOT evidence
that the virus infects or is naturally associated with that host.
This is especially important for model organisms, cell lines,
experimental systems, controls, vectors, and comparison species.

BIOLOGICAL CONTEXT is a PRIOR only. It is not evidence.
VALIDATED TARGET-VIRUS ALIASES are entity-resolution hints only. They
may indicate historical/accepted names for the same taxon; they do not
by themselves establish a host relationship.

You must separately identify two possible relationships:

A. STUDY HOST -> HOST-ASSOCIATED VIRUS
B. HOST-ASSOCIATED VIRUS -> COMPARISON VIRUS

HOST-ASSOCIATED VIRUS can be supported by direct relationship evidence
or by a clearly host-specific virome/transcriptome/HTS discovery study.
In discovery studies, the host and virus may be established in separate
verified passages from the same study.

Valid host relationship types include:
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

Valid virus-to-virus comparison types include:
SEQUENCE_SIMILARITY
PHYLOGENETIC_COMPARISON
RELATED_VIRUS
SAME_FAMILY
SAME_GENUS
REFERENCE_ONLY
DIFFERENT_HOST
NO_RELATION
UNCLEAR

For relationship B, comparison_source_virus_name is the virus on the
left side of the comparison. In a valid TARGET_HOST_RELATED chain it
normally equals host_virus_name. comparison_virus_name is the compared
virus, which may be written using one of the validated target-virus
aliases.

Do not put a comparison virus into host_virus_name merely because it is
the requested TARGET VIRUS.

HUMAN CLINICAL NORMALIZATION:
If the target host is Homo sapiens, clearly human patients, clinical
specimens, nasopharyngeal/respiratory samples, or infected individuals
may be normalized to Homo sapiens. A human cell line by itself is not
natural-host evidence.

PASSAGES MUST BE VERBATIM substrings of the supplied PAPER evidence.
If a required passage does not exist, return an empty string. Never
invent a quote.

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
        {
            "role": "system",
            "content": system_prompt,
        },
        {
            "role": "user",
            "content": user_prompt,
        },
    ]

    print("\nEvidence Agent analyzing:")
    print(paper.get("title", ""))

    extraction_status = "SUCCESS"
    extraction_attempts = 1
    extraction_error = ""

    try:
        extraction = _run_extraction(messages)
    except Exception as error:
        extraction = _empty_extraction(
            f"Primary extraction failed: {error}"
        )
        extraction_error = str(error)

    target_virus_names = [virus]
    target_virus_names.extend(virus_aliases)

    needs_rescue = _extraction_needs_rescue(
        extraction,
        evidence_source,
        target_virus_names,
    )

    if needs_rescue:
        extraction_attempts = 2

        print(
            "Evidence extraction was empty/incomplete. "
            "Running focused rescue extraction."
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
            )

            if _extraction_needs_rescue(
                rescue,
                evidence_source,
                target_virus_names,
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

        except Exception as error:
            extraction_status = "FAILED"
            extraction_error = str(error)
            cleanup_gpu()

            if not extraction:
                extraction = _empty_extraction(
                    f"Rescue extraction failed: {error}"
                )

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
        supporting_quote = extraction.get(
            "host_virus_passage",
            "",
        )
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