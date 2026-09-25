import re


# ============================================================
# RELATIONSHIP LANGUAGE
# ============================================================
# These patterns are used for paper ranking and deterministic
# verification. They deliberately separate:
#   1. direct host-virus relationship language,
#   2. host-specific discovery/virome context,
#   3. virus-virus comparison language.
# ============================================================

DIRECT_PATTERNS = [
    r"\binfect(?:s|ed|ing|ion|ions)?\b",
    r"\bisolat(?:e|ed|es|ing|ion|ions)\b",
    r"\bdetect(?:ed|ion|ing|s)?\b",
    r"\bidentif(?:y|ied|ies|ication)\b",
    r"\brecover(?:ed|ing|y)?\b",
    r"\bsequenc(?:e|ed|ing|es)\b",
    r"\bfound\s+in\b",
    r"\bpresent\s+in\b",
    r"\bobtained\s+from\b",
    r"\bcollected\s+from\b",
    r"\bderived\s+from\b",
    r"\bdiscover(?:ed|y)\b",
    r"\breport(?:ed|ing)\s+(?:in|from)\b",
    r"\bobserv(?:ed|e)\s+in\b",
    r"\bharbor(?:s|ed|ing)?\b",
    r"\bamplif(?:ied|ication)\s+from\b",
    r"\breplicat(?:e|ed|es|ing|ion)\s+in\b",
    r"\btransmi(?:t|ts|tted|tting|ssion)\b",
    r"\bnatural(?:ly)?\s+infected\b",
    r"\bassociated\s+with\b",
    r"\bhosted\s+by\b",
    r"\bhost\s+of\b",
]


# Host-specific study designs can provide valid evidence even when the
# host and virus are established across adjacent sentences rather than
# by a single "virus X was isolated from host Y" sentence.
STUDY_CONTEXT_PATTERNS = [
    r"\bvirom(?:e|ic|ics)\b",
    r"\bviral\s+diversity\b",
    r"\bvirus\s+diversity\b",
    r"\bvirus\s+discovery\b",
    r"\bnovel\s+(?:rna\s+|dna\s+)?viruses?\b",
    r"\bnew\s+(?:rna\s+|dna\s+)?viruses?\b",
    r"\btranscriptom(?:e|es|ic|ics)\b",
    r"\bmetatranscriptom(?:e|es|ic|ics)\b",
    r"\bmetagenom(?:e|es|ic|ics)\b",
    r"\brna[- ]?seq(?:uencing)?\b",
    r"\bhigh[- ]throughput\s+sequencing\b",
    r"\bnext[- ]generation\s+sequencing\b",
    r"\bsequencing\s+(?:of|from)\b",
    r"\bviral\s+sequences?\b",
    r"\bviruses?\s+(?:were\s+)?(?:identified|detected|discovered|characterized)\b",
    r"\b(?:identified|detected|discovered|characterized)\s+viruses?\b",
]


CLINICAL_HUMAN_PATTERNS = [
    r"\bhuman(?:s)?\b",
    r"\bpatient(?:s)?\b",
    r"\bclinical\s+(?:sample|samples|specimen|specimens)\b",
    r"\bnasopharyngeal\s+(?:swab|swabs|sample|samples)\b",
    r"\boropharyngeal\s+(?:swab|swabs|sample|samples)\b",
    r"\brespiratory\s+(?:sample|samples|specimen|specimens)\b",
    r"\bhuman\s+(?:subject|subjects|participant|participants)\b",
    r"\binfected\s+(?:person|persons|individual|individuals|patients)\b",
]


RELATED_PATTERNS = [
    r"\bsequence\s+similarity\b",
    r"\bsimilar\s+to\b",
    r"\bclosely\s+related\b",
    r"\bphylogen(?:y|etic|etically)\b",
    r"\bclosest\s+(?:relative|match)\b",
    r"\bhomolog(?:ous|y)\b",
    r"\bsame\s+(?:family|genus)\b",
    r"\bcluster(?:ed|s|ing)?\s+with\b",
    r"\bgroup(?:ed|s|ing)?\s+with\b",
]


def _matches_any(text, patterns):
    if not text:
        return False

    text = str(text).lower()

    return any(
        re.search(pattern, text)
        for pattern in patterns
    )


def has_direct_interaction_language(text):
    return _matches_any(
        text,
        DIRECT_PATTERNS,
    )


def has_host_specific_study_context(text):
    return _matches_any(
        text,
        STUDY_CONTEXT_PATTERNS,
    )


def has_human_clinical_context(text):
    return _matches_any(
        text,
        CLINICAL_HUMAN_PATTERNS,
    )


def has_related_language(text):
    return _matches_any(
        text,
        RELATED_PATTERNS,
    )


def interaction_proximity_score(
    text,
    target_host,
    host_aliases,
    virus,
    window=700,
):
    """
    Score local passages where the target virus, target host (or a
    validated host alias), and relationship language occur near one
    another.

    Exact scientific-name matches receive more weight than aliases.
    Host-specific virome/transcriptome context also receives a smaller
    bonus so discovery papers are not lost merely because the direct
    relationship is spread across multiple sentences.
    """

    if not text or not virus:
        return 0

    text_lower = str(text).lower()
    virus_lower = str(virus).lower()

    if virus_lower not in text_lower:
        return 0

    exact_host = (
        str(target_host).lower()
        if target_host
        else ""
    )

    aliases = [
        str(alias).lower()
        for alias in (host_aliases or [])
        if alias
    ]

    best_score = 0
    start = 0

    while True:
        position = text_lower.find(
            virus_lower,
            start,
        )

        if position == -1:
            break

        left = max(
            0,
            position - window,
        )

        right = min(
            len(text_lower),
            position + len(virus_lower) + window,
        )

        snippet = text_lower[left:right]

        exact_host_present = (
            bool(exact_host)
            and exact_host in snippet
        )

        alias_present = any(
            alias in snippet
            for alias in aliases
        )

        if exact_host_present:
            if has_direct_interaction_language(snippet):
                best_score = max(best_score, 20)
            elif has_host_specific_study_context(snippet):
                best_score = max(best_score, 14)
            elif has_related_language(snippet):
                best_score = max(best_score, 10)

        elif alias_present:
            if has_direct_interaction_language(snippet):
                best_score = max(best_score, 10)
            elif has_host_specific_study_context(snippet):
                best_score = max(best_score, 7)
            elif has_related_language(snippet):
                best_score = max(best_score, 5)

        start = position + len(virus_lower)

    return best_score


# ============================================================
# TYPE-SPECIFIC VALIDATION
# ============================================================

TYPE_PATTERNS = {
    "DETECTED_IN": [
        r"\bdetect(?:ed|ion|ing|s)?\b",
    ],
    "IDENTIFIED_IN": [
        r"\bidentif(?:y|ied|ies|ication)\b",
    ],
    "DISCOVERED_IN": [
        r"\bdiscover(?:ed|y)\b",
        r"\bnew\s+(?:rna\s+|dna\s+)?viruses?\b",
        r"\bnovel\s+(?:rna\s+|dna\s+)?viruses?\b",
    ],
    "ISOLATED_FROM": [
        r"\bisolat(?:e|ed|es|ing|ion|ions)\b",
    ],
    "RECOVERED_FROM": [
        r"\brecover(?:ed|ing|y)?\b",
    ],
    "SEQUENCED_FROM": [
        r"\bsequenc(?:e|ed|ing|es)\b",
        r"\btranscriptom(?:e|es|ic|ics)\b",
        r"\bmetatranscriptom(?:e|es|ic|ics)\b",
    ],
    "FOUND_IN": [
        r"\bfound\s+in\b",
    ],
    "PRESENT_IN": [
        r"\bpresent\s+in\b",
        r"\bobserv(?:ed|e)\s+in\b",
    ],
    "OBTAINED_FROM": [
        r"\bobtained\s+from\b",
        r"\bcollected\s+from\b",
        r"\bderived\s+from\b",
    ],
    "AMPLIFIED_FROM": [
        r"\bamplif(?:ied|ication)\b",
    ],
    "INFECTION_OF": [
        r"\binfect(?:s|ed|ing|ion|ions)?\b",
    ],
    "NATURAL_INFECTION": [
        r"\bnatural(?:ly)?\s+infect(?:ed|ion)?\b",
    ],
    "EXPERIMENTAL_INFECTION": [
        r"\bexperimental(?:ly)?\s+infect(?:ed|ion)?\b",
        r"\binoculat(?:ed|ion)\b",
        r"\bchallenge(?:d)?\s+with\b",
    ],
    "REPLICATES_IN": [
        r"\breplicat(?:e|ed|es|ing|ion)\b",
    ],
    "TRANSMITTED_TO": [
        r"\btransmi(?:t|ts|tted|tting|ssion)\b",
    ],
    "ASSOCIATED_WITH": [
        r"\bassociated\s+with\b",
        r"\bassociation\s+(?:with|between)\b",
        r"\bharbor(?:s|ed|ing)?\b",
    ],
    "STUDY_ASSOCIATION": STUDY_CONTEXT_PATTERNS,
}


RELATED_TYPE_PATTERNS = {
    "SEQUENCE_SIMILARITY": [
        r"\bsequence\s+similarity\b",
        r"\bsimilar(?:ity)?\s+to\b",
        r"\bidentity\s+(?:to|with)\b",
    ],
    "PHYLOGENETIC_COMPARISON": [
        r"\bphylogen(?:y|etic|etically)\b",
        r"\bcluster(?:ed|s|ing)?\s+with\b",
        r"\bgroup(?:ed|s|ing)?\s+with\b",
    ],
    "RELATED_VIRUS": [
        r"\bclosely\s+related\b",
        r"\brelated\s+to\b",
        r"\bclosest\s+(?:relative|match)\b",
        r"\bhomolog(?:ous|y)\b",
    ],
    "SAME_FAMILY": [
        r"\bsame\s+family\b",
        r"\bfamily\b",
    ],
    "SAME_GENUS": [
        r"\bsame\s+genus\b",
        r"\bgenus\b",
    ],
}


def relationship_type_supported_by_text(
    relationship_type,
    text,
):
    """
    Verify that a claimed direct relationship type has lexical support
    in the supplied evidence text.
    """

    if not text:
        return False

    relationship_type = str(
        relationship_type or ""
    ).upper().strip()

    patterns = TYPE_PATTERNS.get(
        relationship_type,
        [],
    )

    return _matches_any(
        text,
        patterns,
    )


def related_type_supported_by_text(
    relationship_type,
    text,
):
    if not text:
        return False

    relationship_type = str(
        relationship_type or ""
    ).upper().strip()

    patterns = RELATED_TYPE_PATTERNS.get(
        relationship_type,
        [],
    )

    if patterns:
        return _matches_any(
            text,
            patterns,
        )

    return has_related_language(
        text
    )
