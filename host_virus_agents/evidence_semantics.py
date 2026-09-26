"""Conservative, model-independent verification of directed literature assertions.

The extractor proposes endpoints; this module never repairs them by swapping
names. Unsupported syntax is uncertainty, not proof of absence. No taxonomy or
biological host-range rules are learned from benchmark labels.
"""
from dataclasses import asdict, dataclass
import re

VERSION = "directed-evidence-1"


def normalize(value):
    return " ".join(re.findall(r"[a-z0-9]+", str(value or "").lower()))


def entity_equal(reported, target, aliases=()):
    name = normalize(reported)
    return bool(name) and name in {normalize(x) for x in [target, *(aliases or [])] if x}


def name_pattern(name):
    tokens = normalize(name).split()
    return r"(?<!\w)" + r"[\W_]+".join(map(re.escape, tokens)) + r"(?!\w)" if tokens else r"(?!)"


def mentions(text, name):
    return bool(re.search(name_pattern(name), text or "", re.I))


def role_conflict(extraction):
    source = extraction.get("comparison_source_virus_name")
    other = extraction.get("comparison_virus_name")
    virus = extraction.get("host_virus_name")
    if source and other and entity_equal(source, other):
        return "SELF_COMPARISON"
    if source and virus and not entity_equal(source, virus):
        return "COMPARISON_SOURCE_MISMATCH"
    return ""


@dataclass
class Evidence:
    host_name: str
    virus_name: str
    relationship_type: str
    target_host_binding: str = "UNRESOLVED"
    target_virus_binding: str = "UNRESOLVED"
    host_scope: str = "SPECIES"
    virus_scope: str = "SPECIES"
    evidence_scope: str = "AMBIGUOUS"
    natural_vs_experimental: str = "UNSPECIFIED"
    source_context: str = "UNKNOWN"
    supporting_text: str = ""
    source_start: int | None = None
    source_end: int | None = None
    verified: bool = False
    extraction_version: str = VERSION
    verification_reason: str = "No directed assertion verified."
    linking_text: str = ""


def grounded_units(quote, source):
    """Expand a literal quotation only within its original paragraph/sentence.

    A list is a collection of disjoint historical snippets, never a document.
    No concatenation is allowed across snippet boundaries.
    """
    if not quote:
        return []
    sources = source if isinstance(source, list) else [source]
    if sum(doc.count(quote) for doc in sources) != 1:
        return []
    out = []
    for document in sources:
        if quote not in document:
            continue
        before = len(out)
        for match in re.finditer(r"[^\n.!?]+(?:[.!?]|$)", document):
            if quote.strip().rstrip(".!?") in match.group():
                out.append((match.group().strip(), match.start(), match.end()))
        if len(out) == before:
            # A multi-sentence quote is inspected sentence by sentence.
            start = document.find(quote)
            for match in re.finditer(r"[^\n.!?]+(?:[.!?]|$)", quote):
                out.append((match.group().strip(), start + match.start(), start + match.end()))
    return out


def assertion(text, host, virus):
    """Recognize bounded direct constructions, never bag-of-words co-mention."""
    h, v = name_pattern(host), name_pattern(virus)
    # Optional local acronym definitions are part of the entity, not a new edge.
    definition = r"(?:\s*\([A-Za-z][A-Za-z0-9-]{1,15}\))?"
    h, v = h + definition + r"(?!\s+(?:virus|protein|antigen)\b)", v + definition
    aux = r"\s+(?:(?:was|were|is|are|has been|have been)\s+)?"
    detection = r"(?:detected|identified|discovered|isolated|recovered|sequenced|found|obtained|amplified)"
    samples = r"(?:(?:the|a|an|wild|naturally infected)\s+)*(?:(?:samples|tissues|specimens|leaves)\s+(?:of|from)\s+)?"
    patterns = [
        v + aux + r"(?:naturally\s+)?" + detection + r"\s+(?:in|from)\s+" + samples + h,
        h + aux + r"(?:(?:naturally|experimentally)\s+)?infected\s+(?:with|by)\s+" + v,
        v + r"\s+(?:infects|infected|replicates in)\s+" + h,
        r"(?:detected|identified|isolated|found)\s+" + v + r"\s+(?:in|from)\s+" + samples + h,
    ]
    return any(re.search(p, text, re.I) for p in patterns)


def local_names(names, source):
    """Use explicit local definitions, not guessed acronyms or genus initials."""
    documents = source if isinstance(source, list) else [source]
    accepted = [n for n in names if n and not re.match(r"^[A-Za-z]\.\s", n)]
    for doc in documents:
        for name in list(accepted):
            for match in re.finditer(name_pattern(name) + r"\s*\(([A-Z][A-Za-z0-9-]{1,15})\)", doc, re.I):
                acronym = match.group(1)
                # Repeated definitions are ambiguous unless every definition
                # uses an accepted full name. Do not infer from initials alone.
                definitions = list(re.finditer(r"\(" + re.escape(acronym) + r"\)", doc))
                if all(any(re.search(name_pattern(n) + r"\s*$", doc[:d.start()], re.I) for n in accepted) for d in definitions):
                    accepted.append(acronym)
    return list(dict.fromkeys(accepted))


def classify(extraction, host_aliases, target_virus, evidence_source,
             biological_context=None, virus_aliases=None):
    """Enrich the original extraction without altering its proposed endpoints."""
    e = extraction
    host = str(e.get("study_host") or "")
    virus = str(e.get("host_virus_name") or "")
    kind = str(e.get("host_virus_relationship_type") or "UNCLEAR").upper()
    edge = Evidence(host, virus, kind,
                    host_scope=e.get("host_scope", "SPECIES").upper(),
                    virus_scope=e.get("virus_scope", "SPECIES").upper())
    host_match = any(entity_equal(host, alias) for alias in host_aliases)
    virus_match = entity_equal(virus, target_virus, virus_aliases)
    edge.target_host_binding = "TARGET_HOST" if host_match else "OTHER_HOST" if host else "UNRESOLVED"
    edge.target_virus_binding = "TARGET_VIRUS" if virus_match else "OTHER_VIRUS" if virus else "UNRESOLVED"
    quote = str(e.get("host_virus_passage") or "")
    host_names = local_names([host, *host_aliases] if host_match else [host], evidence_source)
    virus_names = local_names([virus, target_virus, *(virus_aliases or [])] if virus_match else [virus], evidence_source)
    units = grounded_units(quote, evidence_source)
    status, reason = "UNCLEAR", "Missing, ungrounded or syntactically ambiguous directed assertion."
    conflict = role_conflict(e)
    if e.get("extraction_status") == "FAILED":
        reason = "Extraction failed; absence of an extracted edge is not non-support."
    elif conflict:
        reason = conflict + "; targeted re-extraction required; endpoints were not swapped."
    elif kind in {"NO_RELATION", "MENTION_ONLY", "REFERENCE_ONLY", "SEQUENCE_SIMILARITY", "PHYLOGENETIC_COMPARISON", "SAME_GENUS", "SAME_FAMILY", "RELATED_VIRUS"}:
        status = "MENTION_ONLY" if kind == "MENTION_ONLY" else "NO_SUPPORT"
        edge.evidence_scope = "SEQUENCE_SIMILARITY_ONLY" if kind == "SEQUENCE_SIMILARITY" else "COMPARISON_CONTEXT_ONLY"
        reason = "Proposed relationship is contextual, not host-virus support."
    else:
        for text, start, end in units:
            lower = text.lower()
            if re.search(r"\b(?:whether|hypothes\w*|if|might|could|possibly)\b", lower) or text.endswith("?"):
                continue
            exclusion = None
            if re.search(r"\b(?:wastewater|sewage|environmental sample|sediment)\b", lower):
                exclusion = "ENVIRONMENTAL_ASSOCIATION"
            elif re.search(r"\b(?:feces|faeces|fecal|faecal|stool|gut contents|blood meal|ingested)\b", lower):
                exclusion = "SOURCE_MATERIAL_ONLY"
            elif re.search(r"\b(?:cell cultures?|cell lines?|cultured cells)\b", lower):
                exclusion = "EX_VIVO_CELL_EVIDENCE"
            elif re.search(r"\b(?:previous studies|previously reported|background|according to|et al)\b", lower):
                exclusion = "BACKGROUND_MENTION"
            elif re.search(r"\b(?:not|neither|without|absence of|no evidence)\b", lower):
                exclusion = "NEGATIVE_OR_ABSENT_INFECTION"
            elif re.search(name_pattern(virus) + r"\s+(?:protein|antigen|peptide|VLP|virus-like particle)\b", text, re.I):
                exclusion = "COMPONENT_OR_EXPOSURE_ONLY"
            if exclusion:
                edge.evidence_scope, status = exclusion, "NO_SUPPORT"
                reason = "The cited assertion is " + exclusion.lower() + "."
                continue
            if not any(assertion(text, h, v) for h in host_names if h for v in virus_names if v):
                continue
            edge.supporting_text, edge.source_start, edge.source_end = text, start, end
            edge.source_context = "DIRECT_ASSERTION"
            edge.natural_vs_experimental = "EXPERIMENTAL" if re.search(r"\b(?:experimentally|inoculated|challenged|cell culture|cell line)\b", lower) else "NATURAL" if re.search(r"\b(?:naturally|wild|natural infection)\b", lower) else "UNSPECIFIED"
            edge.evidence_scope = "EXPERIMENTAL_INFECTION" if edge.natural_vs_experimental == "EXPERIMENTAL" else "NATURAL_INFECTION_DETECTION" if edge.natural_vs_experimental == "NATURAL" else "DETECTION_OR_INFECTION"
            if edge.host_scope != "SPECIES" or edge.virus_scope != "SPECIES":
                status, reason = "NO_SUPPORT", "Higher-rank evidence cannot establish exact species support."
                break
            edge.verified = True
            status = "EXACT_SUPPORT" if host_match and virus_match else "VIRUS_OTHER_HOST" if virus_match else "NO_SUPPORT"
            reason = "Verified directed assertion: " + edge.target_host_binding + " -> " + edge.target_virus_binding + "."
            break
        # A cross-passage edge requires a literal named specimen and two
        # independently grounded assertions. General virome context is not a link.
        specimen = str(e.get("shared_specimen") or "")
        if status == "UNCLEAR" and len(specimen) >= 3 and edge.host_scope == edge.virus_scope == "SPECIES":
            host_units = grounded_units(str(e.get("study_host_passage") or ""), evidence_source)
            for host_text, _, _ in host_units:
                if re.search(r"\b(?:not|environmental|wastewater|previous|mixed species)\b", host_text, re.I):
                    continue
                linked = any(re.search(r"\bspecimen\s+" + name_pattern(specimen) + r"\s+was collected from\s+" + name_pattern(h) + r"\s*[.!]?$", host_text, re.I) for h in host_names)
                if not linked:
                    continue
                for text, start, end in units:
                    if re.search(r"\b(?:not|previous|protein|antigen|similar|environmental|wastewater)\b", text, re.I):
                        continue
                    if any(assertion(text, "specimen " + specimen, v) for v in virus_names):
                        edge.verified, edge.linking_text = True, host_text
                        edge.supporting_text, edge.source_start, edge.source_end = text, start, end
                        edge.source_context = "EXPLICIT_SPECIMEN_LINK"
                        edge.evidence_scope = "DETECTION_OR_INFECTION"
                        status = "EXACT_SUPPORT" if host_match and virus_match else "VIRUS_OTHER_HOST" if virus_match else "NO_SUPPORT"
                        reason = "Two grounded assertions explicitly bind the same named specimen."
                        break
    comparison = dict(verified=False, source_virus=e.get("comparison_source_virus_name"),
                      comparison_virus=e.get("comparison_virus_name"),
                      relationship_type=e.get("comparison_relationship_type"),
                      supporting_text=e.get("comparison_relationship_passage", ""))
    if edge.verified and host_match and not virus_match and not conflict and entity_equal(comparison["comparison_virus"], target_virus, virus_aliases):
        for text, _, _ in grounded_units(comparison["supporting_text"], evidence_source):
            if re.search(r"\b(?:not|no|neither)\b", text, re.I):
                continue
            left = local_names([virus], evidence_source)
            right = local_names([target_virus, *(virus_aliases or [])], evidence_source)
            if any(re.search(name_pattern(a) + r"\s+(?:(?:is|was)\s+)?(?:most\s+)?(?:similar|related|closest)\s+to\s+" + name_pattern(b), text, re.I) for a in left for b in right):
                comparison["verified"] = True
                status, reason = "TARGET_HOST_RELATED", "Verified host-to-other-virus edge and a separate virus comparison; no exact support."
                break
    e["comparison_edge"] = comparison
    edge.verification_reason = reason
    e["structured_evidence"] = asdict(edge)
    e["host_edge"] = dict(host=host, virus=virus, relationship_type=kind,
                          verified=edge.verified, target_host_match=host_match,
                          target_virus_match=virus_match, support_mode=edge.source_context,
                          host_virus_passage=edge.supporting_text or quote)
    e["classification_basis"] = reason
    e["classification"] = status
    return status


def retrieval_sufficient(metadata):
    planned = metadata.get("queries_planned", 0)
    return (metadata.get("retrieval_complete") is True
            and not metadata.get("source_failures")
            and metadata.get("source_successes", 0) >= 2
            and isinstance(planned, int) and planned > 0
            and metadata.get("queries_executed") == planned)


def decide(results, metadata):
    """Return internal classification plus the existing benchmark vocabulary.

    NO_EVIDENCE_FOUND is a search outcome, never a novelty classification.
    """
    insufficient = dict(classification="INSUFFICIENT_EVIDENCE", literature_status="UNCLEAR", confidence="LOW")
    if metadata.get("taxonomy_resolved") is not True:
        return dict(insufficient, reason="Host or virus taxonomy resolution is unresolved or unrecorded.")
    exact = [r for r in results if (r.get("structured_evidence", {}).get("verified") is True
             and r["structured_evidence"].get("target_host_binding") == "TARGET_HOST"
             and r["structured_evidence"].get("target_virus_binding") == "TARGET_VIRUS"
             and r["structured_evidence"].get("host_scope") == "SPECIES"
             and r["structured_evidence"].get("virus_scope") == "SPECIES"
             and r.get("extraction_status") != "FAILED")]
    if exact:
        return dict(classification="KNOWN", literature_status="KNOWN", confidence="MEDIUM",
                    reason="Verified exact directed evidence supports this pair; experimental status is retained per edge.")
    related = [r for r in results if r.get("structured_evidence", {}).get("verified") is True
               and r["structured_evidence"].get("target_host_binding") == "TARGET_HOST"
               and r["structured_evidence"].get("target_virus_binding") == "OTHER_VIRUS"
               and r.get("comparison_edge", {}).get("verified") is True]
    if related:
        return dict(insufficient, literature_status="POSSIBLY_KNOWN",
                    reason="A verified related-virus chain provides context, not exact association evidence.")
    if not retrieval_sufficient(metadata):
        return dict(insufficient, reason="Search failed, was incomplete, or its coverage was not recorded.")
    if any(r.get("classification") in {"UNCLEAR", "EXACT_SUPPORT"} or r.get("extraction_status") == "FAILED" or not r.get("structured_evidence") for r in results):
        return dict(insufficient, reason="Candidate evidence contains unresolved extraction or attribution errors.")
    return dict(insufficient, literature_status="NO_EVIDENCE_FOUND",
                reason="Completed recorded queries found no verified exact association; novelty is not established.")
