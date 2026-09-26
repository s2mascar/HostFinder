"""Conservative, model-independent verification of directed literature assertions.

The extractor proposes endpoints; this module never repairs them by swapping
names. Unsupported syntax is uncertainty, not proof of absence. No taxonomy or
biological host-range rules are learned from benchmark labels.
"""
from dataclasses import asdict, dataclass
import re
import hashlib

VERSION = "directed-evidence-2"


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


def role_conflict(extraction, evidence_source=""):
    source = extraction.get("comparison_source_virus_name")
    other = extraction.get("comparison_virus_name")
    virus = extraction.get("host_virus_name")
    if source and other and entity_equal(source, other, local_names([other], evidence_source)):
        return "SELF_COMPARISON"
    if source and virus and not entity_equal(source, virus, local_names([virus], evidence_source)):
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


def sentence_spans(document):
    """Sentence/line spans with exact offsets, preserving genus initials."""
    starts = [0]
    ends = []
    for boundary in re.finditer(r"(?<=[.!?])\s+(?=[A-Z0-9])|\n+", document):
        ends.append(boundary.start())
        starts.append(boundary.end())
    ends.append(len(document))
    out = []
    for start, end in zip(starts, ends):
        raw = document[start:end]
        left, right = len(raw) - len(raw.lstrip()), len(raw.rstrip())
        if right > left:
            out.append((raw[left:right], start + left, start + right))
    return out


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
        start = document.find(quote)
        end = start + len(quote)
        out.extend(unit for unit in sentence_spans(document) if unit[1] < end and unit[2] > start)
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
        h + r"\s+(?:tested|was|were)\s+positive\s+for\s+" + v,
        r"\bWe\s+passaged\s+" + v + r"\s+in\s+" + h,
        v + r",\s*(?:the|a)\s+[a-z-]*virus of\s+(?:the\s+)?[a-z -]{1,30}\(\s*" + h + r"\s*\)",
    ]
    if any(re.search(p, text, re.I) for p in patterns):
        return True
    population = re.search(r"\bWe found (?:several|a number of) (?:viruses|microorganisms) in\s+" + h + r",\s*(?:including|the most prevalent of which were)\s+(.+?)[.]?$", text, re.I)
    if population:
        entries = re.split(r",\s*|\s+and\s+", population.group(1))
        return any(re.fullmatch(v + r"(?:\s*\(\d+(?:\.\d+)?%\))?", entry.strip().rstrip("."), re.I) for entry in entries)
    return False


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


def adjacent_link(units, host_names, virus_names, source):
    """Resolve only an explicit adjacent specimen or virus-naming anaphor."""
    if not isinstance(source, str):
        return None
    for first, second in zip(units, units[1:]):
        a, start, _ = first
        b, _, end = second
        gap = source[first[2]:second[1]]
        if gap.strip() or "\n\n" in gap:
            continue
        if re.search(r"\b(?:not|other|different|protein|antigen|hypothes\w*|might|could|previous|whether)\b", a + " " + b, re.I):
            continue
        for host in host_names:
            collected = re.fullmatch(r"We collected (?:samples|specimens|tissues) from\s+" + name_pattern(host) + r"\.", a, re.I)
            if collected and any(assertion(b, "these " + noun, virus) for noun in ("samples", "specimens", "tissues") for virus in virus_names):
                return a, b, start, end, "ADJACENT_SPECIMEN_LINK"
            if assertion(a, host, "A novel virus") and any(re.fullmatch(r"This virus was (?:named|designated)\s+" + name_pattern(v) + r"\.", b, re.I) for v in virus_names):
                return a, b, start, end, "EXPLICIT_DISCOVERY_NAMING"
    return None


def assess_materiality(e, source, hosts, viruses):
    """Independently reject an unresolved paper only on inspectable source data.

    Missing names in a snippet/abstract are not exclusion evidence. A failed
    model call is never itself a reason to call a candidate irrelevant.
    """
    status = e.get("classification")
    mapped = {"EXACT_SUPPORT": "EXACT_SUPPORT", "TARGET_HOST_RELATED": "RELATED_SUPPORT",
              "VIRUS_OTHER_HOST": "OTHER_HOST"}
    state, reason = mapped.get(status, "MATERIAL_UNRESOLVED"), e.get("classification_basis", "")
    documents = source if isinstance(source, list) else [source]
    complete = e.get("evidence_source_complete") is True and bool(any(documents))
    local_hosts, local_viruses = local_names(hosts, source), local_names(viruses, source)
    target_assertion_unprocessed = state != "EXACT_SUPPORT" and any(
        not re.search(r"\b(?:not|previous|whether|hypothes\w*|cell lines?|cell cultures?|protein|antigen|wastewater|feces)\b", unit, re.I)
        and any(assertion(unit, h, v) for h in local_hosts for v in local_viruses)
        for doc in documents for unit, _, _ in sentence_spans(doc))
    if target_assertion_unprocessed:
        state, reason = "MATERIAL_UNRESOLVED", "A potentially exact local assertion remains outside the accepted extracted edge."
    elif state == "MATERIAL_UNRESOLVED" and complete:
        names = local_names(viruses, source)
        # A host string embedded inside a virus's name is not an organism mention.
        host_source = "\n".join(documents)
        for name in sorted(names, key=len, reverse=True):
            host_source = re.sub(name_pattern(name), "[virus entity]", host_source, flags=re.I)
        if not any(mentions(host_source, h) for h in local_names(hosts, source)):
            state, reason = "IRRELEVANT_OR_REJECTED", "Complete available source contains no resolved target-host mention outside virus names."
        else:
            relevant = [s for doc in documents for s in re.split(r"(?<=[.!?])\s+|\n+", doc)
                        if any(mentions(s, v) for v in names)]
            excluded = r"\b(?:cell lines?|cell cultures?|receptor assays?|receptor expression|immuniz\w*|VLPs?|virus-like particles|vector-mediated|protein delivery|sequence similarity|taxonomic inventory|taxonomy list|background mention)\b"
            if relevant and all(re.search(excluded, s, re.I) and not any(assertion(s, h, v) for h in hosts for v in names) for s in relevant):
                state, reason = "IRRELEVANT_OR_REJECTED", "Every target-virus unit in complete available source is independently scoped to a non-association experiment/context."
                if status == "MENTION_ONLY" and e.get("extraction_status") != "FAILED":
                    state = "MENTION_ONLY"
    if state == "MATERIAL_UNRESOLVED" and not target_assertion_unprocessed and e.get("extraction_status") != "FAILED":
        edge = e.get("structured_evidence", {})
        if status in {"MENTION_ONLY", "NO_SUPPORT"} and edge.get("evidence_scope") not in {"AMBIGUOUS", "COMPARISON_CONTEXT_ONLY"}:
            state = "SUCCESSFUL_NO_SUPPORT"
        elif status == "MENTION_ONLY" and complete:
            # A model's MENTION_ONLY label alone cannot clear a relevant paper.
            state = "MATERIAL_UNRESOLVED"
    e["evidence_state"] = state
    e["materiality"] = dict(state=state, reason=reason, source_complete=complete,
                            source_sha256=hashlib.sha256("\n".join(documents).encode()).hexdigest(), version=VERSION)


def classify(extraction, host_aliases, target_virus, evidence_source,
             biological_context=None, virus_aliases=None):
    """Enrich the original extraction without altering its proposed endpoints."""
    e = extraction
    host = str(e.get("study_host") or "")
    virus = str(e.get("host_virus_name") or "")
    kind = str(e.get("host_virus_relationship_type") or "UNCLEAR").upper()
    edge = Evidence(host, virus, kind,
                    host_scope=str(e.get("host_scope") or "SPECIES").upper(),
                    virus_scope=str(e.get("virus_scope") or "SPECIES").upper())
    host_match = any(entity_equal(host, alias) for alias in local_names(host_aliases, evidence_source))
    target_names = local_names([target_virus, *(virus_aliases or [])], evidence_source)
    virus_match = entity_equal(virus, target_virus, target_names)
    edge.target_host_binding = "TARGET_HOST" if host_match else "OTHER_HOST" if host else "UNRESOLVED"
    edge.target_virus_binding = "TARGET_VIRUS" if virus_match else "OTHER_VIRUS" if virus else "UNRESOLVED"
    quote = str(e.get("host_virus_passage") or "")
    host_names = local_names([host, *host_aliases] if host_match else [host], evidence_source)
    virus_names = local_names([virus, target_virus, *(virus_aliases or [])] if virus_match else [virus], evidence_source)
    units = grounded_units(quote, evidence_source)
    status, reason = "UNCLEAR", "Missing, ungrounded or syntactically ambiguous directed assertion."
    conflict = role_conflict(e, evidence_source)
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
            edge.natural_vs_experimental = "EXPERIMENTAL" if re.search(r"\b(?:experimentally|passaged|inoculated|challenged|cell culture|cell line)\b", lower) else "NATURAL" if re.search(r"\b(?:naturally|wild|natural infection)\b", lower) else "UNSPECIFIED"
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
        local_units = sorted(set(units + grounded_units(str(e.get("study_host_passage") or ""), evidence_source)), key=lambda unit: unit[1])
        link = adjacent_link(local_units, host_names, virus_names, evidence_source)
        if status == "UNCLEAR" and link and edge.host_scope == edge.virus_scope == "SPECIES":
            a, b, start, end, method = link
            edge.verified, edge.linking_text = True, a
            edge.supporting_text, edge.source_start, edge.source_end = b, end - len(b), end
            edge.source_context, edge.evidence_scope = method, "DETECTION_OR_INFECTION"
            status = "EXACT_SUPPORT" if host_match and virus_match else "VIRUS_OTHER_HOST" if virus_match else "NO_SUPPORT"
            reason = "Explicit adjacent assertion linkage; no document-wide co-occurrence inference."
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
    if edge.verified and host_match and not virus_match and not conflict and entity_equal(comparison["comparison_virus"], target_virus, target_names):
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
    assess_materiality(e, evidence_source, host_aliases, [target_virus, *(virus_aliases or [])])
    return status


def retrieval_sufficient(metadata):
    planned = metadata.get("queries_planned", 0)
    return (metadata.get("retrieval_complete") is True
            and not metadata.get("source_failures")
            and metadata.get("source_successes", 0) >= 2
            and isinstance(planned, int) and planned > 0
            and metadata.get("queries_executed") == planned)


def is_exact_support(result):
    edge = result.get("structured_evidence", {})
    return (edge.get("verified") is True and edge.get("target_host_binding") == "TARGET_HOST"
            and edge.get("target_virus_binding") == "TARGET_VIRUS"
            and edge.get("host_scope") == edge.get("virus_scope") == "SPECIES"
            and edge.get("extraction_version") == VERSION and bool(edge.get("supporting_text"))
            and edge.get("evidence_scope") in {"DETECTION_OR_INFECTION", "NATURAL_INFECTION_DETECTION", "EXPERIMENTAL_INFECTION"}
            and result.get("extraction_status") != "FAILED")


def is_related_support(result):
    edge = result.get("structured_evidence", {})
    return (edge.get("verified") is True and edge.get("target_host_binding") == "TARGET_HOST"
            and edge.get("target_virus_binding") == "OTHER_VIRUS"
            and edge.get("extraction_version") == VERSION
            and result.get("comparison_edge", {}).get("verified") is True
            and result.get("extraction_status") != "FAILED")


def decide(results, metadata):
    """Return internal classification plus the existing benchmark vocabulary.

    NO_EVIDENCE_FOUND is a search outcome, never a novelty classification.
    """
    factors = dict(entity_resolution=metadata.get("taxonomy_resolved") is True,
                   search_complete=retrieval_sufficient(metadata),
                   extraction_complete=all(r.get("extraction_status") in {"SUCCESS", "RESCUED"} for r in results),
                   evidence_grounded=any(r.get("structured_evidence", {}).get("verified") is True for r in results),
                   provenance_available=any((r.get("pmid") or r.get("pmcid") or r.get("doi")) and r.get("verification_source") for r in results),
                   calibrated=False)
    insufficient = dict(classification="INSUFFICIENT_EVIDENCE", literature_status="UNCLEAR", confidence="LOW", confidence_factors=factors)
    if metadata.get("taxonomy_resolved") is not True:
        return dict(insufficient, reason="Host or virus taxonomy resolution is unresolved or unrecorded.")
    exact = [r for r in results if is_exact_support(r)]
    if exact:
        confidence = "MEDIUM" if factors["search_complete"] and factors["provenance_available"] and factors["extraction_complete"] else "LOW"
        return dict(classification="KNOWN", literature_status="KNOWN", confidence=confidence, confidence_factors=factors,
                    reason="Verified exact directed evidence supports this pair; experimental status is retained per edge.")
    related = [r for r in results if is_related_support(r)]
    material_pending = any(r.get("materiality", {}).get("state") == "MATERIAL_UNRESOLVED"
                           or (r.get("extraction_status") == "FAILED" and r.get("materiality", {}).get("state") != "IRRELEVANT_OR_REJECTED")
                           for r in results)
    if related and not material_pending:
        return dict(insufficient, literature_status="POSSIBLY_KNOWN",
                    reason="A verified related-virus chain provides context, not exact association evidence.")
    if not retrieval_sufficient(metadata):
        return dict(insufficient, reason="Search failed, was incomplete, or its coverage was not recorded.")
    def unresolved(r):
        if r.get("classification") == "EXACT_SUPPORT" and r not in exact:
            return True
        assessment = r.get("materiality", {})
        if assessment.get("version") == VERSION and assessment.get("source_sha256"):
            return assessment.get("state") == "MATERIAL_UNRESOLVED"
        return (r.get("classification") in {"UNCLEAR", "EXACT_SUPPORT"}
                or r.get("extraction_status") == "FAILED" or not r.get("structured_evidence"))
    if any(unresolved(r) for r in results):
        return dict(insufficient, reason="Candidate evidence contains unresolved extraction or attribution errors.")
    return dict(insufficient, literature_status="NO_EVIDENCE_FOUND",
                reason="Completed recorded queries found no verified exact association; novelty is not established.")
