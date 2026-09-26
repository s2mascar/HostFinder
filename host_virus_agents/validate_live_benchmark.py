"""Validate a captured live run; this does not replace manual evidence review."""
import argparse
import csv
import json
from pathlib import Path
from benchmark_capture import digest


def validate(sidecar, input_file="test_pairs.csv"):
    records = [json.loads(line) for line in Path(sidecar).read_text(encoding="utf-8").splitlines()]
    pairs = [r for r in records if r.get("kind") == "pair"]
    with open(input_file, encoding="utf-8-sig", newline="") as handle:
        expected = list(csv.DictReader(handle))
    problems, cases = [], []
    manifests = [r for r in records if r.get("kind") == "manifest"]
    if len(manifests) != 1 or manifests[0].get("input_sha256") != digest(input_file):
        problems.append("Missing or mismatched input manifest.")
    if len(pairs) != len(expected):
        problems.append("Output row count does not match the benchmark definition.")
    for i, (saved, definition) in enumerate(zip(pairs, expected), 1):
        row, details = saved["row"], saved["diagnostics"]
        if any(row.get(k) != definition[k] for k in ("host", "virus", "expected_status")):
            problems.append(f"Row {i} differs from the immutable benchmark input.")
        judge = details.get("info", {}).get("judge_diagnostics", details)
        evidence = judge.get("evidence_results", [])
        search = judge.get("search_metadata", details.get("search_metadata", {}))
        predicted = row.get("predicted_status")
        issues = []
        exact = [e for e in evidence if e.get("evidence_state") == "EXACT_SUPPORT"]
        related = [e for e in evidence if e.get("evidence_state") == "RELATED_SUPPORT"]
        unresolved = [e for e in evidence if e.get("evidence_state") == "MATERIAL_UNRESOLVED"]
        if not search.get("taxonomy_resolution") or not search.get("queries_attempted"):
            issues.append("Missing actual taxonomy/query diagnostics.")
        if predicted == "KNOWN" and not exact:
            issues.append("KNOWN lacks structured exact support.")
        if predicted == "POSSIBLY_KNOWN" and not related:
            issues.append("Related classification lacks a verified two-edge chain.")
        if predicted == "NO_EVIDENCE_FOUND" and (unresolved or not search.get("retrieval_complete")):
            issues.append("Negative call retains material uncertainty or incomplete retrieval.")
        if predicted == "NO_EVIDENCE_FOUND":
            analyzed = row.get("papers_analyzed", row.get("papers_examined", 0))
            if len(evidence) != int(analyzed or 0):
                issues.append("Negative call lacks all analyzed candidate records.")
            if any(e.get("extraction_status") == "FAILED" and
                   (e.get("evidence_state") != "IRRELEVANT_OR_REJECTED" or not e.get("materiality", {}).get("source_sha256"))
                   for e in evidence):
                issues.append("A failed extraction was not independently excluded.")
        for e in exact + related:
            if not e.get("verification_source") or not e.get("extraction_trace"):
                issues.append("Supporting evidence lacks source text or extraction trace.")
            if not e.get("structured_evidence", {}).get("supporting_text"):
                issues.append("Missing exact supporting span.")
        correct = predicted == definition["expected_status"]
        cases.append(dict(row=i, host=row["host"], virus=row["virus"], expected=definition["expected_status"],
                          predicted=predicted, correct=correct, issues=issues,
                          exact_count=len(exact), related_count=len(related), unresolved_count=len(unresolved)))
    checks_passed = len(pairs) == len(expected) and not problems and all(c["correct"] and not c["issues"] for c in cases)
    return dict(total=len(expected), completed=len(pairs), correct=sum(c["correct"] for c in cases),
                false_known=sum(c["predicted"] == "KNOWN" and c["expected"] != "KNOWN" for c in cases),
                false_related=sum(c["predicted"] == "POSSIBLY_KNOWN" and c["expected"] != "POSSIBLY_KNOWN" for c in cases),
                problems=problems, cases=cases, automated_checks_passed=checks_passed,
                manual_evidence_audit="REQUIRED", correctness_checkpoint="NOT_FROZEN",
                note="Checks verify arithmetic and diagnostic presence, not biological attribution. Inspect every accepted edge and every rejected material candidate before freezing correctness.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sidecar")
    parser.add_argument("--input", default="test_pairs.csv")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = validate(args.sidecar, args.input)
    with open(args.output, "x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    print(f"Captured score: {result['correct']}/{result['total']}; automated checks: {result['automated_checks_passed']}; manual evidence audit required.")
