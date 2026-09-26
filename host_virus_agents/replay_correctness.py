"""Read-only replay of preserved v0.6 evidence; no network or model calls.

Saved quotations are disjoint historical source fragments, not newly retrieved
papers. Unknown search/identity validation is intentionally not filled in.
"""
import argparse
from collections import Counter
from copy import deepcopy
import hashlib
import json
from pathlib import Path

from evidence_semantics import VERSION, classify, decide


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def replay_paper(paper, host, virus):
    original = paper["original_diagnostic"]
    h, c = original.get("host_edge", {}), original.get("comparison_edge", {})
    e = deepcopy(original)
    e.update(study_host=h.get("host", original.get("study_host", "")),
             host_virus_name=h.get("virus", original.get("host_virus_name", "")),
             study_host_passage=h.get("study_host_passage", ""),
             host_virus_passage=h.get("virus_passage", ""),
             host_virus_relationship_type=h.get("relationship_type", "UNCLEAR"),
             comparison_source_virus_name=c.get("source_virus", ""),
             comparison_virus_name=c.get("comparison_virus", ""),
             comparison_relationship_type=c.get("relationship_type", "UNCLEAR"),
             comparison_relationship_passage=c.get("relationship_passage", ""))
    sources = list(dict.fromkeys(text for text, verified in [
        (h.get("study_host_passage", ""), h.get("study_host_passage_verified")),
        (h.get("virus_passage", ""), h.get("virus_passage_verified")),
        (c.get("relationship_passage", ""), c.get("passage_verified")),
    ] if text and verified is True))
    classify(e, [host], virus, sources, virus_aliases=original.get("target_virus_aliases", []))
    e["replay_provenance"] = "Historical verified-quotation flags; original full source and host aliases unavailable."
    return e


def replay(path):
    raw = Path(path).read_bytes()
    fixture = json.loads(raw.decode("utf-8-sig"))
    cases = []
    for case in fixture["cases"]:
        papers = [replay_paper(p, case["host"], case["virus"]) for p in case["candidate_papers_considered"]]
        metadata = {"taxonomy_resolution": case["taxonomy_resolution"],
                    "taxonomy_resolved": None, "retrieval_complete": None,
                    "coverage_note": "Historical snapshots lack validated rank/name binding and per-source request outcomes."}
        result = decide(papers, metadata)
        cases.append(dict(host=case["host"], virus=case["virus"], expected_status=case["expected_status"],
                          before=case["predicted_status"], after=result["literature_status"],
                          classification=result["classification"], correct=result["literature_status"] == case["expected_status"],
                          confidence=result["confidence"], reason=result["reason"], papers=papers,
                          historical_taxonomy=case["taxonomy_resolution"],
                          source_search_status=case["search_status"]))
    assert Path(path).read_bytes() == raw, "Historical evidence changed during replay"
    return dict(mode="historical_quotation_replay_not_live", evidence_version=VERSION,
                fixture_sha256=hashlib.sha256(raw).hexdigest(), benchmark_sha256=sha("test_pairs.csv"),
                code_sha256={str(p): sha(p) for p in ["evidence_semantics.py", "evidence_agent.py", "judge_agent.py", "replay_correctness.py"]},
                score=sum(c["correct"] for c in cases), total=len(cases),
                before_score=sum(c["before"] == c["expected_status"] for c in cases),
                false_known=sum(c["after"] == "KNOWN" and c["expected_status"] != "KNOWN" for c in cases),
                false_novel=0, false_biologically_implausible=0,
                insufficient=sum(c["classification"] == "INSUFFICIENT_EVIDENCE" for c in cases), cases=cases)


def report(data):
    lines = ["# Correctness benchmark after implementation", "",
             "## Result and limits", "",
             f"Historical benchmark: **{data['before_score']}/{data['total']}**. Conservative preserved-quotation replay: **{data['score']}/{data['total']}** against unchanged legacy labels. **Live score: not measured. The 12/12 regression gate is not satisfied.**", "",
             "All 12 replay decisions abstain because the artifact does not record independently validated taxonomy/rank or complete search outcomes. This is a data-completeness limitation, not a demonstration of improved benchmark accuracy. The two historically correct KNOWN calls and four correct negative calls are not claimed preserved by this replay.", "",
             "The local bundled Python is 3.12.14. It lacks requests, torch and transformers; this workspace lacks models/Qwen3-8B and agent_env. No live retrieval or model inference was attempted. No model, dependency or paper downloads were performed. Historical model outputs were not repaired, relabeled or treated as fresh source documents.", "",
             "| Metric | Replay count |", "|---|---:|",
             f"| Correct legacy labels | {data['score']}/{data['total']} |",
             f"| False KNOWN | {data['false_known']} |", "| False NOVEL_CANDIDATE | 0 |",
             "| False BIOLOGICALLY_IMPLAUSIBLE | 0 |", f"| INSUFFICIENT_EVIDENCE | {data['insufficient']} |", "",
             "Zero false positive calls here reflects abstention, not solved sensitivity. The legacy benchmark has KNOWN, POSSIBLY_KNOWN and NO_EVIDENCE_FOUND labels. Neither POSSIBLY_KNOWN nor NO_EVIDENCE_FOUND is mapped to novelty; the canonical classification remains INSUFFICIENT_EVIDENCE without exact support.", "",
             "## All 12 before/after results", "",
             "| # | Host | Virus | Expected | Before | After (legacy / canonical) | Correct |",
             "|---|---|---|---|---|---|---|---|"]
    for i, c in enumerate(data["cases"], 1):
        lines.append(f"| {i} | {c['host']} | {c['virus']} | {c['expected_status']} | {c['before']} | {c['after']} / {c['classification']} | {c['correct']} |")
    lines += ["", "## Case evidence and reasoning", "",
              "The following records apply the new verifier to the original saved fields. Each source fragment is admitted only if its historical quotation-verification flag was true. Missing host aliases are not guessed. No comparison endpoints are swapped. Paper-level results are diagnostic re-verifications of saved quotations, not independent validations of the papers.", ""]
    for i, c in enumerate(data["cases"], 1):
        lines += [f"### {i}. {c['host']} → {c['virus']}", "",
                  f"Expected **{c['expected_status']}**; before **{c['before']}**; after **{c['classification']}** (legacy {c['after']}); correct: **{c['correct']}**.", "",
                  c["reason"] + " Complete per-source search coverage is also unavailable.", ""]
        for p in c["papers"]:
            edge = p["structured_evidence"]
            quote = (edge["supporting_text"] or p.get("host_virus_passage") or "No usable quotation saved.").replace("\n", " ")
            ident = p.get("pmid") or p.get("pmcid") or p.get("paper_id", "no ID")
            lines += [f"- **{p.get('paper_id', '')} / {ident}** — {p.get('title', '')}. `{p['classification']}`; {p['classification_basis']} Proposed edge: `{edge['host_name']} → {edge['relationship_type']} → {edge['virus_name']}`. Scope: `{edge['evidence_scope']}`. Saved quotation: {json.dumps(quote[:480], ensure_ascii=False)}"]
        lines.append("")
    lines += ["## Remaining work before accepting the regression gate", "",
              "1. Run fresh extraction on the unchanged 12 pairs in the project model environment. Preserve full source text, raw primary/rescue responses, taxonomy resolution and query outcomes in the new sidecar.",
              "2. Review the discovery paper's real specimen/study linkage and repair unresolved comparison roles through extraction, not endpoint substitution. The preserved introductory host description and virus genome-length fragments cannot establish that linkage.",
              "3. Recover complete direct assertions for the positive plant and mouse cases. Verify primary human infection evidence rather than the old background/cell-line passage. The conservative assertion grammar will abstain on unsupported prose; extend it only with generic positive and adversarial tests.",
              "4. Re-run failed extractions and all decoys. A negative benchmark label cannot be restored by manufacturing successful requests, suppressing an extraction error, or interpreting an unparsed assertion as non-support.",
              "5. Compare every previously correct case before acceptance. No scalability refactor or claim of stable 12/12 is justified yet.", "",
              "## Reproduction", "", "```bash", "python -B -m unittest discover -s tests -p 'test*correctness.py'", "python -B replay_correctness.py", "```", "",
              "On Nibi, from the project directory with its existing agent_env and model, run:", "", "```bash", "sbatch run_correctness_benchmark.sh", "```", "",
              "The runner executes both public benchmark paths and records results/correctness_<job>/pairs.csv, pairs.csv.diagnostics.jsonl, env.csv, env.csv.diagnostics.jsonl and run.log. It uses separate correctness cache files and never changes expected labels. See CORRECTNESS_IMPLEMENTATION_REPORT.md for implementation details and limitations.", "",
              "## Provenance", "", f"Fixture SHA-256: `{data['fixture_sha256']}`.", "", f"Unchanged test_pairs.csv SHA-256: `{data['benchmark_sha256']}`.", "",
              "Full replay decisions, original copied extraction fields and code hashes: results/benchmark_after_correctness_offline.json. The source fixture results/benchmark_diagnostic_v06.json remains unchanged."]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", default="results/benchmark_diagnostic_v06.json")
    args = parser.parse_args()
    data = replay(args.fixture)
    Path("results/benchmark_after_correctness_offline.json").write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    Path("BENCHMARK_AFTER_CORRECTNESS.md").write_text(report(data), encoding="utf-8")
    print(json.dumps({k: v for k, v in data.items() if k != "cases"}, indent=2))
