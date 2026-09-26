# Correctness iteration 2 — awaiting Nibi validation

Date: 2026-09-26. This iteration implements the user's requested correction of over-abstention. It does **not** complete the research application or begin the scalability refactor.

## Authoritative baseline and current result

The user supplied the authoritative summary of the latest live run: Qwen3-8B on Nibi, 12 completed episodes, **0/12 correct**, all predictions `UNCLEAR`, mean reward −1.250, mean steps 9, mean searches 3, mean papers analyzed 5, approximately 13.3 minutes. All 42 original offline tests passed there. The new live CSV/log are intentionally not in this workspace. Their detailed quotes, requests and model responses cannot be reconstructed from this summary.

The existing `results/bioresearch_env_v06_results.csv` is the older **6/12** artifact, SHA-256 `a8b0cee60f260d5cc676a21f732ff80baeb11facfe27203458acad203f8b8922`. It was not overwritten. The current implementation has **77 passing local Python tests**, preserving the original 42. The PowerShell historical-fixture integrity test also passes. **The post-change live score is pending. No 12/12 claim is made.**

The new historical quotation replay still abstains on all 12 rows because raw source coverage and validated identity metadata are missing. Its 0/12 legacy-label agreement is not a second live measurement and is not proof of improvement or regression against the supplied Nibi run.

## Implemented changes

### Materiality and uncertainty

`evidence_semantics.py` records a separate `evidence_state` and a versioned materiality assessment, with source digest, completeness flag and reason. States include `MATERIAL_UNRESOLVED`, `IRRELEVANT_OR_REJECTED`, `SUCCESSFUL_NO_SUPPORT`, `EXACT_SUPPORT`, `RELATED_SUPPORT`, `OTHER_HOST`, and `MENTION_ONLY`.

A failed model extraction is still recorded as `FAILED`. It can cease blocking a pair only through independent source inspection: for example, complete available body text has no resolved target-host mention outside virus names, or every target-virus unit is explicitly limited to a non-association context such as receptor expression/cell-line work. A title or absent name in a partial snippet is not sufficient exclusion. An unrecorded materiality decision remains unresolved.

The verifier also checks for a potentially exact local assertion outside the selected extraction. Such a statement prevents a background quotation or an other-host edge from clearing the entire paper. Related evidence cannot hide an unresolved candidate capable of providing exact support. Verified exact evidence still dominates unrelated failures. Successful recorded searches plus no support and no material uncertainty can yield legacy `NO_EVIDENCE_FOUND`; the canonical result remains `INSUFFICIENT_EVIDENCE`, never inferred novelty.

### Grounding and entity binding

Validation now uses the actual contiguous full text, or the available abstract when full text is absent, instead of overlapping prompt snippets. Literal sentence offsets are retained correctly, including leading whitespace. Sentence splitting preserves scientific genus initials and decimal numbers better than splitting at every period.

Bounded accepted constructions include direct detection/infection, tested-positive statements, explicitly governed population-detection lists, explicit virus-of-host definitions, experimental passaging, adjacent specimen anaphora, and explicit discovery naming. Separate host/virus quotations may link only when the original units are adjacent and their text establishes a specimen or naming connection. A generic virome paper followed by a genome-length statement remains insufficient.

Aliases continue to require normalized identity or a consistent explicit local definition. An extractor-reported acronym can now bind to its locally defined canonical virus. Comparison source acronyms are resolved independently. Ambiguous acronyms, numbered-virus substrings and undeclared shortened scientific names remain rejected. No benchmark entity, PMID or expected outcome was added to decision rules.

Environmental material, fecal/source-only evidence, cell-line work, proteins, exposure, negation, background and comparison language retain separate scope. Experimental infection/passaging is not labeled natural infection. Unsupported prose remains unresolved; this is a conservative assertion verifier, not a complete biomedical parser.

### Bounded rescue and structured aggregation

`evidence_agent.py` retains one primary attempt and at most one rescue. Rescue covers empty output, unbound material assertions, invalid comparison types, source mismatch and self-comparison. It does not auto-swap endpoints. Primary and rescue outputs remain separately auditable. A fully source-justified irrelevant candidate does not incur an unnecessary rescue, even when its first extraction failed.

Rescue instructions now request complete verbatim assertions and explicit local linkage rather than implying that Python will manufacture quotations. The output budget increased from 500 to 1,000 tokens to accommodate the structured schema and complete passages. Model-backed efficacy of this change remains to be tested on Nibi.

The judge, supporting-paper lists and strongest-evidence selection share structured support checks. A free-form `EXACT_SUPPORT` label alone cannot populate the exact-paper list. Exact target support is selected ahead of other-host evidence for the final quoted explanation. Materiality, relationship verification, scope and extraction status remain separate fields.

Confidence is LOW or MEDIUM with explicit factors for entity resolution, grounding, extraction completeness, provenance and search completeness. It is not a calibrated probability. Exact classification does not automatically cause HIGH confidence; incomplete provenance or search is visible without erasing valid positive evidence.

### Reproducible live capture

Both evaluators now write an append-only `.diagnostics.jsonl` sidecar after each completed input. They refuse to overwrite an existing result or sidecar. The sidecar includes actual query attempts/source outcomes, retrieved/selected papers, host/virus aliases, taxonomy resolution, original available source text, structured edges, scopes, materiality, extraction attempts, final reason and confidence. Pipeline exceptions become `UNCLEAR` / `INSUFFICIENT_EVIDENCE` with `processing_status=FAILED`; exceptions are retained rather than turned into empty success.

The manifest records input/source hashes, Git HEAD, Python/package versions, model configuration and weight digests, and cache configuration. `GENERATION_TRACE` lines in `run.log` record effective model input, token limits, deterministic generation settings and device/dtype. These are prerequisites for a reproducible run, not proof that the run will pass. Full trace output is intentionally verbose for this small correctness benchmark.

`run_correctness_benchmark.sh` runs both public benchmark paths in a new run directory, with separate cache files. `validate_live_benchmark.py` checks output alignment, score and required diagnostic presence. Its output explicitly requires manual evidence review; passing arithmetic is not biological acceptance.

## Tests and preserved safeguards

The suite covers direct and adjacent evidence, discovery naming, specimen linkage, explicit/ambiguous acronyms, aliases, other-host/other-virus isolation, comparison reversal through re-extraction, related chains, sequence similarity, environmental material, experiments versus natural infection, cell lines/receptor assays, immunization, VLP/protein/vector exposure, negation, taxonomy inventories, malformed/failed extraction, material versus irrelevant failure, incomplete search, zero results, unresolved taxonomy, duplicate input row retention, repeated-paper confidence, capture errors and overwrite protection.

Tests were added and observed failing before their corresponding behavior was implemented. The original 42 tests remain passing. `test_model.py` is now an explicit CLI smoke demonstration with lazy model imports and supplied host/virus arguments, so root-level offline discovery no longer loads a model. It is not counted as a model-backed test.

The historical PowerShell integrity test still enforces unchanged benchmark definitions, result records, original quotations, caches and overwrite protection. Historical **source** digests are now treated as provenance rather than an instruction that future source code must remain byte-identical; input/evidence hash checks were not relaxed.

## Leakage review

Scanning current top-level and `bioresearch_env` Python source for benchmark names/reference hints found only three old human-normalization references in `_legacy_classify_relationship_v06`, which is retained for comparison and has no active caller. Active decision and rescue paths contain no benchmark-specific host, virus or PMID exceptions. The old model demonstration's named benchmark example was replaced by CLI inputs. Benchmark labels and historical evidence are unchanged.

## Exact Nibi command and required return artifacts

From the Nibi project directory, using the existing `agent_env` and Qwen3-8B assets:

```bash
sbatch --account=def-acdoxey --gres=gpu:nvidia_h100_80gb_hbm3_3g.40gb:1 run_correctness_benchmark.sh
```

The resource request matches the existing project launcher; adjust allocation/resource flags if the allocation changes. The script does not navigate to a hard-coded project path. `MODEL_PATH` can select another explicitly staged model directory. No Windows or Codex process is required by this runner.

Return the **entire** `results/correctness_<job-id>_<UTC timestamp>/` directory, specifically:

- `run.log`, including tests and generation traces;
- `env.csv` and `env.csv.diagnostics.jsonl`;
- `pairs.csv` and `pairs.csv.diagnostics.jsonl`;
- `env_acceptance.json` and `pairs_acceptance.json`;
- `host_aliases.json`, `virus_taxonomy.json`, and `biological_context.json` when produced.

If the job terminates early, retain the partial directory and Slurm error output. Completed sidecar records remain available. This runner is a correctness diagnostic, not a restartable production batch engine.

## Remaining limits and next gate

The local runtime cannot execute this live workflow: its bundled Python lacks requests, Torch and Transformers, the Qwen3-8B model and project environment are absent, and Slurm/Bash are unavailable. The user's live summary is sufficient to guide this iteration but not to certify specific repaired paper edges.

The assertion grammar intentionally remains bounded. More natural prose, plural/common-name resolution, genuinely unregistered virus identity, multi-assertion papers and abstract-only materiality can still require re-extraction or adjudication. In particular, explicit taxonomy resolution remains a strong-decision gate; an unregistered publication-defined virus must not acquire a fabricated TaxID or a false resolution flag to satisfy the benchmark.

No new biological-implausibility rule was introduced. Broad clade discordance remains a prior. No novelty classifier or scalability implementation was introduced. All 12 live labels and evidence attribution must be audited after the next Nibi run. If any fail, use the saved source and model traces for the next general repair. **Correctness is not frozen; production/V2 acceptance is blocked on that external run.**
