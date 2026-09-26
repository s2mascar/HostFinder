# Host-virus pipeline architecture review

Audit date: 2026-09-25. Scope: the currently opened workspace only. This is a design review, not a V2 implementation or a validation of the underlying biological literature.

## 1. Executive summary

The prototype contains valuable retrieval, taxonomy, provenance, and relationship-separation ideas, but it is not ready to classify millions of predictions. The most serious problem is that a verified quotation is sometimes treated as a verified biological relationship. Those are different claims. A passage can be copied correctly while assigning its virus to the wrong host or confusing a comparison virus with the virus studied.

Three highest-priority correctness problems are:

1. **Incorrect entity-to-relationship binding.** Whole-paper discovery context can validate a comparison-only passage as a host edge. The saved v0.6 output actually makes this error, despite its stricter relationship-type checks.
2. **Failures and incomplete coverage become apparent absence of evidence.** Failed extractions can become `NO_SUPPORT`; two successful source calls suffice for `NO_EVIDENCE_FOUND`, regardless of unresolved taxonomy, failed sources, pagination, or candidate truncation. That label is not currently called novelty, but converting it to novelty would be unsafe.
3. **Entity identity and evidence scope are insufficiently constrained.** Virus substring matching, first-result taxonomy fallback, broad/common aliases, and unrestricted abbreviated names can merge distinct entities. Detection, natural infection, experimental infection, and source-organism context are not separate decision dimensions.

The three largest scaling problems are pair-specific research and inference, no persistent reusable paper/evidence index, and whole-file caches/output rewrites without transactional resume or distributed request coordination.

**Baseline recommendation:** use the **v0.4 recorded behavior as the regression reference**, because it has the strongest saved benchmark result, 8/12 with zero false `KNOWN` predictions against the supplied labels. Do not blindly copy `backups/v04`: four of its Python files are byte-identical to the v0.5 distribution. There is no authenticated, complete v0.4 executable snapshot established by this audit. Use the current modular code as the incremental migration scaffold, selectively preserving its improvements, while rebuilding a reproducible baseline and fixing evidence correctness before scaling. Neither v0.5 nor v0.6 should become the default merely because it is newer.

The proposed V2 direction is appropriate: normalize and resolve unique entities; evaluate inspectable plausibility; retrieve reusable entity-centered literature; match candidates; extract source-grounded relationships; persist evidence independently of predictions; classify by indexed joins and explicit coverage rules. Expensive work should scale primarily with unique entities, documents, and unresolved evidence units.

Audit method: read `AGENTS.md` and `PROJECT_SPEC.md` completely first; inspect all requested modules, the environment helpers, historical distributions, six result CSVs, and all six supplied cache files; compare historical files and cache hashes; aggregate saved CSV results without running the research pipeline. No benchmark labels or source files were changed. No live searches, model inference, benchmark reruns, dependency installations, or cache-generating imports were performed. Local model/data/environment directories are absent from the workspace inventory, and Python/Bash commands were not found by command discovery. Code findings below are static unless explicitly supported by saved diagnostics. Reported accuracy is historical label agreement, not a newly reproduced accuracy measurement.

## 2. Current architecture

There are two related entry paths:

```text
evaluate_pairs.py -> run_judge_agent
  -> run_evidence_agent
     -> biological context
     -> run_search_agent -> aliases -> query planner -> source calls
     -> merge/rank/top 8 -> analyze_paper for each paper
  -> deterministic judge -> optional model evidence summary

evaluate_env.py -> BaselineResearchAgent -> BioResearchEnv
  reset: aliases + biological context
  SEARCH actions -> pool/merge/rank/top 8
  ANALYZE_PAPER actions -> same analyze_paper
  same judge -> SUBMIT -> benchmark reward and metrics
```

The environment is a research-episode wrapper, not an independent biological classifier. It shares the current search, extraction, and judge implementations. The baseline follows a fixed plan, analyzes all selected papers, and submits once. It is not a learned agent.

`search_agent.py` also owns the tokenizer and Qwen3-8B model, loaded at import time. Consequently, importing biological decision code through the judge/evidence stack loads the model. Network, text processing, inference, decision rules, and command-line printing are coupled through module imports.

## 3. Current execution flow

| Question | Current behavior |
|---|---|
| Input | CSV evaluators materialize every row and require `host`, `virus`, and `expected_status`; direct agent CLIs accept two names. Optional input TaxIDs are not propagated. There is no production CSV/Parquet classification command. |
| Host normalization | Several helpers lowercase, strip non-ASCII-alphanumeric punctuation, and collapse whitespace. Input strings themselves are not replaced with a single canonical entity record. |
| Virus normalization | Similar string normalization; search broadening may remove a trailing number. Taxonomy aliases are introduced during paper analysis, not consistently during search planning/ranking. |
| Taxonomy | Host aliases, pair biological context, and virus aliases use separate NCBI lookup paths. An existing context TaxID is supplied to the virus resolver. |
| Synonyms | Host scientific/common/equivalent names and a generated genus abbreviation; virus scientific/common/equivalent names, acronyms, and additional synonyms. |
| Queries | Exact pair, preferred host alias, model proposals, then deterministic fallbacks; maximum three distinct query pairs. |
| Sources | PubMed, PMC, and Europe PMC, each with quoted host AND virus terms. |
| Retrieval | Up to five records per source per query in the main path. PubMed/PMC use ESearch followed by grouped EFetch; Europe PMC returns core records. |
| Deduplication | Pool all query results; select one preferred identifier per record; merge records with that key and keep longer text. |
| Ranking | Hand-weighted title/abstract/full-text names, co-occurrence and local relationship vocabulary; keep eight globally. |
| Context | NCBI lineage plus an optional local Virus-Host DB snapshot; deterministic broad-group prior, not LLM-generated biological knowledge. |
| Extraction | One primary model call per selected paper and pair, optionally a rescue call; JSON names/types/quotes; deterministic grounding/classification. |
| Host roles | Extracted study host plus host-alias/clinical-language checks on passages. There is no complete organism mention/coreference model. |
| Virus roles | Host-associated virus and comparison source/target fields, with string/alias matching and separate edge dictionaries. |
| Relationship types | Host association vocabulary and virus-comparison vocabulary; environmental/natural/experimental scope is not a complete independent schema. |
| Final label | Any `EXACT_SUPPORT` wins as `KNOWN`; then related support, ambiguity, and apparent absence, in that order. |
| Confidence | Fixed `HIGH`/`MEDIUM`/`LOW` from final label and source-call count, not calibrated probability. |
| Caching | Whole-file JSON entity aliases, virus taxonomy, and pair context; in-process Virus-Host DB indexes and per-episode analyzed-paper reuse. |
| API errors | Bounded retry wrappers; search failures are recorded, taxonomy failures often become cached fallback records. |
| No results | `NO_EVIDENCE_FOUND` if at least two source calls succeeded, otherwise `UNCLEAR`. No implemented `NOVEL_CANDIDATE` decision. |

## 4. Module-by-module responsibilities

| Module | Responsibility and audit assessment |
|---|---|
| `search_agent.py` | Model loading/generation, host query selection, LLM query planning, retry orchestration, ranking, top-eight selection. Useful global pooling; excessive coupling and per-pair inference. |
| `literature_search.py` | Current three-source adapters, XML/text parsing, identifier selection, source/query provenance merge. No persistent retrieval ledger or pagination. |
| `pubmed_search.py` | Older overlapping three-source CLI, despite its name. Current search imports `literature_search`, not this module. Older merger lacks current query-provenance handling. Retain as historical comparison, then retire from production. |
| `taxonomy_aliases.py` | NCBI requests/retries, host aliases, abbreviation generation, JSON cache, exact-name and bounded textual alias helpers. Does not expose resolution uncertainty in the public alias result. |
| `evidence_agent.py` | Snippet construction, virus aliases, extraction/rescue, quote grounding, two-edge verification, paper evidence classes, full evidence runner. Main correctness hotspot. |
| `judge_agent.py` | Deterministic precedence and reasons, count-based retrieval sufficiency/confidence, model strongest-evidence summary. Does not consume extraction failure or taxonomy status as decision gates. |
| `evaluate_pairs.py` | Legacy benchmark loop; status accuracy, timing, counts, repeated full CSV saves. Output excludes most extraction detail. |
| `evaluate_env.py` | Environment benchmark, prior/edge/rescue diagnostics, timing/reward, repeated CSV saves. Excludes `ERROR` rows from headline completed-case accuracy. |
| `bioresearch_env/env.py` | Episode phases/actions, search pooling, paper IDs `P0...`, per-episode analysis cache, deep-copied observations/history, terminal metrics. Not durable execution state. |
| `bioresearch_env/baseline_agent.py` | Executes existing query plan and all selected analyses, calls judge, selects supporting paper IDs. |
| `bioresearch_env/actions.py` | Action validation and old label vocabulary. |
| `bioresearch_env/rewards.py` | Search/analysis costs, correctness reward and evidence bonus based on the extractor's own classes. Not independent evidence verification. |
| `bioresearch_env/relationship_language.py` | Lexical features for ranking and claimed-type validation. Same features are too weak to prove endpoint binding. |
| `bioresearch_env/biological_context_agent.py` | NCBI lineage, Virus-Host DB indexes, environmental/source separation, approximate taxonomic distance, pair cache, safe prompt view. |
| `bioresearch_env/virus_taxonomy_agent.py` | Name/TaxID lookup, aliases and lineage, `resolved` flag, JSON cache. |
| `bioresearch_env/__init__.py` | Empty package marker. |
| `test_model.py` | Standalone Qwen3-0.6B inference demonstration, no assertions, different model from main pipeline. Not a correctness test suite. |
| `test_pairs.csv` | Twelve labeled cases, including reference hints and an unused `expected_evidence` column. |
| Benchmark shell scripts | Site-specific SLURM resource/environment setup, optional data download, syntax checks, evaluator launches. Names do not reliably identify the executed code version. |
| `requirements.txt` | Pinned environment export with many `+computecanada` builds. No portable dependency split, model revision manifest, or included runtime/data assets. |

## 5. Taxonomy and alias handling

Three resolvers duplicate normalization, NCBI XML parsing, and fallback behavior. Host alias lookup and biological-context lookup prefer an exact scientific name, then accept the first result. Virus alias lookup additionally checks synonyms, then also accepts the first result. Iteration over `.//Taxon` includes lineage nodes as well as top-level results, so candidate selection is not explicitly restricted to returned entities. A nonmatching first result can still be presented as resolved; adding the original query to its aliases further obscures the discrepancy.

The host alias API returns only strings. It discards TaxID, rank, match method, candidate ambiguity and failure status. It accepts `Includes` and common names without distinguishing their precision. Generated abbreviations such as `L. noctiluca` are excluded from preferred search aliases but remain eligible for evidence matching without document-local expansion. Common-name occurrence is not necessarily a species identifier.

`same_virus` (`evidence_agent.py:91`) accepts normalized equality or substring containment in either direction when the contained name is at least eight characters. Thus a shortened numbered-virus name can match a longer distinct target. A family/genus record is not explicitly equated by lineage, which is good, but names/ranks are not checked before calling a match exact. Search broadening and entity equivalence must become separate operations.

The supplied virus caches have eight entries each, seven resolved and one unresolved. Mouse hepatitis virus maps to Murine hepatitis virus; several records have rank `no rank`. The Lampyris partiti-like virus 1 entry is unresolved. These facts argue for retaining raw labels, TaxIDs, rank, taxonomy version, and resolution status separately; `no rank` is not automatically invalid, and a missing public TaxID must not silently become novelty.

V2 should use stable internal entity IDs with nullable authoritative TaxIDs, an explicit candidate-resolution table, and source-typed aliases. A verified document-local virus identifier can be retained provisionally; unresolved mapping should remain visible and normally block a strong classification until identity is adjudicated. Conflicting supplied TaxID/name pairs must not be silently repaired. Taxonomic ancestors and ambiguous aliases are retrieval hints, never exact identity by default.

## 6. Literature retrieval

The main path runs at most three searches across three sources: up to nine source calls and 45 raw returned records before deduplication. With nonempty PubMed/PMC responses, a source call contains both search and fetch, so three queries can entail 15 HTTP requests before retries and taxonomy. This is a structural upper bound for that plan, not an observed request count.

Exact quoted host/virus terms are used without query escaping or a reusable entity query representation. Model host terms must be in the alias set, but model virus terms have no equivalent validation. The planner always invokes the model; it does not use resolved virus aliases as a systematic retrieval plan. Removing trailing numbers is search expansion and cannot establish virus equivalence.

PubMed supplies abstracts, PMC supplies flattened body text, and Europe PMC supplies core title/abstract records. A returned PMCID does not trigger a general full-text acquisition stage. There is no systematic citation follow-up, supplementary-table ingestion, pagination, total-hit tracking, or corpus coverage ledger. All three sources overlap substantially; successful calls are not independent confirmations.

Global pooling before ranking is worth preserving. Ranking is a recall heuristic, not evidence, and no score cutoff proves absence. Top-five source limits and top-eight analysis can hide decisive papers. XML `itertext()` joined without separators can also concatenate adjacent text nodes; flattening loses section, table, and citation structure needed for interpretation.

`paper_key` selects PMID, otherwise PMCID, otherwise DOI, otherwise title. It cannot reconcile a PMID-bearing record and a PMCID-only record sharing the same PMCID because they enter under different keys. Missing-identifier empty titles can collapse unrelated records. V2 needs a canonical paper with multiple unique identifier mappings and preserved source-specific versions, not a single preferred-key dictionary.

## 7. Evidence extraction

The model sees title, up to 7,000 abstract characters, four host-centered snippets and six virus-centered snippets, using 900-character flanks and up to 25 occurrences per term. Total generation input is capped at 5,500 tokens. Prompt/context length can truncate evidence; the Python verifier checks its assembled evidence string rather than recording precisely which tokens the model actually received.

The extraction represents one host edge and one comparison edge per paper/pair. It cannot faithfully encode a paper with many hosts, viruses, experimental arms and comparisons. Extracting per requested pair also encourages target-conditioned attribution and repeats work when the same paper appears again.

Current safeguards to preserve include source-substring checking, separate host/comparison types, explicit edge fields, inferred-source diagnostics, selective rescue, and separate extraction status. However:

- `verified_target_host_context` establishes alias presence, not the grammatical host role. It can be true even when the extracted study-host string disagrees with the target.
- Direct support requires host, virus and relationship vocabulary in one passage, but does not resolve subject/object, negation, speculation, or which of several viruses the verb describes.
- `contextual_host_association_supported` allows whole-evidence discovery keywords to validate separate host/virus quotations without proving a shared sample or experimental arm.
- Clinical words such as patients or respiratory specimens are accepted as human cues without a robust species/cell-line/sample-role disambiguator.
- Comparison source linking can pass without that source occurring in the comparison passage; the supplied/extracted name or another verified passage can suffice.
- `DIFFERENT_HOST` comparison handling requires a nonempty other-host string, not independently resolved other-host evidence. Multi-host studies can therefore be reduced incorrectly to non-support.
- Successful JSON parsing is not successful extraction. Empty or role-invalid content can still reach `NO_SUPPORT`; failed extraction flags do not block the judge.

The v0.6 grounding helpers attempt to replace reconstructed quotes with exact source windows while keeping model roles/types. This can fix quotation formatting but cannot validate those roles. The new `_flexible_name_matches` (`evidence_agent.py:1236`) uses raw regex strings with double escapes: `r"(?<!\\w)"`, `r"[\\W_]+"`, and `r"(?!\\w)"`. The separator matches literal backslash/W/underscore rather than ordinary nonword separators. Static inspection identifies a defect; an analogous .NET regex diagnostic failed on an ordinary spaced/hyphenated name. A Python regression test remains necessary before repair. Stored grounding results do not prove the present file generated them.

V2 should extract all relevant assertions per document section/evidence unit, ground endpoints and predicates to spans, and resolve local abbreviations/coreference explicitly. Cross-passage assertions require a recorded linkage such as a sample ID, table row, local virus definition, or unambiguous experimental group. Merely replacing a quote with a nearby window must not convert an unsupported assertion into verified evidence.

## 8. Biological-context logic

The context builder loads a local Virus-Host DB TSV once per process and indexes virus names/TaxIDs. It tries exact virus name, then TaxID, then a whole-database family scan. Known hosts are deduplicated, target-host records removed from the public prior, and remaining hosts compared by shared lineage names and a rank-distance table. Broad groups include plants, animals, fungi, bacteria, archaea and other eukaryotes.

Environmental/root placeholders are excluded from biological host counts. Numeric source-organism IDs are resolved separately and retained as environmental context. This distinction is valuable: the v0.3-to-v0.4 cache difference adds source-organism lineage information, including Odonata and Araneae, without turning them into true hosts.

Each biological cache has twelve pair entries. v0.4, v0.5 and v0.6 biological cache files are byte-identical; v0.5 and v0.6 virus caches are also byte-identical. The pair priors show why broad-group agreement is weak: Aedes aegypti/Mouse hepatitis virus is `BIOLOGICALLY_CONSISTENT` because both target and recorded host fall under ANIMAL. This is not evidence for the association or a demonstrated compatibility rule.

Exact target-host database records are hidden from the model using a safe prompt view; preserve that isolation in benchmark mode. It should not automatically become the production knowledge policy. In production, curated records can be evidence if imported with provenance and reviewed evidence criteria; they must be distinguished from literature independently retrieved by the pipeline. Database-included and database-held-out evaluation should be separate tracks.

The current judge does not directly evaluate the prior. Context influences extraction prompts, which introduces anchoring risk even when described as a prior. V2 plausibility should use versioned, inspectable rules with explicit uncertainty, not treat absence in known host ranges, a broad-group mismatch, or family fallback alone as strong implausibility. Strong evidence and strong plausibility conflicts require review rather than an automatic override.

## 9. Judge/classification logic

`determine_final_status` applies this precedence: `EXACT_SUPPORT` -> `KNOWN`; otherwise `TARGET_HOST_RELATED` -> `POSSIBLY_KNOWN`; otherwise any `UNCLEAR` -> `UNCLEAR`; otherwise empty/all-nonsupport results -> `NO_EVIDENCE_FOUND` when `source_successes >= 2`, else `UNCLEAR`.

`retrieval_sufficient` (`judge_agent.py:19`) ignores source identity, planned/executed queries, failures, truncation, taxonomy, and extraction health. Two calls to the same database on different queries qualify. `retrieval_complete` only means no recorded source-call failures; it is not completeness of the literature. No final `BIOLOGICALLY_IMPLAUSIBLE` or `NOVEL_CANDIDATE` logic exists.

Confidence is LOW for unclear/insufficient calls, HIGH for known after the count threshold, otherwise MEDIUM. It does not quantify identity certainty, assay quality, independent studies, or calibration. The model only summarizes strongest evidence after the deterministic decision, but its summary is not subsequently quote-verified. The generic unclear reason misleadingly claims a verified host/virus passage even when the cause is total retrieval failure.

Recommended V2 decision policy:

| Classification | Required conditions |
|---|---|
| `KNOWN` | At least one credible accepted assertion supports the exact resolved endpoints under a documented association-scope policy; evidence provenance is available. Missing unrelated searches do not erase valid positive evidence, but coverage remains reported. |
| `BIOLOGICALLY_IMPLAUSIBLE` | Strong, documented biological/taxonomic conflict with adequate identity resolution and no credible unresolved contradiction. Missing host-range records alone do not qualify. |
| `NOVEL_CANDIDATE` | Adequate, versioned search coverage as of a stated date; relevant evidence units processed successfully; no credible exact support; sufficient entity resolution; positive plausible assessment. Unknown plausibility is insufficient. |
| `INSUFFICIENT_EVIDENCE` | Unresolved identity, relevant retrieval/extraction failure, incomplete search budget, inaccessible decisive content, ambiguous assertions, unknown plausibility, or conflicting evidence prevents a stronger decision. |

Keep relationship-level relatedness as a separate result field. `POSSIBLY_KNOWN` is not an exact association class and has no safe automatic mapping to `KNOWN` or novelty. Likewise, legacy `NO_EVIDENCE_FOUND` must not be mass-mapped to novelty. Reclassify against V2 gates while retaining legacy labels.

Define coverage operationally: named source/corpus snapshots, approved alias and entity search plans, pagination completion or an audited stopping policy, dated retrieval, relevant content availability, candidate-processing completion, and no unresolved material failure. Zero hits can participate in an adequate coverage record; zero hits alone cannot establish adequacy. Novelty remains bounded by that search policy and date.

## 10. Current caching and persistence

| Cache/state | Key, reuse, and limitation |
|---|---|
| `host_alias_cache.json` | Normalized host -> strings. Referenced by code but absent in this workspace. Entire file read each call, rewritten on misses; API fallback `[host]` cached without expiry/status. |
| Biological context JSON | Normalized host `||` virus -> full context. Pair-keyed, repeating taxonomy/range data across combinations. No dataset digest, TTL, transaction or concurrent writer protection. |
| Virus taxonomy JSON | `taxid:...` or `name:...` -> record. Useful entity reuse, but equivalent keys are not centrally unified; failures remain cached until refresh. Entire file read even during repeated paper analysis. |
| Virus-Host DB indexes | Full TSV and name/TaxID maps per process; family fallback scans rows. Source taxonomy cache is process-local. |
| Environment analysis dictionary | Paper IDs only within an episode; repeated actions can reuse an analysis. No cross-pair or cross-run evidence reuse. |
| Result CSV | Saved after each completed pair/episode by truncating and rewriting all accumulated rows. Partial persistence, not resumability or atomic checkpointing. |

Corrupt/unreadable JSON is generally treated as an empty cache, which loses the diagnostic distinction and may trigger expensive rework. Direct whole-file writes can be interrupted or overwrite another worker's updates. No persistent paper bodies, query results, extraction units, cache access metrics, or stage-level work ledger exists.

## 11. Current benchmark design

The twelve cases comprise four `KNOWN`, four `POSSIBLY_KNOWN`, and four `NO_EVIDENCE_FOUND` labels. Five cases share the Lampyris discovery-paper hint `PMC7093385`; the negative Tribolium case also retrieves that paper in saved diagnostics. Repeated hosts and viruses exist, but no efficiency assertions test reuse. There are no duplicate input rows and no curated true-novel, failure, or unresolved-taxonomy expected outcomes.

Evaluators score only `expected_status`; `expected_evidence` is not asserted. Its `RELATED_ONLY` term differs from current `TARGET_HOST_RELATED`, further showing that it is descriptive rather than an enforced contract. Reference hints are copied, not used to validate retrieved evidence. The environment hides labels from observations; the baseline does not consume expected labels, although the environment object stores them for scoring.

The benchmark can reveal gross label regressions and some comparison-virus/decoy false positives. Per-paper diagnostics in later CSVs help localize errors. It cannot establish general biological accuracy, evidence-span correctness, natural-host validity, novelty search adequacy, scalability, or calibrated confidence. Four related cases from one paper are correlated, not four independent validations.

Every saved run analyzed 60 paper/pair combinations. `papers_retrieved` in the environment output is the retained candidate count, not total hits or all unique retrieved records. The current environment evaluator's completed-case accuracy excludes `ERROR`; report all-input accuracy and failure rate alongside completed-case accuracy. Its reward bonus trusts extractor classes, so a wrong `EXACT_SUPPORT` can earn an evidence bonus.

## 12. Historical version comparison

The following figures are recomputed from the six supplied CSVs using their existing labels. Runtime is the sum of recorded episode seconds, excluding any startup outside episode timing. No causal runtime comparison is established without pinned models, hardware, responses and cache states.

| Saved result | Correct | Label agreement | False `KNOWN` | Episode seconds | Main observation |
|---|---:|---:|---:|---:|---|
| `bioresearch_env_results.csv` | 7/12 | 58.3% | 3 | 439.52 | All four positives recovered; two decoys and one comparison virus promoted to known. |
| `bioresearch_env_v02_results.csv` | 6/12 | 50.0% | 0 | 544.73 | Decoys corrected; Lampyris and human positives lost; all related cases missed. |
| `bioresearch_env_v03_results.csv` | 7/12 | 58.3% | 0 | 539.80 | Human positive restored; Lampyris and all related cases still missed. |
| `bioresearch_env_v04_results.csv` | 8/12 | 66.7% | 0 | 597.98 | Strongest saved agreement; recovers related Hubei partiti-like virus 51 case. |
| `bioresearch_env_v05_results.csv` | 6/12 | 50.0% | 5 | 647.13 | All four related cases and Tribolium negative falsely called known; Lampyris positive still missed. |
| `bioresearch_env_v06_results.csv` | 6/12 | 50.0% | 1 | 663.11 | Lampyris positive recovered and most false-known calls removed; Arabidopsis/Turnip and mouse/hepatitis positives lost. |

Feature history supported by available code/artifacts:

- Early repository history introduces retrieval/taxonomy, agentic search, structured extraction, judge and benchmark modules. Complete independently pinned v0.1-v0.3 environments are not present in the named historical folders.
- Current search and environment pool all planned searches before selection and favor scientific host names and local relationship language. These are worthwhile retrieval safeguards; available artifacts do not establish exactly which early run first used every feature.
- The v0.3 cache already separates environmental records; the v0.4 cache adds source-organism taxonomy fields.
- The v0.5 change note identifies focused rescue, `DISCOVERED_IN`, virus taxonomy aliases, alias-aware snippets/matching, and extraction diagnostics. It says search, ranking, rewards, prior logic, and judge were intentionally unchanged.
- v0.6 adds stricter host/comparison type separation, grounded host checks, normalized relationship labels, role-conflict rescue, source-window replacement, and more grounding diagnostics. These improve inspectability but do not eliminate role errors.

Historical integrity findings: `backups/v04/evidence_agent.py`, `evaluate_env.py`, `relationship_language.py`, and `biological_context_agent.py` are byte-identical to the corresponding v0.5 files. The historical relationship-language file is identical to current. Historical/current biological-context and virus-taxonomy modules differ only in default cache filenames. Historical/current test-pair CSVs are identical. Therefore directory labels and cache suffixes do not authenticate executable versions.

The v0.5 shell script contains a malformed bare path with a trailing quote where directory setup should be. The backup script uses script-relative directory setup; the current script uses a fixed HPC path. `run_benchmark_v2.sh` actually targets v0.4-named outputs while executing current modules; it is not a V2 implementation or a historical-version selector. The current launcher still has a v05 job name with v06 outputs. Preserve these artifacts for comparison, not as authoritative version manifests.

## 13. Regressions identified

**R1: v0.5 comparison types can become host evidence.** In the historical `classify_relationship`, contextual support only excludes a small negative-type set. A `SEQUENCE_SIMILARITY` host type can pass through `STUDY_CONTEXT`, become a verified host edge, and produce exact support. The saved v0.5 Lampyris/Hubei 31 and Tribolium/Hubei 31 diagnostics demonstrate this exact path. For Tribolium, the study-host quote concerns other organisms while the model-provided host name supplies the target match.

**R2: v0.6 narrows types but preserves the contextual attribution hole.** Saved v0.6 Lampyris/Hubei 31 (`P0`, PMID 31900852) is `EXACT_SUPPORT`, `RESCUED`, `SEQUENCED_FROM`, `STUDY_CONTEXT`. Its virus quote describes a BLAST similarity to the target, and its separate comparison passage places that target in spider mix. A permitted host-type label now disguises the same attribution error. The type gate alone cannot solve it.

**R3: v0.6 loses direct positives through under-sized quotes and literal type validation.** Arabidopsis/Turnip has seven `NO_SUPPORT` and one `VIRUS_OTHER_HOST` results. `P0` keeps a host quote and a separate virus-effect quote but has support mode NONE. The grounding pass leaves an exact but incomplete quote alone; stricter direct rules then reject it. Mouse/hepatitis examples include an abbreviation-only virus quote (`MHV3`) and a prevalence sentence labeled `DETECTED_IN` that lacks the literal detection pattern. These support an extraction/verification-boundary explanation; they do not prove retrieval recall was unchanged.

**R4: extraction failures are classified as biological non-support.** All four `FAILED` extractions in saved v0.6 diagnostics become `NO_SUPPORT`: mouse/hepatitis `P1`, Lampyris/chuvirus `P2` and `P4`, and Arabidopsis/hepatitis `P0`. Parse/rescue failure is not negative evidence. The judge never examines `extraction_status`.

**R5: related-edge endpoint reversal remains.** In v0.6 Lampyris/Hubei 51, a verified host edge names Lampyris partiti-like virus 2, but the comparison source is incorrectly Hubei 51 itself. The source mismatch defeats related support. This is a role-extraction error, not a reason to relax exact matching.

**R6: new grounding implementation has escaped-regex defect.** See section 7. Treat this as a current code defect, not a proven explanation for every saved regression: saved outputs and source files lack an immutable run association.

**R7: reproducibility regressed through mislabeled snapshots and launchers.** A supposed v04 backup contains v05 code; a launcher called v2 writes v04 outputs from current source. Newer filenames cannot support controlled version comparisons.

## 14. Recommended V2 baseline

Adopt v0.4's **recorded behavior** as the initial regression comparison, not as a verified production classifier. It outperforms both later results and avoids their benchmark false-known errors, but still misses a known interaction and three related cases. Its 8/12 score is a floor for legacy-label parity, not the acceptance criterion for biological correctness.

Use the current module boundaries as the migration scaffold, preserving global search pooling, taxonomy synonyms, environmental separation, structured edge diagnostics, and deterministic final decisions. Reimplement the unsafe verification rules behind typed interfaces after tests exist. Keep v0.5 rescue/alias functionality and v0.6 type separation selectively, with explicit regression tests for each.

Before claiming a runnable v0.4 baseline, recover or reconstruct it from authenticated source/run evidence and replay frozen inputs. If that is impossible, name the reference accurately as “v0.4 recorded results”; build a new reproducible baseline from available code and fixtures without falsely attributing it to v0.4. Never substitute the v05-identical backup silently.

## 15. Correctness failure modes

| Failure mode | Existing exposure | V2 requirement |
|---|---|---|
| Co-mention | Local name/verb presence can satisfy support. | Bind both endpoints to the same assertion; co-mention is a retrieval feature only. |
| Other host -> target virus | Model target-host name, generic context, multi-host passages, and source-window replacement can misassign host. | Resolve actual host mention, sample and experimental arm; store the other-host assertion independently. |
| Target host -> other virus | Comparison virus can be placed on host edge. | Verify virus role; separate detected entity from comparison target. |
| Related/comparison virus | Same-family/genus language and comparison keywords can influence host inference. | Store virus-virus edges separately; never traverse relatedness to assert exact host association. |
| Environmental/source organism | Context module separates placeholders, but paper extraction lacks equivalent complete role typing. | Separate sample source, biological host, vector, prey/diet, symbiont and contamination hypotheses. |
| Experimental versus natural | Both types yield exact support without separate scope. | Persist setting and assay; experimental evidence must never be reported as natural infection. |
| Sequence similarity | Whole-paper virome context can promote comparison text. | Similarity alone supports only a comparison assertion. |
| Background and citations | Body/abstract flattening and keywords erase rhetorical/source roles. | Mark primary finding, background, cited claim, reference-only and speculation; follow cited evidence where needed. |
| Genus/family -> species | String identity and rank-free matching permit overbroad attribution. | Exact endpoints/rank policy; no downward propagation from ancestors. |
| Broad aliases/abbreviations | Alias strings carry no ambiguity scope. | Source/type/precision and document-local disambiguation. |
| Negation/hypothesis | Lexical matching ignores “not detected,” uncertain claims and negative controls. | Assertion polarity/modality with span-grounded interpretation. |
| Short or invalid quotes | Literal verification loses real evidence; automatic widening can introduce unrelated evidence. | Retrieve coherent sentence/table units and validate their relational meaning. |
| Failed or partial extraction | Empty fields become non-support. | Technical failure and biological non-support are disjoint states. |
| No hits/incomplete taxonomy | Count-based absence verdict. | Coverage and resolution gates; default to insufficient when material prerequisites fail. |

Biological plausibility rules must be justified independently of benchmark outcomes. Do not add a rule for a particular named benchmark pair or assume cross-group rarity proves impossibility.

## 16. Scalability bottlenecks

1. Every pair repeats a model query-planning call, up to nine source calls, up to eight primary extraction calls plus rescues, and usually a summary call. The maximum selected-paper path is 18 model generations per pair with eight rescues. Identical viruses and documents do not amortize this work.
2. Pair-keyed biological context repeats host/virus taxonomy work across new combinations. Host aliases and source taxonomy have some reuse, but resolver implementations and cache identities are fragmented.
3. Papers are downloaded before global deduplication and are not persisted; the same paper is fetched and extracted in different queries, pairs and runs.
4. Whole JSON caches are parsed repeatedly and rewritten on each miss. Growing pair caches create increasing I/O cost and unsafe concurrent updates.
5. CSV inputs/results are retained in memory; output is rewritten after every row. Total written rows grow quadratically with input size. Environment deep copies add bounded per-episode overhead but do not solve dataset-scale accumulation.
6. No durable query/entity/paper/relationship indexes, batch inference, centralized queue, checkpoint resume or shared request deduplication exist.
7. One model instance is loaded per process at import time. Naively launching one process per pair wastes GPU memory and startup time; CPU fallback is not a practical throughput strategy for an 8B research loop.

Let N be input rows, U unique entities, Q reusable query plans, D unique documents and E candidate evidence units. Aim for O(N) streaming/indexed classification work, taxonomy near O(U), network work near O(Q + D), and inference near O(E), with bounded exceptions for unresolved cases. E can still grow with genuinely distinct relationships in a document; “once per paper” is not an unlimited-context promise. Candidate matching must avoid materializing all host-by-virus combinations mentioned in a large review.

## 17. API and rate-limit risks

The code uses four retry attempts, exponential waits, request timeouts, a one-second pause after NCBI source calls, and 0.4-second taxonomy delays. These are per-process controls. Search/EFetch can occur back-to-back within a source call, and multiple jobs multiply traffic. There is no account/IP-wide governor, request ledger, jitter, `Retry-After` handling, or durable retry schedule.

Retrying a whole source operation after fetch failure repeats the successful search. Broad exception handling loses typed distinctions among transport, parsing, missing data, and biological absence. Some taxonomy exhaustion messages can report an uninformative last error after HTTP-status retry branches. Configured email/tool values are fixed in source; API configuration is not centralized.

V2 should centralize source adapters, use persistent request keys and bounded retries with jitter, honor provider responses, coordinate quotas across workers, and record HTTP/parser/fetch outcomes separately. Cache successful empty searches with expiry; do not cache failed queries as empty success. Stage dataset downloads once and record digests. Provider-specific numeric limits and cluster network policies must be verified during implementation; no live policy verification was needed for this workspace-only audit.

## 18. Recommended V2 architecture

```text
stream CSV/Parquet -> input records and canonical pair mapping
  -> unique entity normalization and taxonomy resolution
  -> inspectable plausibility assessment
  -> join existing evidence and reusable coverage
  -> schedule only missing entity/corpus retrieval work
  -> canonical papers + structured text + entity mentions
  -> indexed candidate matching
  -> assertion extraction + endpoint/quote/role validation
  -> persistent evidence knowledge base
  -> deterministic classification + provenance + output partitions
```

The user-proposed order is sound. Add an early evidence/coverage lookup so cached knowledge can answer pairs immediately. Plausibility may prioritize work but should not suppress retrieval of potentially contradictory credible evidence without an explicit reviewed policy.

Suggested boundaries are input/contracts, entity resolver, plausibility rules, retrieval adapters/planner, document parser, entity matcher, evidence extractor/validator, evidence repository, coverage evaluator, classifier, and job coordinator. Pure biological decision functions should not import Torch or make network calls.

Retrieve primarily by unique virus and validated aliases, with reusable host/cohort searches where helpful and targeted pair queries only for unresolved coverage gaps. A virus-centered corpus can serve millions of pairs, but its coverage must be recorded explicitly; a top-five virus search does not confer adequate coverage on every host. Index entity mentions and accepted assertions to match only the actual input pairs.

Use inexpensive deterministic filtering and existing accepted evidence first. Restrict LLM work to ambiguous evidence units; batch by token budget in long-lived workers. Cache results by document content, section/unit identity, model/prompt/extractor version, entity-linker version and relevant taxonomy version. Classifier changes should rejoin stored evidence without redoing extraction; parser or role-model changes may invalidate specific units.

## 19. Recommended evidence model

The proposed fields are a good export view, but a single flat row is insufficient for multiple spans, unresolved names, comparison relations and provenance. Store assertions independently of target predictions, with normalized related tables and a convenient flattened export.

| Field group | Recommended fields |
|---|---|
| Assertion identity | `assertion_id`, `subject_entity_id`, `object_entity_id`, `relationship_type`, assertion revision/status |
| Host-virus export | `host_taxid`, `host_name`, `virus_taxid`, `virus_name`, resolved rank, resolution status; nullable TaxIDs with stable internal IDs |
| Entity grounding | Mention IDs, surface forms, linked entity candidates, resolution method/confidence, document-local abbreviation links |
| Biological meaning | `evidence_type`, `evidence_strength`, `natural_vs_experimental` = natural/experimental/in_vitro/unknown, sample/tissue, source-organism role, host role, assay, replication evidence |
| Assertion semantics | Polarity, uncertainty/speculation, primary finding/background/cited claim, negation, comparison-only flag, claim scope |
| Provenance | `paper_id`, paper-version ID, PMID/PMCID/DOI through identifier table, section/table/cell IDs, exact supporting text and character offsets, content digest, source/provider, retrieval time |
| Cross-passage links | Ordered supporting spans and an explicit sample/experimental-group/coreference link; do not concatenate unrelated snippets into a synthetic quote |
| Extraction | `extraction_method`, `extraction_version`, model/revision, prompt digest, parser/linker versions, timestamp, raw-output reference, attempt status and failure reason |
| Confidence/review | Separate entity, relationship and evidence-quality confidence, calibration version, human-review status/reason; optional combined confidence only under a documented policy |

Examples of relationship types: `DETECTED_IN_SAMPLE`, `ASSOCIATED_WITH_HOST`, `NATURAL_INFECTION_OF`, `EXPERIMENTAL_INFECTION_OF`, `REPLICATES_IN`, and distinct virus-virus `SEQUENCE_SIMILAR_TO`/`PHYLOGENETICALLY_RELATED_TO`. A sample-to-source-organism relation is separate from virus-to-host. Detection from a sample need not establish infection or replication.

One paper can have many assertions and one assertion many supporting spans. Several papers can repeat one cited experiment; retain citation lineage to avoid counting them as independent biological confirmations. Contradictory or negative assertions should be stored, not overwritten. Failed extraction is a processing outcome, never an assertion of absence.

## 20. Recommended storage schema

Use a normalized transactional store for identities/work state and immutable content artifacts for large documents. Begin with SQLite for a single-node prototype and controlled writer: it supports indexed joins, constraints and transactions without a separate service. Do not share a WAL-mode SQLite file among distributed HPC workers over a network filesystem. Use node-local shards plus deterministic merge, or a managed PostgreSQL service when concurrent multi-node writes are required. Large read-heavy analytical outputs should be partitioned Parquet; they do not replace transactional work state.

| Table/artifact | Key and purpose |
|---|---|
| `taxonomy_snapshot` | Snapshot ID, provider, date, digest, schema version. |
| `taxon` / `taxon_lineage` | Snapshot + TaxID unique, canonical name/rank/parent; indexed ancestry and merged/deleted-ID mapping. Hosts and viruses share taxonomy infrastructure. |
| `entity` / `entity_resolution` | Stable internal entity ID; raw query or supplied TaxID; candidate mapping, status/method/version/error/retry time. |
| `alias` | Entity, normalized alias, source/type/language/scope/ambiguity/version; many-to-many alias index, not globally unique alias strings. |
| `query_plan` / `retrieval_attempt` | Stable source/query/filters/version key; request status, timestamp, attempts, total hits, cursor, limits, error, expiry, completeness. |
| `query_paper` | Query/attempt -> canonical paper, retrieval rank and source provenance; unique mapping within attempt. |
| `paper` / `paper_identifier` | Canonical ID; unique identifier namespace/value with collision review and merge redirects. |
| `paper_version` / content files | Source/parser/content digest, access status, raw XML and normalized text locations; immutable versions. |
| `document_unit` / `entity_mention` | Section/paragraph/table/cell, offsets, surface forms, entity candidates and linker version. Index `(entity_id, paper_id)`. |
| `assertion` / `assertion_span` | Typed endpoints/roles/polarity/setting; exact supporting spans and validation outcomes. Index `(host_entity_id, virus_entity_id, acceptance_status)` for host-virus assertions. |
| `extraction_run` | Unique document-unit/content/extractor/model/linker version key; status, token/runtime metrics, raw response reference. |
| `plausibility_assessment` | Pair or reusable taxonomic rule result, rule version, supporting sources, resolution dependencies and reason. |
| `coverage_assessment` | Entity/corpus/search plan/date and pair applicability, required queries/units, unresolved gaps, adequacy rule/version. |
| `input_dataset` / `input_record` / `pair` | Input checksum/schema; immutable row ordinal; raw values; canonical pair mapping. Unique pairs processed once, original rows retained. |
| `classification` / `classification_evidence` | Pair + classifier/dependency signature; label/reason/confidence, coverage/plausibility IDs and accepted/counterevidence references. |
| `work_item` / `stage_attempt` | Idempotency key, stage, dependencies, pending/running/succeeded/retryable_failed/permanent_failed, lease owner/expiry, attempts, next retry, error. |
| `output_partition` / `run_manifest` | Input/config/model/taxonomy/source digests, row counts, checksum, commit marker, run IDs. |
| `cache_event` or metric aggregate | Hit/miss/stale/error by cache layer and version, counts/latencies; distinguish negative success from lookup failure. |

Apply foreign keys, uniqueness constraints and transactions; store timestamps and dependency versions. Keep large bodies and model responses once by digest rather than duplicating them in every classification row. A content digest alone does not identify the source license or retrieval provenance; retain both.

## 21. Batch-processing strategy

Stream CSV chunks and Parquet row groups. Preserve dataset ID and row ordinal; deduplicate canonical pairs in the database rather than an unbounded Python set. Resolve unique entities in batches. Join existing classifications/evidence/coverage before creating missing work.

Partition retrieval by unique entity/search plan and extraction by document unit so repeated pairs cannot independently enqueue the same external work. Use unique task keys with atomic claims. Process documents with many candidate entities through indexed mention matching against the supplied pairs, avoiding all-pairs expansion.

Use bounded CPU queues for parsing/linking, token-budgeted GPU batches for extraction, and backpressure between stages. Write immutable output partitions and reconcile them to original row order on export. Keep CSV support for small jobs and Parquet for large jobs. Pair duplicates must produce separate output rows referencing the same decision/evidence version.

Measure work reuse directly: unique entities, queries, papers and extraction units; external requests; LLM tokens; hit/miss/stale counts; throughput; peak memory; and wait/retry time. A million rows containing one repeated pair should not cause a million research episodes.

## 22. Checkpoint/resume strategy

Use durable stage state rather than “last CSV row.” A run manifest pins input checksum and all meaningful versions. Each stage's idempotency key includes its dependency versions. Commit stage results and success state atomically. Never mark success before durable evidence/content is committed.

Workers claim leased items; expired leases become eligible for retry. Handle retryable failures with a persisted schedule and terminal failures with explicit output status. At-least-once execution with unique result constraints is sufficient; avoid promising impossible exactly-once network calls after crashes.

Write content/output partitions to temporary files, verify checksums/row counts, atomically finalize, then register completion. On resume, reconcile finalized artifacts with manifests, reuse completed stages, and recompute only invalidated dependencies. A classifier update should reuse paper/extraction results; a new taxonomy snapshot should re-resolve impacted entities and dependent decisions.

On scheduler termination, stop claiming work, finish or release a bounded active batch, and persist progress. Tests must kill workers during fetch, extraction, database commit and output finalization, then verify neither skipped rows nor duplicate assertions appear.

## 23. HPC execution strategy

Separate internet-facing retrieval from compute-heavy extraction. Stage versioned taxonomy, Virus-Host DB, papers and model assets in an approved location accessible to compute jobs; use retrieval nodes/services permitted by the site's network policy. Do not launch unrestricted external requests from every array task.

Use SLURM arrays over immutable document/input partitions with explicit resource profiles. Load one model per long-lived GPU worker and batch inference; CPU workers handle normalization and indexing independently. Parameterize model paths, scratch/output directories, account and GPU choices. Resolve workspace-relative configuration rather than embedding one user's cluster path.

For a service-free deployment, give each worker a node-local SQLite shard and write immutable result partitions, then merge through a controlled coordinator with unique keys. For distributed online coordination, use an approved PostgreSQL endpoint and a shared request governor. Avoid concurrent JSON writes and shared network-filesystem SQLite WAL.

Package a portable minimal dependency specification separately from Alliance-specific build locks. Pin model/tokenizer revisions, parser/extractor versions and data checksums. A syntax check does not replace integration validation. Record resource use, warm/cold-cache state and actual stage costs; the current 60-analysis benchmark is too small to project production capacity reliably.

## 24. Testing strategy

Build deterministic offline tests before changing semantics. Biological decision functions should be importable without model initialization or network access. Freeze representative source text, entity mappings and extracted assertions, with provenance and reviewed expected evidence labels. Keep technical fixtures clearly separate from real literature evidence.

Test layers:

1. Unit tests for normalization, alias scope, TaxID conflicts, identifier merging, relationship endpoint/type/setting validation, negation, coverage gates and precedence.
2. Regression replays of saved problematic evidence: v05 comparison types on host edges; v06 comparison-only `SEQUENCED_FROM`; short exact quotes; reversed comparison endpoints; failed extraction -> insufficient; escaped grounding regex.
3. Adapter integration fixtures for successful empty results, malformed XML/JSON, timeout, 429/5xx, partial fetch, pagination, duplicate identifiers and missing full text.
4. Store/coordinator tests for transactional upserts, simultaneous duplicate tasks, cache invalidation, crash recovery, leases, interrupted exports and stable output multiplicity.
5. Small pinned-model extraction evaluation with human-reviewed spans/roles, separate from deterministic classifier evaluation.
6. Performance tests with repeated entities/papers and increasing distinct evidence volume. Assert bounded memory and reuse, not an arbitrary speedup over uncontrolled live runs.

Keep the existing twelve labels intact for historical comparison. Add a separately versioned V2 gold set with evidence-level annotations and four-class labels under the approved policy. Use paper/virus-family/source-disjoint evaluation splits where practical to reduce leakage. Report all-input and completed-case metrics, per-class precision/recall, false-known/novel/implausible counts, abstention/failure rates, coverage, taxonomy resolution, and confidence calibration only when sample size supports it.

## 25. Benchmark expansion

| Required case | Fixture and expected assertion |
|---|---|
| Exact known | Reviewed primary evidence naming the exact endpoints; `KNOWN` with the correct supporting span and relationship scope. |
| True novel candidate | Prospectively reviewed plausible pair with a dated, adequate search record and no credible exact support. Label is a candidate within that scope, not a proof of universal novelty. |
| Unrelated host | Paper associates target virus with another host; no target support. Final label depends on independent plausibility and coverage. |
| Related host | Same-genus or same-family host has evidence; no downward species transfer. |
| Related virus | Target host has another virus; relatedness does not establish exact support. |
| Comparison virus | Similarity passage and separate discovery passage; forbid comparison target on host edge. |
| Environmental detection | Virus in water/soil/pool/source material; retain sample detection, do not infer host. |
| Experimental infection | Exact challenged host and outcome; preserve experimental setting and scope. Include failed challenge/control. |
| Natural infection | Exact naturally infected host with appropriate evidence; distinguish from laboratory/cell culture. |
| Sequence similarity only | Store comparison relation, never infection/host association. |
| Alias/synonym | Verified renamed entity plus ambiguous acronym and shortened genus controls; accept only justified equivalence. |
| Unresolved taxonomy | Ambiguous/no-match/conflicting TaxID; insufficient strong-decision basis, never automatic novelty. |
| Zero search results | Successful empty requests with incomplete coverage -> insufficient; separately test adequate coverage plus plausible biology. |
| API failure | All-source and partial-source failures; no false empty-success or novelty. |
| Ambiguous evidence | Multiple hosts/viruses, pronouns, citations, contradictory claims; uncertainty retained. |
| Duplicate input pair | One internal computation, one output per original row, shared evidence/version. |
| Repeated host | Entity resolution and aliases reused across viruses. |
| Repeated virus | Retrieval/coverage reused across hosts, without claiming more coverage than established. |
| Repeated paper | One canonical paper/content unit and versioned extraction reused across searches/pairs/runs. |
| Genus/family evidence | Higher-rank assertion retained without species-level `KNOWN`. |
| Background/earlier study | Cited claim distinguished from current experiment; original-study linkage and no double counting. |
| Technical extraction failure | Empty/malformed/role-invalid response and exhausted rescue -> material processing failure, not biological non-support. |
| Pagination/candidate cap | Decisive paper beyond initial limits; incomplete budget must prevent novelty. |
| Grounding boundaries | Hyphens/spaces, inline XML, tables, abbreviations, negation and separate experimental arms. |

For categories that do not determine a single final label, annotate the expected evidence outcome and prerequisites separately. Do not force every unrelated-host or comparison case into novelty or implausibility. Preserve benchmark labels; corrections to future gold annotations require independent evidence review and explicit versioning, never an accuracy-driven edit.

## 26. Migration plan

| Disposition | Components/behavior |
|---|---|
| Preserve | Current source/query provenance; PMID/PMCID/DOI; all-query pooling before selection; taxonomy synonyms/lineage; environmental/source-organism separation; prior-versus-evidence distinction; separate host/comparison edges; extraction attempt diagnostics; deterministic final decisions; original benchmarks and results. |
| Refactor | Pure normalization/resolution interfaces; shared network adapters; parser and passage indexing; alias policies; selective rescue; versioned evidence validation; benchmark runners and output contracts; model loading into an explicit service/worker. |
| Replace | Pair-at-a-time research orchestration, count-only sufficiency, substring identity, whole-paper keyword host verification, automatic quote repair as semantic validation, fixed confidence rules, JSON caches, repeated CSV rewrites, model-generated summary on every pair. |
| Retire but retain temporarily | `pubmed_search.py` duplicate path; legacy four-label judge and environment/reward harness as production entry points; `test_model.py` as a purported test; mislabeled historical launchers and snapshots. Keep them read-only for comparison until replacements are validated. |

After approval, first pin the reproducibility record and fixtures. Then introduce typed data/failure contracts alongside existing code, without wholesale replacement. Move pure decisions behind testable interfaces; fix endpoint/coverage correctness; introduce persistent entities/papers/assertions; add batch/resume; finally migrate CLI/HPC workflows. Shadow-run old and new decisions against frozen inputs and explain every changed outcome. Keep a separate legacy-label comparison rather than distorting V2 semantics to match old related/absence labels.

Existing caches may be imported as historical records with unknown snapshot dates/status, not as fully trusted current truth. Do not backfill missing provenance or mark cached unresolved results successful. Re-extract affected documents when the old evidence cannot substantiate a structured assertion.

## 27. Implementation milestones

| Milestone | Deliverable after approval | Acceptance gate |
|---|---|---|
| **M1: Reproducible evidence baseline and contracts** | Immutable source/data/model manifests where available; honest v0.4 recorded-reference designation; offline evidence fixtures; typed resolution/retrieval/extraction/classification contracts and proposed association-scope policy. | Existing labels untouched; fixtures expose comparison promotion, failure-to-nonsupport and direct-positive losses; decision tests require no model/network. Missing historical assets are explicitly listed. |
| M2: Correctness boundary | Endpoint-grounded assertion validation, distinct processing failures, coverage gates and four-class policy behind an isolated interface. | No false exact support on reviewed adversarial fixtures; known direct/cross-sentence evidence retained; failures/unresolved taxonomy cannot produce novelty. |
| M3: Reusable knowledge store | SQLite schema, entity/alias records, canonical papers/identifiers, immutable document units, assertions and versioned provenance. | Repeated paper/entity work reuses stored results; collisions and invalidation tested. |
| M4: Retrieval and extraction jobs | Reusable entity queries, pagination/coverage ledger, coordinated retries, bounded model batches. | Requests and extraction count scale with distinct work; truncated/inaccessible evidence remains visible. |
| M5: Production batch/resume CLI | CSV/Parquet ingestion/export, input multiplicity, durable work state and restartable partitions. | Crash/restart yields complete outputs without duplicate assertions or lost input rows; bounded memory. |
| M6: HPC pilot and expanded benchmark | Site-configurable jobs, pinned assets, controlled network workers and representative scaling runs. | Reviewed accuracy/error gates plus measured throughput, resource use and reuse; no rollout to millions before these pass. |

M1 is the first implementation milestone. It prevents another undocumented “better version” from hiding worse entity attribution. The audit itself does not implement any milestone.

## 28. Risks and open questions

- **Association scope:** does `KNOWN` include credible discovery/detection and experimental infection, or must it mean natural biological host? The specification allows association/discovery evidence but requires distinctions. Recommended policy is to preserve all scopes and make the accepted scope explicit; never label experimental-only evidence natural. Approve this before final gold labels.
- **Unregistered viruses:** how should a publication-defined virus without a public TaxID be adjudicated? Preserve a provisional stable identity and evidence, while withholding unjustified taxonomy certainty and novelty.
- **Historical reproducibility:** authenticated v0.4 source/model/query responses are missing or not linked. Available backup names are unreliable. Saved result changes are observed; not all causal attributions can be proved.
- **Coverage:** define the minimum corpus/search plan and refresh interval by biological domain. No finite search proves absence everywhere. Abstract-only, inaccessible supplements and unresolved citations need explicit limitations.
- **Plausibility:** broad host-group agreement/disagreement is insufficient. Strong rules require biological review, citations and versioning; unknown mechanisms must not be treated as impossibility.
- **Curated databases and leakage:** distinguish production knowledge import from held-out literature evaluation; the current hidden target-record mechanism is only one benchmark control.
- **Evidence dependence/conflicts:** repeated papers and citations may describe one experiment. Retractions, corrections and contradictory studies require versioned reconsideration.
- **Storage operations:** choose node-local shards versus a managed database based on actual HPC permissions and concurrency, not projected row count alone.
- **Model behavior:** source material can contain misleading instructions or text artifacts. Treat it as data; accept only schema-valid, source-grounded assertions. Record truncation, model/prompt versions and unresolved role warnings.
- **Data access:** persist permitted content and provenance, and record inaccessible documents as coverage gaps. Raw source snapshots/model assets are absent here and will be needed for reproducible live evaluation.

Approval of this review should precede V2 implementation. No source code, benchmark labels, project instructions, historical files or caches were modified during this audit.
