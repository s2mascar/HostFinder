# Twelve-pair correctness baseline and diagnostics

## Result and limits

**The latest available benchmark artifact scores 6/12 (50%). This is the saved v0.6 result, not a fresh run of the current source.** There are two false-negative known cases, three missed related cases, and one related case incorrectly promoted to `KNOWN`. All four negative labels match. Label agreement does not establish that the underlying evidence attribution is correct.

The project instructions and architecture review were read before this phase. No classifier, extractor, retrieval rule, expected label, or biological rule was changed. The historical v0.4 score remains 8/12; its purported backup is not an authenticated v0.4 executable snapshot. The objective remains 12/12 through correct attribution, not through pair-specific exceptions.

Environment preflight found no `python`, `python3`, `py`, or `bash` command, and no workspace `models/Qwen3-8B`, `agent_env`, `data`, or `logs` directory. The required model cannot load here. There are no raw-paper or model-response snapshots that permit end-to-end offline replay. No external requests, model inference, or live benchmark were attempted after this decisive preflight. Existing JSON taxonomy/context caches cannot substitute for missing literature and model assets.

Instead, the diagnostic phase **executed a historical artifact export and validated it**. This reproduces the saved score and preserves the evidence records, not the historical model computation. It would be misleading to report the current executable as freshly scoring 6/12, or to claim that all historical errors have been reproduced in Python.

## New tooling and reproduction

Files added:

- `export_benchmark_diagnostics.ps1`: parameterized, model-free exporter. Validates input/result alignment, unchanged expected labels, stored correctness flags, unique per-episode paper IDs, paper counts, and exact/related counts. Records SHA-256 dependencies and refuses output overwrites. It does not call or reimplement classification.
- `results/benchmark_diagnostic_v06.json`: all 12 cases, all 60 original per-paper diagnostics, full original result rows, historical taxonomy/context snapshots, explicit unknown fields, and dependency manifest.
- `test_benchmark_diagnostics.ps1`: verifies all rows and all original paper JSON, quote/edge fidelity, expected/predicted labels, missing-field honesty, hashes, and overwrite refusal.
- This report.

Run from the workspace root in PowerShell:

```powershell
./test_benchmark_diagnostics.ps1
# Choose a NEW output path for a fresh export; existing snapshots are protected.
./export_benchmark_diagnostics.ps1 -OutputJson results/benchmark_diagnostic_v06_reexport.json
```

The second command is an example for later reproduction; the delivered export is `results/benchmark_diagnostic_v06.json`. Inputs default to the unchanged `test_pairs.csv`, saved v06 result, and v06 context/virus caches. No benchmark names or per-pair exceptions appear in the exporter. Alignment is deliberately strict by input ordinal: a shuffled or relabeled result file is rejected rather than silently paired to the wrong row. Future exports of other schemas may need an explicit adapter; this utility requires complete per-paper diagnostics.

Executed validation result: **PASS: 12 rows, 60 unaltered paper records, 6 correct; hashes, unknown fields and overwrite protection verified.** Four saved failed extractions are retained as failed, even though their historical classifications are `NO_SUPPORT`.

Frozen input SHA-256: `be0f564702c1d8ba2654e7f2aa920cac8532060df79a61dac7a89ae8a550106b`.

Frozen v06 result SHA-256: `a8b0cee60f260d5cc676a21f732ff80baeb11facfe27203458acad203f8b8922`.

The manifest records current source hashes for future comparison. **It does not assert that those source hashes or the supplied caches generated the historical CSV.** Original run date, model revision, tokenizer revision, raw responses, and runtime environment were not saved. Model behavior cannot be made retroactively reproducible by adding a manifest now.

## Definitions, evaluators and output availability

The CSV defines four positive exact cases, four related-virus cases, and four negative cases. `expected_evidence` is descriptive and is not checked by either evaluator. Its `RELATED_ONLY` corresponds conceptually to the extractor's `TARGET_HOST_RELATED`, not exact association. `reference_hint` is not asserted as a retrieved/supporting paper. Purposes below are inferred from the provided label, test type and reference hint, not newly adjudicated biological gold labels.

`evaluate_pairs.py` invokes `run_judge_agent` and saves host, virus, expected/predicted status, correctness, paper counts, confidence, reason, runtime and benchmark descriptors. It omits individual extraction diagnostics. `evaluate_env.py` invokes the baseline environment policy and additionally saves reward/steps/search counts, biological-prior summaries, taxonomy family/genus fields, environmental context, extraction-status counts, edge counts, and per-paper JSON. Both write the full accumulated CSV after each episode. Neither pins models, caches or retrieved texts, nor resumes a run from a manifest. Environment headline accuracy excludes `ERROR` rows; this report scores every input row, and none of these 12 saved rows is `ERROR`.

`run_benchmark.sh` executes current `evaluate_env.py` into a v06-named CSV and optionally the current legacy path. Its SLURM job name still says v05. `run_benchmark_v2.sh` executes current code into v04-named files; it is not a version selector or a V2 pipeline. Both embed an HPC path and assume a project environment and GPU. Do not use their output suffixes as proof of source identity. No saved logs were found in this workspace. All six historical result CSVs remain intact.

| Requested diagnostic | Historical availability / export treatment |
|---|---|
| Host, virus, expected/predicted, correctness | Complete, validated against the definition CSV. |
| Taxonomy resolution | Per-paper target-virus TaxID/name/aliases plus separately labeled historical context/cache snapshots. Original resolver attempts and cache identity unknown. |
| Host aliases used | Not recorded in v06 CSV; host alias cache absent. Null, not reconstructed from names in quotes. |
| Virus aliases used | Preserved separately for every paper; no invented union that obscures which paper used which aliases. |
| Literature queries attempted | Search count exists, query strings do not. Null strings; never infer exact/model queries from planner code. |
| Search success/failure | HTTP/source statuses and `source_successes` absent from saved CSV. Unknown, even when a final label implies the old sufficiency check passed. |
| Number of papers retrieved | Saved count is selected candidates, not raw hits. Raw pool and preselection count unknown. |
| Candidate papers | All selected/analyzed paper records preserved; discarded candidates unavailable. |
| Exact/related support | Paper IDs derived solely from saved classes and validated against saved counts. |
| Edges and bindings | Host/comparison objects, source/target names, passage verification and match flags retained in full. These are model/verifier claims, not independent validation. |
| Evidence type/scope | Claimed relationship types and support mode retained. Natural/experimental/sample scope was not separately recorded and remains null. Do not infer natural infection from `STUDY_CONTEXT`. |
| Confidence/reason/runtime | Original values retained, no recalibration. |
| Cache usage | Cache snapshots available; actual hits/misses unknown. A file's presence does not prove a cache hit. |

The caches resolve all six host names in the supplied context records. The Lampyris partiti-like virus 1 is unresolved in the virus cache; other benchmark virus records have TaxIDs, with some ranked `no rank`. Mouse hepatitis virus has Murine hepatitis virus aliases. These are useful context, but unresolved/public-taxonomy identity must remain distinct from a paper-defined virus name. There is no demonstrated cache-alias root cause for all six wrong labels.

## Scorecard

`NEF` below means the legacy label `NO_EVIDENCE_FOUND`, never novelty. `PK` means `POSSIBLY_KNOWN`. Candidate counts are selected/analyzed papers. Runtime is historical episode time, not the exporter runtime.

| # | Host | Virus | Expected | Actual | Correct | Papers | Confidence | Seconds |
|---:|---|---|---|---|---|---:|---|---:|
| 1 | Lampyris noctiluca | Lampyris noctiluca partiti-like virus 1 | KNOWN | KNOWN | Yes | 1 | HIGH | 42.14 |
| 2 | Arabidopsis thaliana | Turnip mosaic virus | KNOWN | NEF | No | 8 | MEDIUM | 70.12 |
| 3 | Mus musculus | Mouse hepatitis virus | KNOWN | NEF | No | 8 | MEDIUM | 72.26 |
| 4 | Homo sapiens | SARS-CoV-2 | KNOWN | KNOWN | Yes | 8 | HIGH | 85.47 |
| 5 | Lampyris noctiluca | Hubei partiti-like virus 31 | PK | KNOWN | No | 1 | HIGH | 31.45 |
| 6 | Lampyris noctiluca | Hubei partiti-like virus 51 | PK | NEF | No | 1 | MEDIUM | 34.27 |
| 7 | Lampyris noctiluca | Hubei chuvirus-like virus 3 | PK | NEF | No | 6 | MEDIUM | 66.76 |
| 8 | Lampyris noctiluca | Hubei toti-like virus 16 | PK | NEF | No | 2 | MEDIUM | 43.25 |
| 9 | Tribolium castaneum | Hubei partiti-like virus 31 | NEF | NEF | Yes | 1 | MEDIUM | 24.00 |
| 10 | Arabidopsis thaliana | Mouse hepatitis virus | NEF | NEF | Yes | 8 | MEDIUM | 70.64 |
| 11 | Mus musculus | Turnip mosaic virus | NEF | NEF | Yes | 8 | MEDIUM | 59.69 |
| 12 | Aedes aegypti | Mouse hepatitis virus | NEF | NEF | Yes | 8 | MEDIUM | 63.06 |

Total recorded episode runtime: 663.11 seconds. Historical comparison: initial environment 7/12; v02 6/12; v03 7/12; v04 8/12; v05 6/12; v06 6/12. The same aggregate score can hide very different failures: v05 had five false `KNOWN` predictions; v06 has one and loses two additional positive cases.

## Case-by-case diagnosis

Paper IDs `P0`, `P1`, etc. are local to each case. PMID/PMCID and title identify the source artifact record, not a new external lookup. All selected papers and full edges are available in the JSON; this report highlights decision-driving records. Quotes below describe stored passages; the full original articles were not independently fetched.

### 1. Lampyris noctiluca -> Lampyris noctiluca partiti-like virus 1

**Expected/actual: KNOWN / KNOWN; label correct.** Positive discovery case, reference hint PMC7093385. Tests an exact host-associated virus described across discovery-study passages, including a paper-defined name without resolved public taxonomy.

The sole paper, P0 / PMID 31900852, *Identification and characterisation of common glow-worm RNA viruses*, was RESCUED and classified `EXACT_SUPPORT`. Host edge: Lampyris -> Lampyris partiti-like virus 1, `IDENTIFIED_IN`, `STUDY_CONTEXT`, verified. The virus span describes LnoPLV1 genome length/protein and then a comparison to Hubei 31. The separate comparison edge has an incorrect source assignment and is unverified, but exact-support precedence does not require it. HIGH confidence comes from the deterministic judge.

This is a necessary guard against an overstrict single-sentence-only fix. However, the host span is a general introductory description; the current broad-context verification is weaker than an explicitly linked specimen/discovery assertion. Preserve the correct exact conclusion by recovering the real study/sample linkage, not by allowing every virus mentioned in this paper to inherit its host. Evidence scope: discovery/sequence association; the saved result does not establish natural infection or replication separately.

### 2. Arabidopsis thaliana -> Turnip mosaic virus

**Expected/actual: KNOWN / NO_EVIDENCE_FOUND; incorrect.** Positive exact interaction, reference hint PMID17427806. Tests recovery of ordinary host-virus evidence without confusing other plant hosts, viruses or experimental scope.

Seven papers become `NO_SUPPORT`, one `VIRUS_OTHER_HOST`; no verified host edge remains. The result is not explained by zero retrieved papers.

- P0, *Turnip mosaic virus drives selective filtering and community reassembly...*: target-host and target-virus matches are true, but the host and virus-effect quotes are separate and the virus quote does not contain the claimed infection wording. Support mode is NONE. The exact but short quote is not widened by the grounding helper.
- P2 / PMID42372481: the quoted Arabidopsis/TuMV pathosystem includes both names but lacks the literal `INFECTION_OF` pattern. It becomes `RELATED_LANGUAGE_BUT_INCOMPLETE_EDGE_CHAIN` because the model also labels a virus synonym as `SAME_GENUS` comparison.
- P1 / PMID42265533: the stored larger comparison passage explicitly discusses both Nicotiana and Arabidopsis, including TuMV propagation in Arabidopsis. The host-virus quote is only a short infection phrase. A `DIFFERENT_HOST` comparison label with Nicotiana yields `VIRUS_OTHER_HOST`, incorrectly reducing a multi-host record.
- P3 (experimental evolution) extracts only the virus name for its experimental-infection quote; P6 / PMID42671063 uses TuMV and A. thaliana in a longer study statement but labels it natural infection. P5 / PMID42347223 uses TuMV-GFP and Potyvirus rapae nomenclature. These expose local acronym/name and relationship-scope problems rather than proving that broad aliases should be accepted.
- P7 / PMID25018765 correctly distinguishes the Brassica/TuMV table entry from Arabidopsis/Turnip crinkle virus. It must remain a negative evidence fixture for the exact target pair.

**Observable root causes:** incomplete extracted assertion spans, overly literal type validation, no correction of otherwise exact-but-insufficient spans, and incorrect multi-host/comparison labeling. All-nonsupport aggregation then produces absence. The reference-hint paper is not among these selected records; historical query/rank logs are missing, so why it was omitted is unknown. Enough potentially relevant material is already present to locate failures after retrieval, but sufficiency of overall coverage cannot be certified.

**Categories:** extraction failure; entity binding error; other-host contamination; judge aggregation failure; other (type/quote mismatch). Alias/acronym handling is a contributor in specific records, not a demonstrated wrong NCBI mapping. **Responsible:** `evidence_agent.py` `analyze_paper`, `_ground_extraction_passages`, `classify_relationship` direct support and `explicit_comparison_other_host`; `bioresearch_env/relationship_language.py`; `judge_agent.py` aggregation. Search selection in `search_agent.py` is a coverage question, not a proven primary cause.

### 3. Mus musculus -> Mouse hepatitis virus

**Expected/actual: KNOWN / NO_EVIDENCE_FOUND; incorrect.** Exact positive and synonym/abbreviation case, reference hint PMID8677422. Tests mouse/Mus and Mouse/Murine hepatitis naming, ordinary infection evidence, and distinction from receptor or cell-line experiments.

All eight papers become `NO_SUPPORT`.

- P0 / PMID8677422 is the actual reference-hint paper. Both target matches and both quote-verification flags are true. The short quote identifies the virus as the mouse coronavirus, but the model labels it `NATURAL_INFECTION`; the type validator requires explicit natural-infection wording. Direct support fails; natural infection is not in the contextual fallback set. A synonym is mislabeled `SAME_GENUS`, resulting in an incomplete comparison-chain class. This demonstrates post-retrieval loss directly.
- P3 / PMID15953191 extracts only `MHV3` as the virus quote but reports Mouse hepatitis virus as the entity and natural infection as the type. Full-name passage matching and relationship wording fail; the title discusses experimental infection. P4 / PMID15638127 and P5 / PMID10803365 have unverified reconstructed virus quotes.
- P6 / PMID25502736 retains a virus prevalence fragment with a mouse count. P7 / PMID21508117 retains a sentence describing organisms found in pet-shop mice. Both are labeled `DETECTED_IN`, but the lexical validator does not accept these particular formulations. Do not fix this by equating every prevalence or detection signal to natural infection; retain assay scope.
- P1 / PMID1719235 has a JSON extraction failure, empty virus fields and `NO_SUPPORT`. Its receptor-expression experiment in human/hamster cells is not a substitute for mouse natural-host evidence.

**Observable root causes:** quote/type mismatch on the retrieved reference, acronym grounding failures, reconstructed quotes, and extraction failure collapsed into non-support. The cached Mouse/Murine aliases and target-match flags work in several failed records, so failed taxonomy lookup is not the explanation. **Categories:** extraction failure; entity binding error; judge aggregation failure; other (literal type validation). **Responsible:** `evidence_agent.py` span grounding, type validation and failed-extraction path; `bioresearch_env/relationship_language.py` `NATURAL_INFECTION`/`DETECTED_IN`; `judge_agent.py` failure-blind aggregation.

### 4. Homo sapiens -> SARS-CoV-2

**Expected/actual: KNOWN / KNOWN; label correct, evidence quality not certified.** Positive human normalization case, reference hint PMID32699094. Tests human clinical language and target-virus aliases while excluding unrelated organisms and viruses.

Only P0 / PMID42743972, *PhoSARte: identification of SARS-CoV-2 phosphorylation sites using contrastive learning and protein language models*, is `EXACT_SUPPORT`, RESCUED, `DETECTED_IN`, `STUDY_CONTEXT`. Its host passage combines cited work on human cell lines with Vero E6/C. sabaeus and other experimental systems. Its virus quote is background about kinase signaling. The comparison passage discusses adenovirus type 2 in human IMR-90 cells. The remaining seven candidates contain no accepted exact support (three mention-only, four no-support).

The known label agrees with the benchmark, but the driving record does not independently demonstrate natural human infection. The reference hint is not the supporting paper. A binding/scope repair may correctly reject this current support and temporarily lower label accuracy until sound evidence is recovered. P6 / PMID42645657, about a human blood-brain barrier model in the absence of infection, must not become positive merely because it contains human cells and the target virus. No label is changed here. Risk categories: co-mention false positive, other-host contamination, evidence-scope error (other), confidence/calibration failure. This is a warning about evidence acceptance, not a claim that the biological pair is unknown.

### 5. Lampyris noctiluca -> Hubei partiti-like virus 31

**Expected/actual: POSSIBLY_KNOWN / KNOWN; incorrect.** Related-virus case, reference PMC7093385. Tests the distinction between a virus discovered in the target host and the comparison virus from another source.

Sole P0 / PMID31900852 is RESCUED and classified exact. The host edge wrongly names Hubei 31 with `SEQUENCED_FROM`, `STUDY_CONTEXT`. Its virus quote states BLAST similarity to Hubei 31, not sequencing from Lampyris. The comparison edge names Lampyris partiti-like virus 1 as source and Hubei 31 as target, with spider mix as the comparison source context, but is unverified because the host-edge virus and comparison source differ.

**Exact decision cause visible in the artifact:** the verifier accepts global study context plus a permitted host-type label and a correctly copied comparison passage. `host_edge_verified && host_virus_matches_target` triggers exact support before the contradictory comparison edge matters. The judge then gives `KNOWN/HIGH`. The correct chain is Lampyris -> its study virus -> comparison Hubei 31; this must not become Lampyris -> Hubei 31.

**Categories:** entity binding error; co-mention false positive; comparison-virus contamination; sequence-similarity error; confidence/calibration failure. **Responsible:** `evidence_agent.py:850` contextual support, rescue validation and exact-support branch; `judge_agent.py` exact precedence and fixed HIGH confidence. This is not fixed by banning only `SEQUENCE_SIMILARITY` as a host type: the recorded wrong type is already `SEQUENCED_FROM`.

### 6. Lampyris noctiluca -> Hubei partiti-like virus 51

**Expected/actual: POSSIBLY_KNOWN / NO_EVIDENCE_FOUND; incorrect.** Related-virus case, reference PMC7093385; tests the direction of the two-edge chain.

Sole P0 / PMID31900852 is RESCUED. Its host edge is verified: Lampyris -> Lampyris partiti-like virus 2, `IDENTIFIED_IN`, `STUDY_CONTEXT`. Its comparison passage contains similarity/phylogenetic language and the target Hubei 51, with Chinese land snail mix as the target's source context. However, the model sets **comparison source = Hubei 51**, rather than the host-associated LnoPLV2. `source_matches_host_virus` is false; the comparison edge is unverified; the paper becomes `RELATED_LANGUAGE_BUT_INCOMPLETE_EDGE_CHAIN` and the judge returns absence.

**Observable root cause:** reversed/self-referential comparison source, with no rescue triggered merely by this endpoint conflict. Verification correctly refuses the malformed edge; loosening the comparison-source equality check would be an unsafe fix. **Categories:** entity binding error; extraction failure. **Responsible:** `evidence_agent.py` extraction/rescue and comparison-source linkage (`:912`), with downstream expected aggregation in `judge_agent.py`. The expected repair requires an explicitly grounded local abbreviation/coreference chain, not an unconditional swap based on the requested target.

### 7. Lampyris noctiluca -> Hubei chuvirus-like virus 3

**Expected/actual: POSSIBLY_KNOWN / NO_EVIDENCE_FOUND; incorrect.** Related-virus case, reference PMC7093385; tests target/reference-virus role separation and distinguishing Odonata source context from the target host.

- P0 / PMID31900852, the discovery paper, names Hubei chuvirus-like virus 3 on the host edge and Lampyris chuvirus-like virus 1 as the comparison target: the biological roles are reversed. Its host quote describes two Finnish populations and tissues without naming the host, so `target_host_match=false`; the comparison target also fails to match the requested virus. It becomes `NO_SUPPORT`.
- P1 / PMID40512168 and P3 / PMID37622664 are taxonomic updates, classified no-support. A host name occurs inside a virus name in a taxonomy list, which is not host-association evidence. Their quoted virus text concerns Hubei chuvirus-like virus 1, despite the requested virus 3. These should not be rescued by generic genus/name matching.
- P2 / PMID34463877 fails JSON extraction; P4 / PMID36437428 returns no usable rescue entities. Both become `NO_SUPPORT`. P5 / PMID35108077 is mention-only with unverified target passages.

**Observable root causes:** P0's reversed roles and unresolved cross-passage host context prevent the related chain; two independent extraction failures are treated as biological non-support. The paper is retrieved, so its absence is not the cause. The original source/model inputs are unavailable, preventing a claim about whether truncation, prompt bias or model generation originally caused the role reversal.

**Categories:** entity binding error; comparison-virus contamination; other-host contamination; extraction failure; judge aggregation failure. **Responsible:** `evidence_agent.py` host grounding, extraction/rescue and comparison linkage; `judge_agent.py` failure handling. The separate source-organism cache correctly treats Odonata context as distinct from a resolved biological host; no cache environmental-source bug is established here.

### 8. Lampyris noctiluca -> Hubei toti-like virus 16

**Expected/actual: POSSIBLY_KNOWN / NO_EVIDENCE_FOUND; incorrect.** Related-virus case, reference PMC7093385; tests that an other-host statement about a comparison virus can coexist with a valid related chain from the target host.

P0 / PMID31900852 becomes `VIRUS_OTHER_HOST`. The host quote is the same unnamed Finnish-population fragment as case 7. The reported host-edge virus is incorrectly Hubei toti-like virus 16 with `SEQUENCED_FROM`. Its passage describes sequence/phylogenetic similarity and isolation from spiders. The comparison source is Lampyris totivirus-like virus 1, target Hubei 16; because the host-edge virus is wrong, `source_matches_host_virus=false`. With no verified target-host context, the direct-other-host branch fires. The existence of comparison-virus evidence from spiders is not itself wrong; losing the separate Lampyris/LnoTLV1 chain is wrong.

P1 / PMID35337056 is a canegrub transcriptomics study. Its RESCUED extraction retains `SEQUENCE_SIMILARITY` on the host edge and an unverified virus quote; it is `NO_SUPPORT`. This non-target-host paper must not be promoted to target-host evidence to repair P0.

**Observable root causes:** wrong host-associated-virus role, missing host binding across study passages, comparison-source mismatch, and consequent loss of a related chain. **Categories:** entity binding error; comparison-virus contamination; other-host contamination; sequence-similarity error; extraction failure. **Responsible:** `evidence_agent.py` extraction/rescue, host/context verification and non-target-host branch; downstream `judge_agent.py` absence aggregation. Do not remove the legitimate other-host distinction; represent both relationships correctly.

### 9. Tribolium castaneum -> Hubei partiti-like virus 31

**Expected/actual: NO_EVIDENCE_FOUND / NO_EVIDENCE_FOUND; label correct.** HostFinder negative candidate, tests target-host substitution and comparison contamination when a highly relevant virus paper is about another host.

Sole P0 / PMID31900852 is `NO_SUPPORT`. The model claims Tribolium, but the host quote describes another virus/Odonata context and is unverified; target-host binding is false. The Hubei 31 similarity quote is verified, but this cannot establish the target edge. The comparison target is a Lampyris virus, not the requested target. This is a critical regression guard: v05 falsely called the case known by trusting the model host name and study context. MEDIUM confidence is the old fixed absence confidence, not proof of adequate search. No novel-candidate conclusion is warranted.

### 10. Arabidopsis thaliana -> Mouse hepatitis virus

**Expected/actual: NO_EVIDENCE_FOUND / NO_EVIDENCE_FOUND; label correct with an unresolved technical failure.** Cross-group decoy, tests plant/mammalian-virus co-mentions and biological-prior leakage.

Seven no-support and one mention-only records. P0 / PMID39605984, a fibrillarin review, has a JSON extraction failure and empty edges. Other papers discuss plant pathways and mouse-virus biology in distinct contexts: P3 / PMID38599165 discusses NLR sensors; P4 / PMID40828816 concerns tomato bushy stunt virus with coronavirus comparisons; P6 / PMID32553580 is mention-only. P7 / PMID32699064 includes model-generated “quotes” taken from the biological prior rather than the paper; verification rejects them.

The target edge is not established, which matches the legacy label. However, one failed extraction cannot certify absence. A proposed failure-aware judge could temporarily return `UNCLEAR` until that relevant candidate is successfully processed or independently excluded. Do not add an exception for this decoy to preserve its score. Keep prior text from becoming literature evidence. The cached discordant prior is not the judge's proof of absence or implausibility.

### 11. Mus musculus -> Turnip mosaic virus

**Expected/actual: NO_EVIDENCE_FOUND / NO_EVIDENCE_FOUND; label correct.** Cross-group decoy, tests exposure/immunization, recombinant viral components, VLPs and therapeutic delivery versus infection by the target virus.

Seven no-support and one mention-only records. P2 / PMID36159839 discusses TuMV-derived VLP allergy treatments; P3 / PMID31658770 is mention-only for plant-virus nanoparticles in an in-vivo mouse model. P4 / PMID1343828 extracts an `INFECTION_OF` label for splenocytes from mice immunized with TuMV, but has no verified target infection edge. P5 / PMID42100486 discusses a TuMV protease delivered with another virus (AAV); the delivery vector is not evidence that TuMV infects mice. Plant infection papers P1 / PMID39091990 and P7 / PMID40241733 also fail target-host binding.

This is the strongest guard against fixing positives by accepting immunization, viral protein expression, any experimental exposure, or study-wide keywords as infection. MEDIUM confidence remains uncalibrated and search completeness unknown.

### 12. Aedes aegypti -> Mouse hepatitis virus

**Expected/actual: NO_EVIDENCE_FOUND / NO_EVIDENCE_FOUND; label correct.** Decoy within the broad ANIMAL group, tests mosquito/mouse/virus role separation where a coarse plausibility prior cannot help.

All eight records are no-support. P0 / PMID39074957 describes dengue-infected mosquitoes and mice; extracted Mouse hepatitis comparison content does not establish the target edge. P2 / PMID42043253, a virology meeting report, separately discusses Aedes and naturally mouse-infecting MHV-Y; the mouse-virus passage cannot bind the mosquito. P3 / PMID34957305 discusses dengue transmission and in-vitro antiviral activity against MHV and DENV-2. P7 / PMID33525547 describes mouse-coronavirus replication organelles but fails target-host binding.

The cached prior is `BIOLOGICALLY_CONSISTENT` because known hosts and the target share ANIMAL. The correct rejection illustrates why that prior cannot establish the interaction. Whole-document host binding or loose common-name/acronym repairs would endanger this case.

## Failure taxonomy and causal certainty

| Incorrect case | Primary failure categories | Root decision failure | Responsible files |
|---|---|---|---|
| 2 Arabidopsis/TuMV | extraction failure; entity binding error; other-host contamination; judge aggregation failure; other | Short/mislabeled assertions rejected; multi-host evidence reduced to other-host; all remaining classes treated as absence. | `evidence_agent.py`, `relationship_language.py`, `judge_agent.py` |
| 3 mouse/hepatitis | extraction failure; entity binding error; judge aggregation failure; other | Reference retrieved but exact quote/type mismatch; acronym-only/reconstructed quotes; failed parse hidden by non-support. | Same three files |
| 5 Lampyris/Hubei 31 | entity binding error; co-mention false positive; comparison-virus contamination; sequence-similarity error; confidence/calibration failure | Study-wide context validates a comparison passage as host evidence; exact support wins. | `evidence_agent.py`, `judge_agent.py` |
| 6 Lampyris/Hubei 51 | entity binding error; extraction failure | Comparison source wrongly points to target itself; chain fails. | `evidence_agent.py` |
| 7 Lampyris/chuvirus 3 | entity binding error; comparison-virus contamination; other-host contamination; extraction failure; judge aggregation failure | Reversed roles and missing local host binding; failed extractions become non-support. | `evidence_agent.py`, `judge_agent.py` |
| 8 Lampyris/toti 16 | entity binding error; comparison-virus contamination; other-host contamination; sequence-similarity error; extraction failure | Wrong host-edge virus and missing host linkage lose valid related chain while other-host branch fires. | `evidence_agent.py` |

No historical HTTP/query ledger exists, so **API/retrieval failure cannot be assigned as a proven root cause**. Search coverage is inadequate to audit fully, but key papers are demonstrably selected in cases 3 and 5–8. Generic taxonomy alias error and environmental-source error are risks, not demonstrated universal causes of these wrong labels. This report distinguishes identifiable verifier branches from unknown upstream causes of model output. “Exact root cause” can be established at that recorded decision boundary; it cannot honestly be extended to missing prompts or network events.

The current `_flexible_name_matches` in `evidence_agent.py:1236` contains double-escaped raw regex components, including `r"[\\W_]+"`, so ordinary spaces/hyphens are not the intended separators. This is a concrete current-code defect identified during the architecture audit, but **not proof that it caused the saved v06 outputs**. Those outputs include grounded windows and lack source-version authentication. A real Python regression test is required before changing it. No regex repair was implemented here.

## Minimal ordered fixes proposed for approval

These are narrow changes to existing modules, not a V2 architecture migration. Do not weaken rejection checks to force the related or negative labels. Run all twelve cases plus focused misleading-evidence fixtures at each step; compare both labels and accepted evidence.

| Order | Small change proposed | Cases likely helped | Existing correct cases at risk / required guard |
|---:|---|---|---|
| 1 | Add live diagnostic capture around the existing evaluator: exact queries/source statuses, raw/selected papers, host aliases, full extraction input/output and rescue attempts, effective model/tokenizer/config hashes, cache files before/after and outcomes. Use new run paths; replay saved inputs without network. | Enables trustworthy attribution and future comparison for every case. | All cases: instrumentation must not consume extra RNG, modify prompts, reorder queries/papers, or mutate cached results. Compare diagnostic-on/off outputs on identical frozen input. |
| 2 | Fix and unit-test flexible-name regex; record original offsets and complete sentence/paragraph spans. For an exact but incomplete quote, request a coherent assertion span rather than treating substring validity as sufficient. Resolve abbreviations only from local definitions. | 2, 3, 7, 8; allows proper linked evidence instead of literal fragments. | 9, 10, 12: larger windows must not borrow another host/virus. 11: exposure/protein/VLP evidence must not become infection. 1: preserve valid multi-passage discovery. |
| 3 | Validate host-edge roles against the predicate and comparison context. Reject a comparison-only virus span as a host assertion even with a permitted host-type label. Require explicit shared study/sample linkage for cross-passage support. Trigger bounded rescue on role conflicts, not just empty fields or invalid type labels. | 5 primarily; 7, 8; evidence quality in 4. | 1 and 4: their current `STUDY_CONTEXT` positives may be rejected. Recover grounded discovery/clinical evidence; do not preserve unsafe acceptance for score. 9 and 12 must remain resistant to wrong-host promotion. |
| 4 | Ground comparison source and target separately; use explicit document-local acronym/coreference links. Detect self-comparisons and source mismatches for targeted re-extraction. Keep source equality checks; never auto-swap from benchmark expectations. | 6, 7, 8, and produces the correct related chain for 5 after rejecting false exact support. | 1: comparison repair must not displace legitimate exact support. 9, 10, 12: no related chain without a genuine target-host edge. |
| 5 | Correct relationship type selection and minimal semantic validation for complete assertions (infection versus natural infection, occurrence/prevalence versus detection assay). Re-extract unsupported type labels instead of adding broad keyword OR clauses. Handle multi-host statements as separate candidate assertions. | 2, 3; cleaner evidence for 4. | 11: immunization, antigen exposure and AAV delivery remain non-support for TuMV infection. 4: in-vitro/no-infection models remain scoped. 1: sequence discovery remains distinct from demonstrated replication. |
| 6 | Preserve material failed/ambiguous extraction as `UNCLEAR` under the existing label vocabulary until retry or independently justified exclusion succeeds. Make the judge consume extraction/search health; keep credible exact support precedence where valid. Do not translate missing evidence to novelty. | 3, 7; prevents overstated absence. | 10 currently matches its label despite a failed extraction and may temporarily become `UNCLEAR`. Resolve that failure honestly; do not suppress it. Other negative cases may change if live retrieval fails. |
| 7 | Give reasons/confidence their actual basis: invalid/mixed-scope support cannot earn HIGH; unresolved material coverage/extraction cannot earn unqualified absence confidence. Avoid claiming calibrated probability. Only add targeted retrieval escalation if corrected extraction still lacks adequate evidence. | 5, 2, 3, 7; all diagnostics. | 1 and 4 may lose HIGH confidence until support quality is verified; negative labels must not rely on an arbitrary source-call count. |

Files implicated: `evidence_agent.py` for steps 2–5 and failure outcomes; `bioresearch_env/relationship_language.py` for narrowly tested type validation; `judge_agent.py` for step 6–7; `evaluate_env.py`, `bioresearch_env/env.py` and search adapters only for step 1 instrumentation. Taxonomy alias lookup should not be broadly rewritten to solve these recorded failures. Retrieval expansion is conditional, not an excuse to conceal extraction defects.

Step 1 is still needed for full **live-run reproducibility**; the delivered artifact exporter cannot recover historical information that was never recorded. It is now possible to reproduce the historical diagnostic baseline locally and verify later exports against preserved evidence. When the model/runtime and source assets are available, capture a new run before applying any proposed semantic fixes. A fresh score may differ because the current source is not authenticated to the old CSV.

## Acceptance criteria and stopping point

Keep all twelve expected labels unchanged. The eventual 12/12 requirement must be accompanied by correct paper/host/virus/relationship attribution, not just labels. In particular, a 12/12 score obtained by accepting case 5's comparison-only span or case 4's mixed-scope support is not a sound baseline.

Before approving any semantic fix, require focused tests for: comparison promotion; swapped comparison endpoints; short-but-verbatim assertions; defined versus ambiguous acronyms; negation; multi-host sentences; immunization/VLP/protein exposure; other-organism sample context; failed extraction; and valid cross-passage discovery. The exporter tests are data-integrity tests, not biological unit tests, and do not claim to validate proposed fixes.

This phase stops at diagnosis. Only diagnostic/export/test files and the new JSON snapshot were added. Original benchmarks, historical results/caches, project instructions, architecture review and production source remain unchanged. The proposed fixes have not been implemented.
