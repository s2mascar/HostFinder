# Update — correctness iteration 2 (2026-09-26)

The authoritative latest **live baseline is 0/12, all UNCLEAR**, supplied by the user for Qwen3-8B on Nibi. The old live CSV/log are not present locally. The original 42 tests passed on Nibi before that run. This iteration now has **77 passing local Python tests**, with all original tests preserved.

**Post-change live benchmark: PENDING.** All 12 final classifications and correctness judgments require the new Nibi run. No repaired live label or evidence edge is claimed. Zero prior false KNOWN/related calls reflected abstention; it was not a passing benchmark.

| # | Host | Virus | Expected | Latest live before | New live after | Correct after |
|---|---|---|---|---|---|---|
| 1 | Lampyris noctiluca | Lampyris noctiluca partiti-like virus 1 | KNOWN | UNCLEAR | PENDING | Not measured |
| 2 | Arabidopsis thaliana | Turnip mosaic virus | KNOWN | UNCLEAR | PENDING | Not measured |
| 3 | Mus musculus | Mouse hepatitis virus | KNOWN | UNCLEAR | PENDING | Not measured |
| 4 | Homo sapiens | SARS-CoV-2 | KNOWN | UNCLEAR | PENDING | Not measured |
| 5 | Lampyris noctiluca | Hubei partiti-like virus 31 | POSSIBLY_KNOWN | UNCLEAR | PENDING | Not measured |
| 6 | Lampyris noctiluca | Hubei partiti-like virus 51 | POSSIBLY_KNOWN | UNCLEAR | PENDING | Not measured |
| 7 | Lampyris noctiluca | Hubei chuvirus-like virus 3 | POSSIBLY_KNOWN | UNCLEAR | PENDING | Not measured |
| 8 | Lampyris noctiluca | Hubei toti-like virus 16 | POSSIBLY_KNOWN | UNCLEAR | PENDING | Not measured |
| 9 | Tribolium castaneum | Hubei partiti-like virus 31 | NO_EVIDENCE_FOUND | UNCLEAR | PENDING | Not measured |
| 10 | Arabidopsis thaliana | Mouse hepatitis virus | NO_EVIDENCE_FOUND | UNCLEAR | PENDING | Not measured |
| 11 | Mus musculus | Turnip mosaic virus | NO_EVIDENCE_FOUND | UNCLEAR | PENDING | Not measured |
| 12 | Aedes aegypti | Mouse hepatitis virus | NO_EVIDENCE_FOUND | UNCLEAR | PENDING | Not measured |

Current changes and tests are documented in `CORRECTNESS_IMPLEMENTATION_REPORT.md`. Run `run_correctness_benchmark.sh` on Nibi and return the entire generated run directory. Full source, extraction and request traces are now captured. `CORRECTNESS_FINAL_REPORT.md` explicitly marks the checkpoint blocked, not frozen.

Latest offline quotation replay: `results/correctness_iteration2_offline.json`; 0/12 due unavailable identity/coverage, not a live substitute. The original historical fixture remains unchanged. Machine-readable user summary: `results/nibi_baseline_user_summary.json`.

---

The remainder below preserves the previous iteration's historical replay report. Its 'live score not measured' statement refers to that earlier report, not the newly supplied live baseline above.

# Correctness benchmark after implementation

## Result and limits

Historical benchmark: **6/12**. Conservative preserved-quotation replay: **0/12** against unchanged legacy labels. **Live score: not measured. The 12/12 regression gate is not satisfied.**

All 12 replay decisions abstain because the artifact does not record independently validated taxonomy/rank or complete search outcomes. This is a data-completeness limitation, not a demonstration of improved benchmark accuracy. The two historically correct KNOWN calls and four correct negative calls are not claimed preserved by this replay.

The local bundled Python is 3.12.14. It lacks requests, torch and transformers; this workspace lacks models/Qwen3-8B and agent_env. No live retrieval or model inference was attempted. No model, dependency or paper downloads were performed. Historical model outputs were not repaired, relabeled or treated as fresh source documents.

| Metric | Replay count |
|---|---:|
| Correct legacy labels | 0/12 |
| False KNOWN | 0 |
| False NOVEL_CANDIDATE | 0 |
| False BIOLOGICALLY_IMPLAUSIBLE | 0 |
| INSUFFICIENT_EVIDENCE | 12 |

Zero false positive calls here reflects abstention, not solved sensitivity. The legacy benchmark has KNOWN, POSSIBLY_KNOWN and NO_EVIDENCE_FOUND labels. Neither POSSIBLY_KNOWN nor NO_EVIDENCE_FOUND is mapped to novelty; the canonical classification remains INSUFFICIENT_EVIDENCE without exact support.

## All 12 before/after results

| # | Host | Virus | Expected | Before | After (legacy / canonical) | Correct |
|---|---|---|---|---|---|---|---|
| 1 | Lampyris noctiluca | Lampyris noctiluca partiti-like virus 1 | KNOWN | KNOWN | UNCLEAR / INSUFFICIENT_EVIDENCE | False |
| 2 | Arabidopsis thaliana | Turnip mosaic virus | KNOWN | NO_EVIDENCE_FOUND | UNCLEAR / INSUFFICIENT_EVIDENCE | False |
| 3 | Mus musculus | Mouse hepatitis virus | KNOWN | NO_EVIDENCE_FOUND | UNCLEAR / INSUFFICIENT_EVIDENCE | False |
| 4 | Homo sapiens | SARS-CoV-2 | KNOWN | KNOWN | UNCLEAR / INSUFFICIENT_EVIDENCE | False |
| 5 | Lampyris noctiluca | Hubei partiti-like virus 31 | POSSIBLY_KNOWN | KNOWN | UNCLEAR / INSUFFICIENT_EVIDENCE | False |
| 6 | Lampyris noctiluca | Hubei partiti-like virus 51 | POSSIBLY_KNOWN | NO_EVIDENCE_FOUND | UNCLEAR / INSUFFICIENT_EVIDENCE | False |
| 7 | Lampyris noctiluca | Hubei chuvirus-like virus 3 | POSSIBLY_KNOWN | NO_EVIDENCE_FOUND | UNCLEAR / INSUFFICIENT_EVIDENCE | False |
| 8 | Lampyris noctiluca | Hubei toti-like virus 16 | POSSIBLY_KNOWN | NO_EVIDENCE_FOUND | UNCLEAR / INSUFFICIENT_EVIDENCE | False |
| 9 | Tribolium castaneum | Hubei partiti-like virus 31 | NO_EVIDENCE_FOUND | NO_EVIDENCE_FOUND | UNCLEAR / INSUFFICIENT_EVIDENCE | False |
| 10 | Arabidopsis thaliana | Mouse hepatitis virus | NO_EVIDENCE_FOUND | NO_EVIDENCE_FOUND | UNCLEAR / INSUFFICIENT_EVIDENCE | False |
| 11 | Mus musculus | Turnip mosaic virus | NO_EVIDENCE_FOUND | NO_EVIDENCE_FOUND | UNCLEAR / INSUFFICIENT_EVIDENCE | False |
| 12 | Aedes aegypti | Mouse hepatitis virus | NO_EVIDENCE_FOUND | NO_EVIDENCE_FOUND | UNCLEAR / INSUFFICIENT_EVIDENCE | False |

## Case evidence and reasoning

The following records apply the new verifier to the original saved fields. Each source fragment is admitted only if its historical quotation-verification flag was true. Missing host aliases are not guessed. No comparison endpoints are swapped. Paper-level results are diagnostic re-verifications of saved quotations, not independent validations of the papers.

### 1. Lampyris noctiluca → Lampyris noctiluca partiti-like virus 1

Expected **KNOWN**; before **KNOWN**; after **INSUFFICIENT_EVIDENCE** (legacy UNCLEAR); correct: **False**.

Host or virus taxonomy resolution is unresolved or unrecorded. Complete per-source search coverage is also unavailable.

- **P0 / 31900852** — Identification and characterisation of common glow-worm RNA viruses. `UNCLEAR`; SELF_COMPARISON; targeted re-extraction required; endpoints were not swapped. Proposed edge: `Lampyris noctiluca → IDENTIFIED_IN → Lampyris noctiluca partiti-like virus 1`. Scope: `AMBIGUOUS`. Saved quotation: "Lampyris noctiluca partiti-like virus 1 (LnoPLV1) genome was 1462 nt long and codes for a protein of 377 aa. According to Blastp search, the protein is similar to RdRP of Hubei partiti-like virus 31 (1,923,038, 93% coverage and 72% identity)."

### 2. Arabidopsis thaliana → Turnip mosaic virus

Expected **KNOWN**; before **NO_EVIDENCE_FOUND**; after **INSUFFICIENT_EVIDENCE** (legacy UNCLEAR); correct: **False**.

Host or virus taxonomy resolution is unresolved or unrecorded. Complete per-source search coverage is also unavailable.

- **P0 / P0** — Turnip mosaic virus drives selective filtering and community reassembly in the <i>Arabidopsis thaliana</i> root microbiome in a genotype-specific manner. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Arabidopsis thaliana → INFECTION_OF → Turnip mosaic virus`. Scope: `AMBIGUOUS`. Saved quotation: "Turnip mosaic virus (TuMV) alters root-associated bacterial and fungal communities"
- **P1 / 42265533** — Turnip mosaic virus utilizes the lipid droplet biogenesis machinery to facilitate its propagation in plants.. `UNCLEAR`; SELF_COMPARISON; targeted re-extraction required; endpoints were not swapped. Proposed edge: `Arabidopsis thaliana → INFECTION_OF → Turnip mosaic virus`. Scope: `AMBIGUOUS`. Saved quotation: "turnip mosaic virus (TuMV) infection"
- **P2 / 42372481** — A viral infection reshapes Arabidopsis water management via root hydraulics, aquaporin downregulation and osmotic adjustment.. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Arabidopsis thaliana → INFECTION_OF → Turnip mosaic virus`. Scope: `AMBIGUOUS`. Saved quotation: "Using a hydroponic Arabidopsis thaliana-Turnip mosaic virus (TuMV) pathosystem"
- **P3 / P3** — Evolution of virulence of a plant RNA virus in age-diverse host populations. `UNCLEAR`; SELF_COMPARISON; targeted re-extraction required; endpoints were not swapped. Proposed edge: `Arabidopsis thaliana → EXPERIMENTAL_INFECTION → Turnip mosaic virus`. Scope: `AMBIGUOUS`. Saved quotation: "turnip mosaic virus (TuMV)"
- **P4 / 42710218** — Turnip mosaic virus alters phosphorus metabolism and shoot-root allocation without resource competition.. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Arabidopsis thaliana → INFECTION_OF → Turnip mosaic virus`. Scope: `AMBIGUOUS`. Saved quotation: "Turnip mosaic virus (TuMV) drawed significant P internal pools leading to P competition"
- **P5 / 42347223** — Contrasting Effects of Tagging Turnip Mosaic Virus Proteins.. `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Arabidopsis thaliana → MENTION_ONLY → Turnip mosaic virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "TuMV-GFP"
- **P6 / 42671063** — Evolution of virulence of a plant RNA virus in developmental stage-structured host populations.. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Arabidopsis thaliana → NATURAL_INFECTION → Turnip mosaic virus`. Scope: `AMBIGUOUS`. Saved quotation: "We used populations of A. thaliana that differed in demographic composition (Fig. 1) to evaluate the evolutionary dynamics of their natural parasite, TuMV (species Potyvirus rapae, genus Potyvirus, family Potyviridae)."
- **P7 / 25018765** — Dominant resistance against plant viruses. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Arabidopsis thaliana → DIFFERENT_HOST → Turnip mosaic virus`. Scope: `AMBIGUOUS`. Saved quotation: "Brassica campestris BcTuR3 TIR-NBS-LRR TuMV [Turnip mosaic virus] Potyvirus Unknown 17, 18"

### 3. Mus musculus → Mouse hepatitis virus

Expected **KNOWN**; before **NO_EVIDENCE_FOUND**; after **INSUFFICIENT_EVIDENCE** (legacy UNCLEAR); correct: **False**.

Host or virus taxonomy resolution is unresolved or unrecorded. Complete per-source search coverage is also unavailable.

- **P0 / 8677422** — [Mouse hepatitis virus].. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Mus musculus → NATURAL_INFECTION → Mouse hepatitis virus`. Scope: `AMBIGUOUS`. Saved quotation: "Mouse hepatitis virus (MHV), the coronavirus of the mouse (mus musculus),"
- **P1 / 1719235** — Cloning of the mouse hepatitis virus (MHV) receptor: expression in human and hamster cell lines confers susceptibility to MHV.. `UNCLEAR`; Extraction failed; absence of an extracted edge is not non-support. Proposed edge: ` → UNCLEAR → `. Scope: `AMBIGUOUS`. Saved quotation: "No usable quotation saved."
- **P2 / 41860231** — A novel five-plex digital PCR assay for the simultaneous detection of murine pathogens: Sendai virus, reovirus, mouse parvoviruses, pneumonia virus of mice, and mouse hepatitis virus.. `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Mus musculus → MENTION_ONLY → Mouse hepatitis virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "five experimental animal viruses-Sendai virus (SeV), reovirus 3 (REO3), mouse parvovirus (MPV), pneumonia virus of mice (PVM), and mouse hepatitis virus (MHV)-that infect laboratory mice and rats relatively easily"
- **P3 / 15953191** — Mouse hepatitis virus 3 binding to macrophages correlates with resistance to experimental infection.. `UNCLEAR`; SELF_COMPARISON; targeted re-extraction required; endpoints were not swapped. Proposed edge: `Mus musculus → NATURAL_INFECTION → Mouse hepatitis virus`. Scope: `AMBIGUOUS`. Saved quotation: "MHV3"
- **P4 / 15638127** — Arginine metabolism during macrophage autocrine activation and infection with mouse hepatitis virus 3.. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Mus musculus → NATURAL_INFECTION → Mouse hepatitis virus 3`. Scope: `AMBIGUOUS`. Saved quotation: "anti-mouse hepatitis virus 3 (MHV3) state, infection with mouse hepatitis virus 3 (MHV3)"
- **P5 / 10803365** — Sequence analysis of major structural proteins of newly isolated mouse hepatitis virus.. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Mus musculus → INFECTION_OF → Mouse hepatitis virus`. Scope: `AMBIGUOUS`. Saved quotation: "mouse hepatitis virus (MHV) epidemic"
- **P6 / 25502736** — Microbiological survey of mice (Mus musculus) purchased from commercial pet shops in Kanagawa and Tokyo, Japan.. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Mus musculus → DETECTED_IN → mouse hepatitis virus`. Scope: `AMBIGUOUS`. Saved quotation: "mouse hepatitis virus (12 mice; 42.8%)"
- **P7 / 21508117** — Infectious microorganisms in mice (Mus musculus) purchased from commercial pet shops in Germany.. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Mus musculus → DETECTED_IN → Mouse hepatitis virus`. Scope: `AMBIGUOUS`. Saved quotation: "We found a number of microorganisms in these pet shop mice, the most prevalent of which were Helicobacter species (92.9%), mouse parvovirus (89.3%), mouse hepatitis virus (82.7%), Pasteurella pneumotropica (71.4%) and Syphacia species (57.1%)."

### 4. Homo sapiens → SARS-CoV-2

Expected **KNOWN**; before **KNOWN**; after **INSUFFICIENT_EVIDENCE** (legacy UNCLEAR); correct: **False**.

Host or virus taxonomy resolution is unresolved or unrecorded. Complete per-source search coverage is also unavailable.

- **P0 / 42743972** — PhoSARte: identification of SARS-CoV-2 phosphorylation sites using contrastive learning and protein language models.. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Homo sapiens → DETECTED_IN → SARS-CoV-2`. Scope: `AMBIGUOUS`. Saved quotation: "SARS-CoV-2 hijacks host kinase signaling networks to enhance replication and suppress antiviral responses [3]. This virus-driven reprogramming is extensive, with approximately 70 phosphorylation sites identified in viral proteins and more than 15 000 phosphorylation events mapped in host proteins during infection [4]."
- **P1 / 40354427** — Infectious potential and circulation of SARS-CoV-2 in wild rats.. `UNCLEAR`; SELF_COMPARISON; targeted re-extraction required; endpoints were not swapped. Proposed edge: `Rattus norvegicus → NO_RELATION → SARS-CoV-2`. Scope: `AMBIGUOUS`. Saved quotation: "We studied the circulation of SARS-CoV-2 in wild Rattus norvegicus (n = 401) captured in urban areas and sewage systems of several French cities."
- **P2 / 36317085** — Molecular Mimicry of SARS-CoV-2 Spike Protein in the Nervous System: A Bioinformatics Approach.. `UNCLEAR`; SELF_COMPARISON; targeted re-extraction required; endpoints were not swapped. Proposed edge: `Homo sapiens → MENTION_ONLY → SARS-CoV-2`. Scope: `AMBIGUOUS`. Saved quotation: "severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)"
- **P3 / 42212595** — FKBP8 inhibits influenza a virus infection by degrading viral M2 protein in lysosomes.. `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Homo sapiens → MENTION_ONLY → Influenza A virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "Influenza A virus (IAV) remains a major threat to global public health, causing seasonal epidemics and occasional pandemics with significant morbidity and mortality."
- **P4 / 40268235** — BEAGLE 2.0: A Web Server for RNA Secondary Structure Similarity Detection Leveraging SHAPE-directed RNA Structure Determination.. `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Homo sapiens → MENTION_ONLY → SARS-CoV-2`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "various viruses, including SARS-CoV-2"
- **P5 / 42751539** — Slowdown of synonymous substitution rate preceding the emergence of multiple SARS-CoV-2 variants.. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Homo sapiens → NATURAL_INFECTION → SARS-CoV-2`. Scope: `AMBIGUOUS`. Saved quotation: "The stem branches of variants have the same pattern of accelerated nonsynonymous substitutions, and so persistent infections have been proposed as a possible source of variants."
- **P6 / 42645657** — SARS-CoV-2 disrupts the integrity of a human blood-brain barrier model in the absence of infection.. `UNCLEAR`; SELF_COMPARISON; targeted re-extraction required; endpoints were not swapped. Proposed edge: `Homo sapiens → EXPERIMENTAL_INFECTION → SARS-CoV-2`. Scope: `AMBIGUOUS`. Saved quotation: "exposure of the model to SARS-CoV-2 or spike protein"
- **P7 / 34146538** — Peptides of H. sapiens and P. falciparum that are predicted to bind strongly to HLA-A*24:02 and homologous to a SARS-CoV-2 peptide.. `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Homo sapiens → MENTION_ONLY → SARS-CoV-2`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "SARS-CoV-2 peptide with single letter amino acid code CFLGYFCTCYFGLFC has the highest identity to P. vivax. Its YFCTCYFGLF part is predicted to bind strongly to HLA-A*24:02."

### 5. Lampyris noctiluca → Hubei partiti-like virus 31

Expected **POSSIBLY_KNOWN**; before **KNOWN**; after **INSUFFICIENT_EVIDENCE** (legacy UNCLEAR); correct: **False**.

Host or virus taxonomy resolution is unresolved or unrecorded. Complete per-source search coverage is also unavailable.

- **P0 / 31900852** — Identification and characterisation of common glow-worm RNA viruses. `UNCLEAR`; COMPARISON_SOURCE_MISMATCH; targeted re-extraction required; endpoints were not swapped. Proposed edge: `Lampyris noctiluca → SEQUENCED_FROM → Hubei partiti-like virus 31`. Scope: `AMBIGUOUS`. Saved quotation: "According to Blastp search, the protein is similar to RdRP of Hubei partiti-like virus 31 (1,923,038, 93% coverage and 72% identity)."

### 6. Lampyris noctiluca → Hubei partiti-like virus 51

Expected **POSSIBLY_KNOWN**; before **NO_EVIDENCE_FOUND**; after **INSUFFICIENT_EVIDENCE** (legacy UNCLEAR); correct: **False**.

Host or virus taxonomy resolution is unresolved or unrecorded. Complete per-source search coverage is also unavailable.

- **P0 / 31900852** — Identification and characterisation of common glow-worm RNA viruses. `UNCLEAR`; SELF_COMPARISON; targeted re-extraction required; endpoints were not swapped. Proposed edge: `Lampyris noctiluca → IDENTIFIED_IN → Lampyris noctiluca partiti-like virus 2`. Scope: `AMBIGUOUS`. Saved quotation: "Lampyris noctiluca partiti-like virus 2 (LnoPLV2) genome segment was 1461 nt long and encoded a 436 aa protein."

### 7. Lampyris noctiluca → Hubei chuvirus-like virus 3

Expected **POSSIBLY_KNOWN**; before **NO_EVIDENCE_FOUND**; after **INSUFFICIENT_EVIDENCE** (legacy UNCLEAR); correct: **False**.

Host or virus taxonomy resolution is unresolved or unrecorded. Complete per-source search coverage is also unavailable.

- **P0 / 31900852** — Identification and characterisation of common glow-worm RNA viruses.. `NO_SUPPORT`; The cited assertion is background_mention. Proposed edge: `Lampyris noctiluca → SEQUENCED_FROM → Hubei chuvirus-like virus 3`. Scope: `BACKGROUND_MENTION`. Saved quotation: "According to phylogenetic analysis, LnoCLV1 was most similar to Hubei chuvirus-like virus 3 (1,922,858) (Online Resource 9), which has been isolated from Odonata mix and has a monopartite genome [16]."
- **P1 / 40512168** — Annual (2024) taxonomic update of RNA-directed RNA polymerase-encoding negative-sense RNA viruses (realm Riboviria: kingdom Orthornavirae: phylum Negarnaviricota). `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Lampyris noctiluca → MENTION_ONLY → Hubei chuvirus-like virus 3`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "Scarabeuvirus hubeiense Húběi chuvirus-like virus 1 (HbCLV1)"
- **P2 / 34463877** — 2021 TAXONOMIC UPDATE OF PHYLUM Negarnaviricota (Riboviria: Orthornavirae), INCLUDING THE LARGE ORDERS Bunyavirales AND Mononegavirales. `UNCLEAR`; Extraction failed; absence of an extracted edge is not non-support. Proposed edge: ` → UNCLEAR → `. Scope: `AMBIGUOUS`. Saved quotation: "No usable quotation saved."
- **P3 / 37622664** — Annual (2023) taxonomic update of RNA-directed RNA polymerase-encoding negative-sense RNA viruses (realm Riboviria: kingdom Orthornavirae: phylum Negarnaviricota). `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Lampyris noctiluca → MENTION_ONLY → Hubei chuvirus-like virus 3`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "Scarabeuvirus hubeiense Húběi chuvirus-like virus 1 (HbCLV1)"
- **P4 / 36437428** — 2022 TAXONOMIC UPDATE OF PHYLUM Negarnaviricota (Riboviria: Orthornavirae), INCLUDING THE LARGE ORDERS Bunyavirales AND Mononegavirales. `UNCLEAR`; Extraction failed; absence of an extracted edge is not non-support. Proposed edge: ` → NO_RELATION → `. Scope: `AMBIGUOUS`. Saved quotation: "No usable quotation saved."
- **P5 / 35108077** — Jingchuvirales: a New Taxonomical Framework for a Rapidly Expanding Order of Unusual Monjiviricete Viruses Broadly Distributed among Arthropod Subphyla. `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Lampyris noctiluca → MENTION_ONLY → Hubei chuvirus-like virus 3`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "Hubei chuvirus-like virus 3"

### 8. Lampyris noctiluca → Hubei toti-like virus 16

Expected **POSSIBLY_KNOWN**; before **NO_EVIDENCE_FOUND**; after **INSUFFICIENT_EVIDENCE** (legacy UNCLEAR); correct: **False**.

Host or virus taxonomy resolution is unresolved or unrecorded. Complete per-source search coverage is also unavailable.

- **P0 / 31900852** — Identification and characterisation of common glow-worm RNA viruses. `UNCLEAR`; COMPARISON_SOURCE_MISMATCH; targeted re-extraction required; endpoints were not swapped. Proposed edge: `Lampyris noctiluca → SEQUENCED_FROM → Hubei toti-like virus 16`. Scope: `AMBIGUOUS`. Saved quotation: "According to Blastp search, the shorter ORF was similar to hypothetical protein 2 of Hubei toti-like virus 16 (99% coverage and 31% identity), whereas HHPred found no similar protein sequences. Hubei toti-like virus 16 is isolated form spiders [16]. According to phylogenetic analysis, LnoTLV1 was most similar to Hubei toti-like virus 16 and Beihai sea slater virus 3, isolated from wharf roach (1,922,659) (Online Resource 11)."
- **P1 / 35337056** — Transcriptomics Reveal Several Novel Viruses from Canegrubs (Coleoptera: Scarabaeidae) in Central Queensland, Australia. `NO_SUPPORT`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Lampyris noctiluca → SEQUENCE_SIMILARITY → Hubei toti-like virus 16`. Scope: `SEQUENCE_SIMILARITY_ONLY`. Saved quotation: "The first ORF is 3984 nt, and encodes a 1328 amino acids (aa) protein that shares 29.4% aa sequence identity (query coverage, 84%) with the hypothetical protein 2 of Hubei toti-like virus 16, which was isolated from spiders in China [48]."

### 9. Tribolium castaneum → Hubei partiti-like virus 31

Expected **NO_EVIDENCE_FOUND**; before **NO_EVIDENCE_FOUND**; after **INSUFFICIENT_EVIDENCE** (legacy UNCLEAR); correct: **False**.

Host or virus taxonomy resolution is unresolved or unrecorded. Complete per-source search coverage is also unavailable.

- **P0 / 31900852** — Identification and characterisation of common glow-worm RNA viruses. `NO_SUPPORT`; The cited assertion is background_mention. Proposed edge: `Tribolium castaneum → SEQUENCED_FROM → Hubei partiti-like virus 31`. Scope: `BACKGROUND_MENTION`. Saved quotation: "According to Blastp search, the protein is similar to RdRP of Hubei partiti-like virus 31 (1,923,038, 93% coverage and 72% identity)."

### 10. Arabidopsis thaliana → Mouse hepatitis virus

Expected **NO_EVIDENCE_FOUND**; before **NO_EVIDENCE_FOUND**; after **INSUFFICIENT_EVIDENCE** (legacy UNCLEAR); correct: **False**.

Host or virus taxonomy resolution is unresolved or unrecorded. Complete per-source search coverage is also unavailable.

- **P0 / 39605984** — Advances in the structure and function of the nucleolar protein fibrillarin. `UNCLEAR`; Extraction failed; absence of an extracted edge is not non-support. Proposed edge: ` → UNCLEAR → `. Scope: `AMBIGUOUS`. Saved quotation: "No usable quotation saved."
- **P1 / 39525088** — Tiny but mighty: Diverse functions of uORFs that regulate gene expression. `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Arabidopsis thaliana → MENTION_ONLY → Mouse hepatitis virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "The 5′-UTR of the coronavirus mouse hepatitis virus (MHV) contains four stem-loop (SL) structures that regulate the expression of the downstream ORF1"
- **P2 / 39364891** — The role of structure in regulatory RNA elements. `NO_SUPPORT`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Arabidopsis thaliana → NO_RELATION → Mouse hepatitis virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "No usable quotation saved."
- **P3 / 38599165** — The NLR family of innate immune and cell death sensors. `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Arabidopsis thaliana → MENTION_ONLY → Mouse hepatitis virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "During murine coronavirus mouse hepatitis virus (MHV) infection, NLRP6 binds to dsRNA and undergoes liquid–liquid phase separation dependent on its disordered poly-lysine sequence (K350-354) to drive inflammasome formation and innate immune signaling."
- **P4 / 40828816** — Mobilization of nuclear antiviral factors by exportin XPO1 via the actin network inhibits RNA virus replication. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Arabidopsis thaliana → IDENTIFIED_IN → Tomato bushy stunt virus`. Scope: `AMBIGUOUS`. Saved quotation: "Previous genome- and proteome-wide approaches have identified numerous nuclear proteins, including restriction factors that affect replication of tomato bushy stunt virus (TBSV)."
- **P5 / 36928641** — Current research on viral proteins that interact with fibrillarin.. `NO_SUPPORT`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Arabidopsis thaliana → NO_RELATION → `. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "No usable quotation saved."
- **P6 / 32553580** — An 'Arms Race' between the Nonsense-mediated mRNA Decay Pathway and Viral Infections.. `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Arabidopsis thaliana → MENTION_ONLY → Mouse hepatitis virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "In this review we highlight the reciprocal interactions between the host NMD pathway and viral pathogens, which have shaped both the host antiviral defense and viral pathogenesis."
- **P7 / 32699064** — Viral subversion of nonsense-mediated mRNA decay.. `NO_SUPPORT`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Arabidopsis thaliana → NO_RELATION → Mouse hepatitis virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "The target virus is associated with a host in the ANIMAL broad group."

### 11. Mus musculus → Turnip mosaic virus

Expected **NO_EVIDENCE_FOUND**; before **NO_EVIDENCE_FOUND**; after **INSUFFICIENT_EVIDENCE** (legacy UNCLEAR); correct: **False**.

Host or virus taxonomy resolution is unresolved or unrecorded. Complete per-source search coverage is also unavailable.

- **P0 / 21187975** — The nuclear inclusion a (NIa) protease of turnip mosaic virus (TuMV) cleaves amyloid-β.. `UNCLEAR`; SELF_COMPARISON; targeted re-extraction required; endpoints were not swapped. Proposed edge: `Mus musculus → EXPERIMENTAL_INFECTION → Turnip mosaic virus`. Scope: `AMBIGUOUS`. Saved quotation: "lentiviral-mediated expression of NIa in APP(sw)/PS1 transgenic mice"
- **P1 / 39091990** — Adaptation of turnip mosaic virus to <i>Arabidopsis thaliana</i> involves rewiring of VPg-host proteome interactions.. `UNCLEAR`; SELF_COMPARISON; targeted re-extraction required; endpoints were not swapped. Proposed edge: `Mus musculus → MENTION_ONLY → Turnip mosaic virus`. Scope: `AMBIGUOUS`. Saved quotation: "The outcome of a viral infection depends on a complex interplay between the host physiology and the virus, mediated through numerous protein-protein interactions."
- **P2 / 36159839** — Suitability of potyviral recombinant virus-like particles bearing a complete food allergen for immunotherapy vaccines.. `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Mus musculus → MENTION_ONLY → Turnip mosaic virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "VLPs derived from Turnip mosaic virus (TuMV)"
- **P3 / 31658770** — Elongated Flexuous Plant Virus-Derived Nanoparticles Functionalized for Autoantibody Detection.. `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Mus musculus → MENTION_ONLY → Turnip mosaic virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "Nanoparticles derived from the elongated flexuous capsids of Turnip mosaic virus (TuMV)"
- **P4 / 1343828** — Establishment of hybridoma cell line secreting specific monoclonal antibodies against turnip mosaic virus and analysis of properties of the McAb.. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Mus musculus → INFECTION_OF → Turnip mosaic virus`. Scope: `AMBIGUOUS`. Saved quotation: "immunized by TuMV"
- **P5 / 42100486** — Secretory form of viral protease NIa ameliorates amyloid-β pathology and cognitive deficits in a mouse model of Alzheimer’s disease. `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Mus musculus → MENTION_ONLY → Turnip mosaic virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "Nuclear inclusion a (NIa) is a cytosolic protease of turnip mosaic virus (TuMV) from the Potyviridae family, where it processes the viral polyprotein during maturation"
- **P6 / 40283196** — Transient Expression of Hen Egg White Lysozyme (EWL) in &lt;i&gt;Nicotiana benthamiana&lt;/i&gt; Influences Plant Pathogen Infection.. `NO_SUPPORT`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Mus musculus → NO_RELATION → Turnip mosaic virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "No usable quotation saved."
- **P7 / 40241733** — Resveratrol synthase homologs participate in infection of &lt;i&gt;Nicotiana benthamiana&lt;/i&gt; by pathogenic plant viruses and fungi.. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Nicotiana benthamiana → ASSOCIATED_WITH → Turnip mosaic virus`. Scope: `AMBIGUOUS`. Saved quotation: "The results showed that RS expression in plants significantly contributed to infection by turnip mosaic virus (TuMV) and slightly contributed to viral infection of tobacco mosaic virus (TMV)."

### 12. Aedes aegypti → Mouse hepatitis virus

Expected **NO_EVIDENCE_FOUND**; before **NO_EVIDENCE_FOUND**; after **INSUFFICIENT_EVIDENCE** (legacy UNCLEAR); correct: **False**.

Host or virus taxonomy resolution is unresolved or unrecorded. Complete per-source search coverage is also unavailable.

- **P0 / 39074957** — Electropenetrography with Alternating Current Reveals In Situ Changes of <i>Aedes aegypti</i> Probing Behaviors Associated with Dengue Virus Infection.. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Aedes aegypti → INFECTION_OF → Dengue virus`. Scope: `AMBIGUOUS`. Saved quotation: "DENV-infected Aedes aegypti mosquitoes feeding on uninfected mice and uninfected A. aegypti feeding on DENV-infected mice"
- **P1 / 38086541** — AC-DC Electropenetrography as a Tool to Quantify Probing and Ingestion Behaviors of the Yellow Fever Mosquito (<i>Aedes aegypti</i>) on Mice in Biocontainment.. `NO_SUPPORT`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Aedes aegypti → NO_RELATION → Mouse hepatitis virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "No usable quotation saved."
- **P2 / 42043253** — 25th Annual Meeting of the Rocky Mountain Virology Association. `UNCLEAR`; Missing, ungrounded or syntactically ambiguous directed assertion. Proposed edge: `Aedes aegypti → NATURAL_INFECTION → Mouse hepatitis virus Y`. Scope: `AMBIGUOUS`. Saved quotation: "mouse hepatitis virus Y (MHV-Y), that naturally infects the intestinal tract of mice"
- **P3 / 34957305** — Production, Transmission, Pathogenesis, and Control of Dengue Virus: A Literature-Based Undivided Perspective. `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Aedes aegypti → MENTION_ONLY → Mouse hepatitis virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "The ethyl acetate fraction of H. cordata and quercetin showed in vitro activity against mouse hepatitis virus (MHV) and DENV-2 with IC50 0.98 and 125 μg/mL for MHV while 7.50 and 176.76 μg/mL for DENV-2 [157]."
- **P4 / 39513877** — Acute Chikungunya Infection Induces Vascular Dysfunction by Directly Disrupting Redox Signaling in Endothelial Cells.. `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Aedes aegypti → MENTION_ONLY → Chikungunya virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "Chikungunya virus (CHIKV) infection is characterized by febrile illness, severe joint pain, myalgia, and cardiovascular complications."
- **P5 / 41933731** — SARS-CoV-2 envelope protein mitochondrial localization reveals host metabolic disruption. `UNCLEAR`; COMPARISON_SOURCE_MISMATCH; targeted re-extraction required; endpoints were not swapped. Proposed edge: `Aedes aegypti → MENTION_ONLY → Mouse hepatitis virus`. Scope: `AMBIGUOUS`. Saved quotation: "mouse hepatitis virus, a coronavirus discovered decades before SARS-CoV-2"
- **P6 / 42334213** — Double trouble: how co- and superinfections shape viral dynamics and host responses. `MENTION_ONLY`; Proposed relationship is contextual, not host-virus support. Proposed edge: `Aedes aegypti → MENTION_ONLY → Mouse hepatitis virus`. Scope: `COMPARISON_CONTEXT_ONLY`. Saved quotation: "However, different observations were obtained in murine lung epithelial cells using the model coronavirus murine hepatitis virus strain 1 (MHV-1)."
- **P7 / 33525547** — Multiscale Electron Microscopy for the Study of Viral Replication Organelles. `NO_SUPPORT`; The cited assertion is negative_or_absent_infection. Proposed edge: `Aedes aegypti → STUDY_ASSOCIATION → Mouse hepatitis virus`. Scope: `NEGATIVE_OR_ABSENT_INFECTION`. Saved quotation: "In order to avoid the use of fixatives and preserve the samples as close as possible to their native conditions, the study focused on the DMVs induced by mouse hepatitis virus (MHV), a well-established model coronavirus that does not pose serious biosafety constraints for cryo-EM sample preparation and imaging."

## Remaining work before accepting the regression gate

1. Run fresh extraction on the unchanged 12 pairs in the project model environment. Preserve full source text, raw primary/rescue responses, taxonomy resolution and query outcomes in the new sidecar.
2. Review the discovery paper's real specimen/study linkage and repair unresolved comparison roles through extraction, not endpoint substitution. The preserved introductory host description and virus genome-length fragments cannot establish that linkage.
3. Recover complete direct assertions for the positive plant and mouse cases. Verify primary human infection evidence rather than the old background/cell-line passage. The conservative assertion grammar will abstain on unsupported prose; extend it only with generic positive and adversarial tests.
4. Re-run failed extractions and all decoys. A negative benchmark label cannot be restored by manufacturing successful requests, suppressing an extraction error, or interpreting an unparsed assertion as non-support.
5. Compare every previously correct case before acceptance. No scalability refactor or claim of stable 12/12 is justified yet.

## Reproduction

```bash
python -B -m unittest discover -s tests -p 'test*correctness.py'
python -B replay_correctness.py
```

On Nibi, from the project directory with its existing agent_env and model, run:

```bash
sbatch run_correctness_benchmark.sh
```

The runner executes both public benchmark paths and records results/correctness_<job>/pairs.csv, pairs.csv.diagnostics.jsonl, env.csv, env.csv.diagnostics.jsonl and run.log. It uses separate correctness cache files and never changes expected labels. See CORRECTNESS_IMPLEMENTATION_REPORT.md for implementation details and limitations.

## Provenance

Fixture SHA-256: `06f496862dcf8cc13d435e721ac0184a2db77e1b3c831a33dd9e3fd293c1f5ec`.

Unchanged test_pairs.csv SHA-256: `be0f564702c1d8ba2654e7f2aa920cac8532060df79a61dac7a89ae8a550106b`.

Full replay decisions, original copied extraction fields and code hashes: results/benchmark_after_correctness_offline.json. The source fixture results/benchmark_diagnostic_v06.json remains unchanged.
