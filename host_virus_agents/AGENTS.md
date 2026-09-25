# Host-Virus Association Validation Pipeline

## Project purpose

This repository contains a research pipeline for determining whether predicted host-virus associations are already known in the scientific literature, potentially novel, or biologically implausible.

The predictions originate from HostFinder, which produces large numbers of candidate host-virus associations from large-scale SRA co-occurrence analysis.

The long-term goal is to process millions of host-virus predictions efficiently.

## Required classifications

The pipeline should ultimately support:

1. `KNOWN`
   Credible literature evidence supports the exact host-virus association.

2. `NOVEL_CANDIDATE`
   No credible prior evidence was found after adequate search coverage, and the association remains biologically plausible.

3. `BIOLOGICALLY_IMPLAUSIBLE`
   Taxonomic or biological evidence strongly conflicts with the proposed association.

4. `INSUFFICIENT_EVIDENCE`
   Search coverage, taxonomy resolution, or evidence quality is inadequate to confidently classify the pair.

A failed literature query or zero search results must never automatically be interpreted as novelty.

## Existing code

Important modules include:

- `search_agent.py`
- `literature_search.py`
- `pubmed_search.py`
- `evidence_agent.py`
- `judge_agent.py`
- `taxonomy_aliases.py`
- `evaluate_pairs.py`
- `evaluate_env.py`

Previous implementations and experiments are also present in:

- `bioresearch_env/`
- `BioResearchEnv_v0.5/`
- `backups/`

Benchmark and test files include:

- `test_pairs.csv`
- `benchmark_results.csv`
- `run_benchmark.sh`
- `run_benchmark_v2.sh`

Inspect previous versions before redesigning the project because later changes may contain regressions relative to earlier versions.

## Biological semantics

Simple co-occurrence of a host name and virus name in the same paper is NOT sufficient evidence of an association.

The system must distinguish between:

- target host
- target virus
- other hosts
- comparison viruses
- natural hosts
- experimental hosts
- environmental source organisms
- sequence similarity statements
- background mentions
- discovery/detection evidence
- actual host-virus association evidence

Evidence for:

`OTHER_HOST -> TARGET_VIRUS`

must not be interpreted as evidence for:

`TARGET_HOST -> TARGET_VIRUS`

Likewise, evidence involving a related or comparison virus must not be attributed to the target virus.

## Scalability requirement

The production architecture must eventually support millions of predicted host-virus pairs.

Therefore:

- Do not design around one PubMed query per pair.
- Do not design around one LLM invocation per pair.
- Reuse work across repeated hosts and viruses.
- Cache taxonomy records.
- Cache aliases.
- Cache literature queries.
- Cache papers.
- Cache extracted evidence.
- Deduplicate external requests.
- Prefer joins and indexed lookups over repeated inference.
- Support batch processing.
- Support checkpointing and resume.
- Keep expensive reasoning restricted to ambiguous candidate evidence.

The intended architecture should resemble:

input pairs
→ normalization
→ taxonomy resolution
→ biological plausibility
→ reusable literature retrieval
→ candidate paper filtering
→ evidence extraction
→ evidence storage
→ classification

Literature-derived associations should become a reusable structured knowledge base rather than being rediscovered independently for every input pair.

## Input contract

At minimum support:

`host,virus`

Preferably support:

`host_taxid,host,virus_taxid,virus`

CSV support is required.

Parquet support should be added for large datasets.

## Output contract

Each input association should eventually produce:

- host
- virus
- host_taxid
- virus_taxid
- classification
- classification_reason
- biological_plausibility
- plausibility_reason
- literature_search_status
- papers_retrieved
- evidence_count
- best_supporting_paper
- confidence

Evidence provenance must be retained.

## Engineering requirements

The system should:

- remain modular
- have unit tests
- have integration tests
- preserve useful existing behavior
- fail gracefully when APIs fail
- support checkpoints
- avoid duplicate external requests
- use structured logging
- separate API/network code from biological decision logic
- use structured data models where appropriate
- avoid hard-coded benchmark cases
- remain runnable on Linux HPC systems such as Digital Research Alliance Canada clusters

Do not change biological semantics merely to improve benchmark accuracy.

## Development procedure

Before large changes:

1. Inspect the complete current implementation.
2. Inspect previous implementations.
3. Inspect benchmarks and previous outputs.
4. Understand the current execution flow.
5. Identify correctness failures.
6. Identify scalability failures.
7. Propose a replacement architecture.
8. Preserve or add tests.
9. Refactor incrementally.
10. Run benchmarks after major milestones.

Do not delete previous useful behavior until its replacement has been validated.
