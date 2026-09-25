# Host-Virus Pipeline V2 Specification

## Objective

Convert the existing host-virus literature research prototype into a robust pipeline capable of classifying large numbers of predicted host-virus associations.

The intended downstream use is to prioritize novel HostFinder predictions for further genomic and experimental validation.

## Primary question

For each proposed host-virus association:

`HOST -> VIRUS`

determine whether the best available evidence indicates:

- `KNOWN`
- `NOVEL_CANDIDATE`
- `BIOLOGICALLY_IMPLAUSIBLE`
- `INSUFFICIENT_EVIDENCE`

## Critical interpretation rule

Absence of retrieved literature is not equivalent to biological novelty.

A `NOVEL_CANDIDATE` result requires both:

1. adequate literature-search coverage with no credible evidence supporting the exact association;
2. no strong biological evidence making the association implausible.

## Current prototype

The existing system uses literature search, taxonomy normalization, evidence extraction, biological context, and an evidence judge.

Previous versions introduced improvements such as:

- multiple literature query strategies
- taxonomy-aware aliases
- host and virus normalization
- distinction between target and comparison viruses
- evidence extraction from papers
- biological context caching
- explicit host-virus relationship types
- benchmark-driven development

Some later versions also introduced regressions.

The existing project should therefore be audited rather than replaced blindly.

## Target architecture

### Stage 1: Input normalization

Read host-virus associations in batches.

Canonicalize:

- host name
- virus name
- optional host TaxID
- optional virus TaxID

Deduplicate internal processing while preserving one output row per input record.

### Stage 2: Taxonomy resolution

Resolve and cache:

- canonical scientific names
- TaxIDs
- synonyms
- taxonomic lineage
- relevant ranks
- virus family/genus/species where available

Do this once per unique entity, not once per pair.

### Stage 3: Biological plausibility

Use inspectable rules and structured biological information to determine whether a proposed interaction is:

- plausible
- uncertain
- strongly implausible

Do not rely solely on free-form LLM reasoning.

### Stage 4: Literature retrieval

Retrieve literature using reusable searches based on unique hosts, viruses, taxonomic aliases, and entities.

Avoid executing independent literature searches for every host-virus pair.

Cache:

- queries
- results
- papers
- retrieval status

### Stage 5: Candidate-paper matching

Determine which retrieved papers may contain evidence relevant to a given pair.

Use inexpensive filtering before expensive evidence extraction.

### Stage 6: Evidence extraction

Extract explicit structured biological relationships.

The extraction system must distinguish:

- target virus detected in target host
- target virus detected in another host
- another virus detected in target host
- experimental infection
- natural infection
- environmental detection
- sequence similarity only
- taxonomic comparison
- background mention
- discovery evidence

Positive evidence must retain provenance to the source paper and supporting text.

### Stage 7: Evidence store

Store extracted host-virus relationships independently of prediction jobs.

Suggested logical schema:

`host_taxid`
`virus_taxid`
`host_name`
`virus_name`
`relationship_type`
`evidence_strength`
`paper_id`
`supporting_text`
`source`
`extraction_version`

This evidence should be reusable for future HostFinder predictions.

### Stage 8: Classification

Use taxonomy, plausibility, literature-search status, and structured evidence to produce the final classification.

### Stage 9: Batch processing

Provide a command approximately like:

`python classify_pairs.py --input predictions.csv --output classifications.parquet`

Support:

- CSV
- Parquet
- batching
- resumability
- persistent caches
- partial output recovery

## Scaling philosophy

A dataset may contain tens of millions of predicted associations but far fewer unique:

- hosts
- viruses
- papers
- taxonomy entities

The number of expensive operations should therefore scale primarily with unique entities and candidate evidence, not directly with the number of pair rows.

## Benchmarking

Preserve the current benchmark and improve it.

Track:

- total accuracy
- false KNOWN calls
- false NOVEL calls
- false BIOLOGICALLY_IMPLAUSIBLE calls
- failed searches
- unresolved taxonomy
- papers retrieved
- evidence extraction failures
- runtime
- cache hits
- cache misses

Benchmark expected labels must not be altered simply to improve performance.

## Near-term success criterion

Before testing millions of associations, V2 should be able to:

1. accept a list of host-virus pairs;
2. resolve taxonomy consistently;
3. reuse cached retrieval work;
4. retrieve candidate literature;
5. distinguish exact association evidence from misleading co-mentions;
6. classify each pair;
7. produce a transparent machine-readable explanation;
8. reproduce or improve the strongest previous benchmark behavior without introducing obvious biological errors.
