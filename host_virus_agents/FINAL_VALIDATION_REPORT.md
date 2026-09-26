# Validation status — project NOT complete

Date: 2026-09-26. This is a status report, not a production acceptance certificate.

| Definition of Done item | Status | Evidence / blocker |
|---|---|---|
| 1. All unit tests pass | PASS | Current local unittest suite passes; original 42 retained. |
| 2. All integration tests pass | PASS | Current offline/mocked integration tests pass; live model validation is separately gated below. Total local Python tests: 77. |
| 3. Live Nibi benchmark 12/12 | FAIL | Latest authoritative user summary is 0/12. Post-change run requires Nibi and is pending. |
| 4. All 12 evidence attributions manually auditable | BLOCKED | New full diagnostic run and manual review required. The old live raw artifacts are not local. |
| 5. Zero known benchmark-specific hacks | PASS | Active source scan finds no benchmark-name/PMID decision exceptions. Inactive legacy human-normalization references are documented. Labels unchanged. |
| 6. Failed retrieval/extraction does not become novelty | PASS | Failure/coverage/materiality tests pass; no novelty promotion exists in this correctness iteration. This does not validate a future novelty implementation. |
| 7. classify_pairs.py exists | BLOCKED | Production batch implementation intentionally gated on correctness. |
| 8. Production CSV input works | BLOCKED | Legacy benchmark CSV capture is tested; production batch CLI is not implemented. |
| 9. Parquet input works | BLOCKED | Not implemented before correctness freeze. |
| 10. Persistent production caching works | BLOCKED | Existing JSON caches are not a validated production evidence store. |
| 11. Duplicate expensive work reused | BLOCKED | Duplicate-row capture/repeated-evidence behavior is tested; cross-pair persistent reuse is not implemented. |
| 12. Checkpoint/resume works | BLOCKED | Per-row diagnostic persistence exists; production resume is not implemented. |
| 13. Large synthetic batch passes | BLOCKED | Scale tests gated; no throughput result claimed. |
| 14. Linux/Nibi execution works | BLOCKED | Previous code ran on Nibi per user. Updated runner has not been executed there; Bash/Slurm unavailable locally. |
| 15. Exact installation/execution documentation | BLOCKED | README documents local tests and existing Nibi-environment execution; a production installation/deployment guide awaits V2. |
| 16. Old experimental code separated from production | BLOCKED | Inactive legacy verifier retained; smoke model demo is import-safe. Full production separation is gated. |
| 17. Output filterable for NOVEL_CANDIDATE | BLOCKED | Production four-class policy/CLI not implemented. No apparent novelty is fabricated. |
| 18. Fresh non-benchmark end-to-end smoke test | BLOCKED | Synthetic non-benchmark unit/integration inputs pass; no fresh biological production run exists. |

The current milestone is a tested correctness iteration ready for an external Nibi run. It is not stable 12/12, production-ready, or ready for millions of predictions. Execute the command in `README.md`, return the complete run directory, then continue repairs using actual evidence traces. Do not freeze correctness or begin scaling until its live acceptance gates pass.
