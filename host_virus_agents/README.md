# Host-virus association research pipeline

**Status: correctness iteration awaiting Nibi validation. Not production-ready.**

This repository tests literature evidence for predicted host-virus associations. Directed host-to-virus evidence, related-virus comparisons, biological priors and search failures are kept separate. The latest user-reported live benchmark is 0/12; the updated implementation has not yet been run live. See `CORRECTNESS_IMPLEMENTATION_REPORT.md` and `FINAL_VALIDATION_REPORT.md`.

## Local correctness tests

From the repository root with Python 3.10 or newer:

```bash
python -B -m unittest discover -v
```

The 77 current tests require only the standard library and use synthetic fixtures/mocked external boundaries. They do not download a model or query literature. Python 3.12.14 was used locally.

The historical artifact test additionally runs in PowerShell:

```powershell
./test_benchmark_diagnostics.ps1
```

It validates the original 12 rows and 60 paper records. It does not prove biological correctness.

## Run the updated live benchmark on Nibi

Use the existing project `agent_env` and staged `models/Qwen3-8B` used by the previous live run. The checked-in `requirements.txt` describes the existing Alliance environment and is not advertised as a portable Windows installation lock. This iteration introduces no new runtime dependency for the decision logic or diagnostic tooling.

Submit from the project root:

```bash
sbatch --account=def-acdoxey --gres=gpu:nvidia_h100_80gb_hbm3_3g.40gb:1 run_correctness_benchmark.sh
```

The runner executes offline tests, then both benchmark entry points with the unchanged `test_pairs.csv`. It creates `results/correctness_<job-id>_<UTC timestamp>/`, containing CSVs, full diagnostic JSONL, logs, validation summaries and separate cache files. `MODEL_PATH` can override the staged model location. Hashing weights for the reproducibility manifest adds startup I/O.

Return the entire run directory for review. The validator checks label arithmetic and diagnostic presence; it does not certify biological attribution. All 12 supporting/related/rejected decisions must be inspected before freezing correctness.

For a single exploratory pair in that model environment:

```bash
source agent_env/bin/activate
python judge_agent.py "HOST NAME" "VIRUS NAME"
```

This is the legacy research path, not a scalable batch application. Its `NO_EVIDENCE_FOUND` label is not novelty. Canonical uncertainty is `INSUFFICIENT_EVIDENCE`.

## Production and scaling gate

The intended production labels are `KNOWN`, `NOVEL_CANDIDATE`, `BIOLOGICALLY_IMPLAUSIBLE` and `INSUFFICIENT_EVIDENCE`. Novelty must require adequate dated coverage and defensible plausibility; missing papers, failed extraction and unresolved taxonomy cannot establish it.

`classify_pairs.py`, persistent evidence indexing, Parquet, restartable batch jobs and scale tests are intentionally not implemented before the live correctness gate. Do not submit millions of rows to the current pair-by-pair research workflow. The proposed migration is in `ARCHITECTURE_REVIEW.md`; implementation remains gated as recorded in `V2_IMPLEMENTATION_REPORT.md`.
