#!/bin/bash
#SBATCH --job-name=hostvirus_agents
#SBATCH --account=def-acdoxey
#SBATCH --time=00:30:00
#SBATCH --gres=gpu:nvidia_h100_80gb_hbm3_3g.40gb:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48000M
#SBATCH --mail-user=s2mascar@uwaterloo.ca
#SBATCH --mail-type=ALL

cd /project/6002943/smascar/HostFinder/host_virus_agents

source agent_env/bin/activate

mkdir -p logs

echo "============================================================"
echo "NODE"
echo "============================================================"

hostname

echo
echo "============================================================"
echo "GPU"
echo "============================================================"

nvidia-smi

echo
echo "============================================================"
echo "STARTING 12-PAIR BENCHMARK"
echo "============================================================"

python evaluate_pairs.py \
    test_pairs.csv \
    benchmark_results.csv


echo
echo "============================================================"
echo "BENCHMARK FINISHED"
echo "STARTING EVIDENCE DIAGNOSTICS"
echo "============================================================"

DIAG_FILE="logs/evidence_diagnostics_${SLURM_JOB_ID}.txt"

python - <<'PY' | tee "$DIAG_FILE"

from evidence_agent import run_evidence_agent


# ============================================================
# CASES TO INSPECT IN DETAIL
# ============================================================

cases = [

    # --------------------------------------------------------
    # Positive that keeps failing
    # --------------------------------------------------------
    (
        "Lampyris noctiluca",
        "Lampyris noctiluca partiti-like virus 1",
        "KNOWN"
    ),

    # --------------------------------------------------------
    # Positive control that has sometimes changed behaviour
    # --------------------------------------------------------
    (
        "Homo sapiens",
        "SARS-CoV-2",
        "KNOWN"
    ),

    # --------------------------------------------------------
    # Related-literature cases
    # --------------------------------------------------------
    (
        "Lampyris noctiluca",
        "Hubei partiti-like virus 31",
        "POSSIBLY_KNOWN"
    ),

    (
        "Lampyris noctiluca",
        "Hubei partiti-like virus 51",
        "POSSIBLY_KNOWN"
    ),

    (
        "Lampyris noctiluca",
        "Hubei chuvirus-like virus 3",
        "POSSIBLY_KNOWN"
    ),

    (
        "Lampyris noctiluca",
        "Hubei toti-like virus 16",
        "POSSIBLY_KNOWN"
    ),

    # --------------------------------------------------------
    # Important negative specificity control
    # --------------------------------------------------------
    (
        "Tribolium castaneum",
        "Hubei partiti-like virus 31",
        "NO_EVIDENCE_FOUND"
    )
]


def short(text, limit=1800):

    if not text:
        return ""

    text = str(text)

    if len(text) <= limit:
        return text

    return text[:limit] + "\n...[TRUNCATED]..."


for case_number, (host, virus, expected) in enumerate(
    cases,
    start=1
):

    print("\n")
    print("=" * 100)
    print(f"DIAGNOSTIC CASE {case_number}/{len(cases)}")
    print("=" * 100)

    print(f"HOST:     {host}")
    print(f"VIRUS:    {virus}")
    print(f"EXPECTED: {expected}")

    print("=" * 100)

    try:

        output = run_evidence_agent(
            host,
            virus
        )

    except Exception as error:

        print(
            f"\nERROR RUNNING CASE:\n{error}"
        )

        continue

    print("\nHOST ALIASES:")

    for alias in output.get(
        "host_aliases",
        []
    ):
        print(f"  - {alias}")

    search_metadata = output.get(
        "search_metadata",
        {}
    )

    print("\nSEARCH METADATA:")

    print(
        "Queries executed:",
        search_metadata.get(
            "queries_executed"
        )
    )

    print(
        "Candidate papers:",
        search_metadata.get(
            "candidate_papers"
        )
    )

    print(
        "Retrieval complete:",
        search_metadata.get(
            "retrieval_complete"
        )
    )

    print(
        "Source failures:",
        len(
            search_metadata.get(
                "source_failures",
                []
            )
        )
    )

    results = output.get(
        "evidence_results",
        []
    )

    print(
        "\nPAPERS ANALYZED:",
        len(results)
    )

    # --------------------------------------------------------
    # Summary of evidence classes
    # --------------------------------------------------------

    counts = {}

    for result in results:

        classification = result.get(
            "classification",
            "MISSING"
        )

        counts[classification] = (
            counts.get(
                classification,
                0
            )
            + 1
        )

    print("\nEVIDENCE CLASS COUNTS:")

    for classification, count in counts.items():

        print(
            f"  {classification}: {count}"
        )

    # --------------------------------------------------------
    # Paper-by-paper diagnostics
    # --------------------------------------------------------

    for paper_number, result in enumerate(
        results,
        start=1
    ):

        print("\n")
        print("-" * 100)
        print(
            f"PAPER {paper_number}/{len(results)}"
        )
        print("-" * 100)

        print(
            "TITLE:",
            result.get(
                "title"
            )
        )

        print(
            "PMID:",
            result.get(
                "pmid"
            )
        )

        print(
            "PMCID:",
            result.get(
                "pmcid"
            )
        )

        print()
        print(
            "CLASSIFICATION:",
            result.get(
                "classification"
            )
        )

        print(
            "RELATIONSHIP TYPE:",
            result.get(
                "relationship_type"
            )
        )

        print()
        print(
            "STUDY HOST:",
            result.get(
                "study_host"
            )
        )

        print(
            "STUDY HOST MATCHES TARGET:",
            result.get(
                "study_host_matches_target"
            )
        )

        print(
            "STRONG STUDY HOST CONTEXT:",
            result.get(
                "strong_study_host_context"
            )
        )

        print(
            "STUDY HOST PASSAGE VERIFIED:",
            result.get(
                "study_host_passage_verified"
            )
        )

        print()
        print(
            "REPORTED VIRUS:",
            result.get(
                "reported_virus_name"
            )
        )

        print(
            "REPORTED VIRUS MATCHES TARGET:",
            result.get(
                "reported_virus_matches_target"
            )
        )

        print(
            "REPORTED VIRUS PASSAGE VERIFIED:",
            result.get(
                "reported_virus_passage_verified"
            )
        )

        print()
        print(
            "REPORTED HOST/SAMPLE:",
            result.get(
                "reported_host_or_sample"
            )
        )

        print(
            "REPORTED HOST MATCHES TARGET:",
            result.get(
                "reported_host_matches_target"
            )
        )

        print()
        print(
            "RELATIONSHIP PASSAGE VERIFIED:",
            result.get(
                "relationship_passage_verified"
            )
        )

        print(
            "RELATIONSHIP MENTIONS TARGET VIRUS:",
            result.get(
                "relationship_mentions_target_virus"
            )
        )

        print(
            "DIRECT RELATION EXPLICIT:",
            result.get(
                "direct_relation_explicit"
            )
        )

        print()
        print("MODEL REASON:")

        print(
            result.get(
                "reason"
            )
        )

        print()
        print("STUDY-HOST PASSAGE:")

        print(
            short(
                result.get(
                    "study_host_passage"
                )
            )
        )

        print()
        print("REPORTED-VIRUS PASSAGE:")

        print(
            short(
                result.get(
                    "reported_virus_passage"
                )
            )
        )

        print()
        print("RELATIONSHIP PASSAGE:")

        print(
            short(
                result.get(
                    "relationship_passage"
                )
            )
        )


print("\n")
print("=" * 100)
print("ALL DIAGNOSTIC CASES COMPLETE")
print("=" * 100)

PY


echo
echo "============================================================"
echo "JOB COMPLETE"
echo "============================================================"

echo "Benchmark:"
echo "  benchmark_results.csv"

echo
echo "Detailed evidence diagnostics:"
echo "  ${DIAG_FILE}"