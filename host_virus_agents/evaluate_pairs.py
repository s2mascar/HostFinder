import csv
import sys
import time

from judge_agent import run_judge_agent
import benchmark_capture


def run_benchmark(input_file, output_file):
    sidecar = benchmark_capture.begin(input_file, output_file)

    with open(
        input_file,
        "r",
        encoding="utf-8"
    ) as f:

        pairs = list(
            csv.DictReader(f)
        )

    print("=" * 80)
    print("HOST-VIRUS LITERATURE BENCHMARK")
    print("=" * 80)

    print(f"Pairs to test: {len(pairs)}")

    results = []

    correct = 0

    benchmark_start = time.time()

    for i, row in enumerate(
        pairs,
        start=1
    ):

        host = row["host"]
        virus = row["virus"]
        expected = row["expected_status"]

        print("\n" + "#" * 80)

        print(
            f"PAIR {i}/{len(pairs)}"
        )

        print("#" * 80)

        print(f"Host:     {host}")
        print(f"Virus:    {virus}")
        print(f"Expected: {expected}")

        start = time.time()
        result = {}

        try:

            result = run_judge_agent(
                host,
                virus
            )

            predicted = result.get(
                "literature_status",
                "ERROR"
            )

            runtime = time.time() - start

            is_correct = (
                predicted == expected
            )

            if is_correct:
                correct += 1

            print(
                f"\nExpected:  {expected}"
            )

            print(
                f"Predicted: {predicted}"
            )

            print(
                f"Correct:   {is_correct}"
            )

            print(
                f"Runtime:   {runtime:.1f} sec"
            )

            results.append({
                "host": host,
                "virus": virus,

                "expected_status":
                    expected,

                "predicted_status":
                    predicted,

                "correct":
                    is_correct,

                "papers_examined":
                    result.get(
                        "papers_examined",
                        0
                    ),

                "exact_supporting_papers":
                    len(
                        result.get(
                            "exact_supporting_papers",
                            []
                        )
                    ),

                "related_papers":
                    len(
                        result.get(
                            "related_papers",
                            []
                        )
                    ),

                "confidence":
                    result.get(
                        "confidence",
                        ""
                    ),

                "reason":
                    result.get(
                        "reason",
                        ""
                    ),

                "runtime_seconds":
                    round(runtime, 2),

                "test_type":
                    row.get(
                        "test_type",
                        ""
                    ),

                "reference_hint":
                    row.get(
                        "reference_hint",
                        ""
                    ),
            })

        except Exception as error:
            result = {"classification": "INSUFFICIENT_EVIDENCE", "processing_status": "FAILED", "error": str(error)}

            runtime = time.time() - start

            print(
                f"\nERROR: {error}"
            )

            results.append({
                "host": host,
                "virus": virus,
                "expected_status": expected,
                "predicted_status": "UNCLEAR",
                "correct": False,
                "papers_examined": 0,
                "exact_supporting_papers": 0,
                "related_papers": 0,
                "confidence": "",
                "reason": str(error),
                "runtime_seconds":
                    round(runtime, 2),
                "test_type":
                    row.get(
                        "test_type",
                        ""
                    ),
                "reference_hint":
                    row.get(
                        "reference_hint",
                        ""
                    ),
            })

        results[-1]["classification"] = result.get("classification", "INSUFFICIENT_EVIDENCE")
        results[-1]["processing_status"] = result.get("processing_status", "COMPLETED")
        benchmark_capture.append(sidecar, results[-1], result)
        # Save after every pair so progress is not lost
        with open(
            output_file,
            "w",
            newline="",
            encoding="utf-8"
        ) as f:

            writer = csv.DictWriter(
                f,
                fieldnames=results[0].keys()
            )

            writer.writeheader()

            writer.writerows(
                results
            )

    total_time = (
        time.time()
        - benchmark_start
    )

    accuracy = (
        correct / len(pairs)
        if pairs
        else 0
    )

    print("\n" + "=" * 80)
    print("BENCHMARK COMPLETE")
    print("=" * 80)

    print(
        f"Correct: {correct}/{len(pairs)}"
    )

    print(
        f"Accuracy: {accuracy:.1%}"
    )

    print(
        f"Total runtime: "
        f"{total_time / 60:.1f} minutes"
    )

    print(
        f"Results written to: "
        f"{output_file}"
    )


if __name__ == "__main__":

    if len(sys.argv) != 3:

        print(
            "Usage:\n"
            "python evaluate_pairs.py "
            "test_pairs.csv "
            "benchmark_results.csv"
        )

        sys.exit(1)

    run_benchmark(
        sys.argv[1],
        sys.argv[2]
    )
