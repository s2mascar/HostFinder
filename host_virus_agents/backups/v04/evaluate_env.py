import csv
import json
import statistics
import sys
import time

from bioresearch_env.env import (
    BioResearchEnv,
)

from bioresearch_env.baseline_agent import (
    BaselineResearchAgent,
)


# ============================================================
# CONTEXT OUTPUT HELPERS
# ============================================================


def summarize_context(context):
    context = context or {}

    target_host = context.get(
        "target_host",
        {},
    )

    target_virus = context.get(
        "target_virus",
        {},
    )

    host_range = context.get(
        "host_range",
        {},
    )

    biological_prior = context.get(
        "biological_prior",
        {},
    )

    closest_hosts = host_range.get(
        "closest_known_hosts",
        [],
    )

    closest = (
        closest_hosts[0]
        if closest_hosts
        else {}
    )

    return {
        "biological_prior": biological_prior.get(
            "status",
            "UNKNOWN",
        ),
        "virus_host_match_level": host_range.get(
            "match_level",
            "",
        ),
        "target_host_group": target_host.get(
            "broad_group",
            "UNKNOWN",
        ),
        "target_host_family": target_host.get(
            "family",
            "",
        ),
        "target_host_order": target_host.get(
            "order",
            "",
        ),
        "target_virus_family": target_virus.get(
            "family",
            "",
        ),
        "target_virus_genus": target_virus.get(
            "genus",
            "",
        ),
        "known_host_groups": json.dumps(
            host_range.get(
                "broad_group_counts",
                {},
            ),
            sort_keys=True,
        ),
        "known_non_target_hosts": host_range.get(
            "known_non_target_hosts",
            0,
        ),
        "closest_known_host": closest.get(
            "host_name",
            "",
        ),
        "closest_host_lca_rank": closest.get(
            "lca_rank",
            "",
        ),
        "closest_host_lca_name": closest.get(
            "lca_name",
            "",
        ),
        "environmental_record_count": host_range.get(
            "environmental_record_count",
            0,
        ),
        "environmental_contexts": json.dumps(
            host_range.get(
                "environmental_contexts",
                [],
            ),
            sort_keys=True,
        ),
    }


# ============================================================
# EVIDENCE EDGE DIAGNOSTICS
# ============================================================


def summarize_evidence(analyzed_results):
    analyzed_results = analyzed_results or {}

    classification_counts = {}
    extraction_status_counts = {}
    paper_diagnostics = []
    verified_host_edges = 0
    verified_comparison_edges = 0
    rescued_extractions = 0
    failed_extractions = 0
    virus_taxonomy_resolved_papers = 0

    for paper_id, result in analyzed_results.items():
        classification = result.get(
            "classification",
            "",
        )

        classification_counts[classification] = (
            classification_counts.get(
                classification,
                0,
            )
            + 1
        )

        extraction_status = result.get(
            "extraction_status",
            "UNKNOWN",
        )

        extraction_status_counts[extraction_status] = (
            extraction_status_counts.get(
                extraction_status,
                0,
            )
            + 1
        )

        if extraction_status == "RESCUED":
            rescued_extractions += 1

        if extraction_status == "FAILED":
            failed_extractions += 1

        if result.get("virus_taxonomy_resolved"):
            virus_taxonomy_resolved_papers += 1

        host_edge = result.get(
            "host_edge",
            {},
        ) or {}

        comparison_edge = result.get(
            "comparison_edge",
            {},
        ) or {}

        if host_edge.get("verified"):
            verified_host_edges += 1

        if comparison_edge.get("verified"):
            verified_comparison_edges += 1

        paper_diagnostics.append({
            "paper_id": paper_id,
            "title": result.get("title"),
            "pmid": result.get("pmid"),
            "pmcid": result.get("pmcid"),
            "classification": classification,
            "classification_basis": result.get(
                "classification_basis"
            ),
            "extraction_status": result.get(
                "extraction_status"
            ),
            "extraction_attempts": result.get(
                "extraction_attempts"
            ),
            "extraction_error": result.get(
                "extraction_error"
            ),
            "target_virus_scientific_name": result.get(
                "target_virus_scientific_name"
            ),
            "target_virus_tax_id": result.get(
                "target_virus_tax_id"
            ),
            "target_virus_aliases": result.get(
                "target_virus_aliases",
                [],
            ),
            "study_host": result.get(
                "study_host"
            ),
            "host_virus_name": result.get(
                "host_virus_name"
            ),
            "comparison_source_virus_name": result.get(
                "comparison_source_virus_name"
            ),
            "comparison_virus_name": result.get(
                "comparison_virus_name"
            ),
            "host_edge": host_edge,
            "comparison_edge": comparison_edge,
        })

    return {
        "classification_counts": json.dumps(
            classification_counts,
            sort_keys=True,
        ),
        "extraction_status_counts": json.dumps(
            extraction_status_counts,
            sort_keys=True,
        ),
        "rescued_extractions": rescued_extractions,
        "failed_extractions": failed_extractions,
        "virus_taxonomy_resolved_papers": virus_taxonomy_resolved_papers,
        "verified_host_edges": verified_host_edges,
        "verified_comparison_edges": verified_comparison_edges,
        "paper_diagnostics": json.dumps(
            paper_diagnostics,
            sort_keys=True,
            ensure_ascii=False,
        ),
    }


# ============================================================
# EVALUATION
# ============================================================


def run_environment_benchmark(
    input_file,
    output_file,
):
    with open(
        input_file,
        "r",
        encoding="utf-8",
    ) as handle:
        tasks = list(
            csv.DictReader(handle)
        )

    print("=" * 80)
    print("BIORESEARCHENV v0.5 BENCHMARK")
    print("=" * 80)
    print(f"Episodes: {len(tasks)}")

    results = []
    benchmark_start = time.time()

    agent = BaselineResearchAgent()

    for episode_number, task in enumerate(
        tasks,
        start=1,
    ):
        print("\n" + "#" * 80)
        print(
            f"EPISODE {episode_number}/{len(tasks)}"
        )
        print("#" * 80)
        print(f"Host:  {task['host']}")
        print(f"Virus: {task['virus']}")

        episode_start = time.time()

        try:
            env = BioResearchEnv(
                task=task,
                max_steps=15,
                max_papers=8,
            )

            info = agent.run(env)

            runtime = (
                time.time()
                - episode_start
            )

            context_fields = summarize_context(
                env.get_biological_context()
            )

            evidence_fields = summarize_evidence(
                env.get_analyzed_results()
            )

            row = {
                "host": task["host"],
                "virus": task["virus"],
                "expected_status": info.get(
                    "expected_status",
                    task["expected_status"],
                ),
                "predicted_status": info.get(
                    "predicted_status"
                ),
                "correct": info.get(
                    "correct",
                    False,
                ),
                "episode_reward": round(
                    info.get(
                        "episode_reward",
                        0.0,
                    ),
                    4,
                ),
                "steps": info.get(
                    "steps",
                    0,
                ),
                "searches": info.get(
                    "searches",
                    0,
                ),
                "papers_retrieved": info.get(
                    "papers_retrieved",
                    0,
                ),
                "papers_analyzed": info.get(
                    "papers_analyzed",
                    0,
                ),
                "exact_supporting_papers": info.get(
                    "exact_supporting_papers",
                    0,
                ),
                "related_papers": info.get(
                    "related_papers",
                    0,
                ),
                "confidence": info.get(
                    "confidence",
                    "",
                ),
                "reason": info.get(
                    "reason",
                    "",
                ),
                **context_fields,
                **evidence_fields,
                "runtime_seconds": round(
                    runtime,
                    2,
                ),
                "test_type": task.get(
                    "test_type",
                    "",
                ),
                "reference_hint": task.get(
                    "reference_hint",
                    "",
                ),
            }

        except Exception as error:
            runtime = (
                time.time()
                - episode_start
            )

            print(
                f"\nEPISODE ERROR: {error}"
            )

            row = {
                "host": task["host"],
                "virus": task["virus"],
                "expected_status": task["expected_status"],
                "predicted_status": "ERROR",
                "correct": False,
                "episode_reward": 0.0,
                "steps": 0,
                "searches": 0,
                "papers_retrieved": 0,
                "papers_analyzed": 0,
                "exact_supporting_papers": 0,
                "related_papers": 0,
                "confidence": "",
                "reason": str(error),
                "biological_prior": "ERROR",
                "virus_host_match_level": "",
                "target_host_group": "",
                "target_host_family": "",
                "target_host_order": "",
                "target_virus_family": "",
                "target_virus_genus": "",
                "known_host_groups": "{}",
                "known_non_target_hosts": 0,
                "closest_known_host": "",
                "closest_host_lca_rank": "",
                "closest_host_lca_name": "",
                "environmental_record_count": 0,
                "environmental_contexts": "[]",
                "classification_counts": "{}",
                "extraction_status_counts": "{}",
                "rescued_extractions": 0,
                "failed_extractions": 0,
                "virus_taxonomy_resolved_papers": 0,
                "verified_host_edges": 0,
                "verified_comparison_edges": 0,
                "paper_diagnostics": "[]",
                "runtime_seconds": round(
                    runtime,
                    2,
                ),
                "test_type": task.get(
                    "test_type",
                    "",
                ),
                "reference_hint": task.get(
                    "reference_hint",
                    "",
                ),
            }

        results.append(row)

        print()
        print(
            "Predicted:",
            row["predicted_status"],
        )
        print(
            "Expected:",
            row["expected_status"],
        )
        print(
            "Correct:",
            row["correct"],
        )
        print(
            "Reward:",
            row["episode_reward"],
        )
        print(
            "Biological prior:",
            row["biological_prior"],
        )
        print(
            "Papers analyzed:",
            row["papers_analyzed"],
        )

        # Save after every episode so interrupted SLURM jobs retain
        # completed results.
        with open(
            output_file,
            "w",
            newline="",
            encoding="utf-8",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=results[0].keys(),
            )
            writer.writeheader()
            writer.writerows(results)

    # ========================================================
    # SUMMARY
    # ========================================================

    total_runtime = (
        time.time()
        - benchmark_start
    )

    completed = [
        row
        for row in results
        if row["predicted_status"] != "ERROR"
    ]

    correct_count = sum(
        1
        for row in completed
        if row["correct"]
    )

    accuracy = (
        correct_count / len(completed)
        if completed
        else 0.0
    )

    mean_reward = (
        statistics.mean(
            row["episode_reward"]
            for row in completed
        )
        if completed
        else 0.0
    )

    mean_steps = (
        statistics.mean(
            row["steps"]
            for row in completed
        )
        if completed
        else 0.0
    )

    mean_searches = (
        statistics.mean(
            row["searches"]
            for row in completed
        )
        if completed
        else 0.0
    )

    mean_papers = (
        statistics.mean(
            row["papers_analyzed"]
            for row in completed
        )
        if completed
        else 0.0
    )

    print("\n" + "=" * 80)
    print("BIORESEARCHENV BENCHMARK COMPLETE")
    print("=" * 80)
    print(
        f"Completed: {len(completed)}/{len(tasks)}"
    )
    print(
        f"Correct: {correct_count}/{len(completed)}"
    )
    print(
        f"Accuracy: {accuracy:.1%}"
    )
    print(
        f"Mean episode reward: {mean_reward:.3f}"
    )
    print(
        f"Mean steps: {mean_steps:.2f}"
    )
    print(
        f"Mean searches: {mean_searches:.2f}"
    )
    print(
        f"Mean papers analyzed: {mean_papers:.2f}"
    )
    print(
        f"Total runtime: {total_runtime / 60:.1f} minutes"
    )
    print(
        f"Results written to: {output_file}"
    )


# ============================================================
# CLI
# ============================================================


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(
            "Usage:\n"
            "python evaluate_env.py "
            "test_pairs.csv "
            "results/bioresearch_env_results.csv"
        )
        sys.exit(1)

    run_environment_benchmark(
        sys.argv[1],
        sys.argv[2],
    )
