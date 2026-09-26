from search_agent import (
    generate_search_queries,
)

from judge_agent import (
    judge_interaction,
)

from bioresearch_env.actions import (
    SEARCH,
    ANALYZE_PAPER,
    SUBMIT,
)


class BaselineResearchAgent:
    """
    Baseline policy that reuses the current HostFinder search planner,
    evidence analyzer, and deterministic judge while interacting through
    BioResearchEnv actions.

    The agent never sees expected_status.
    """

    def run(self, env):
        observation = env.reset()

        host = observation["host"]
        virus = observation["virus"]

        # ----------------------------------------------------
        # EXISTING HOSTFINDER SEARCH PLAN
        # ----------------------------------------------------

        queries = generate_search_queries(
            host,
            virus,
            env.host_aliases,
        )

        env.search_metadata[
            "queries_planned"
        ] = len(queries)

        for query in queries:
            if env.done:
                break

            action = {
                "type": SEARCH,
                "host_term": query["host_term"],
                "virus_term": query["virus_term"],
            }

            observation, reward, done, info = env.step(
                action
            )

            if done:
                return info

        # ----------------------------------------------------
        # GLOBAL PAPER SELECTION, THEN ANALYZE ALL TOP PAPERS
        # ----------------------------------------------------

        paper_ids = env.get_paper_ids()

        for paper_id in paper_ids:
            if env.done:
                break

            action = {
                "type": ANALYZE_PAPER,
                "paper_id": paper_id,
            }

            observation, reward, done, info = env.step(
                action
            )

            if done:
                return info

        # ----------------------------------------------------
        # EXISTING DETERMINISTIC JUDGE
        # ----------------------------------------------------

        analyzed_map = env.get_analyzed_results()

        evidence_results = [
            result
            for result in analyzed_map.values()
        ]

        search_metadata = env.get_search_metadata()

        judge_result = judge_interaction(
            host,
            virus,
            evidence_results,
            search_metadata,
        )

        predicted_status = judge_result.get(
            "literature_status",
            "UNCLEAR",
        )

        # ----------------------------------------------------
        # SELECT SUPPORTING PAPERS
        # ----------------------------------------------------

        supporting_paper_ids = []

        if predicted_status == "KNOWN":
            for paper_id, result in analyzed_map.items():
                if result.get("classification") == "EXACT_SUPPORT":
                    supporting_paper_ids.append(
                        paper_id
                    )

        elif predicted_status == "POSSIBLY_KNOWN":
            for paper_id, result in analyzed_map.items():
                if result.get("classification") == "TARGET_HOST_RELATED":
                    supporting_paper_ids.append(
                        paper_id
                    )

        # ----------------------------------------------------
        # SUBMIT
        # ----------------------------------------------------

        action = {
            "type": SUBMIT,
            "status": predicted_status,
            "supporting_paper_ids": supporting_paper_ids,
        }

        observation, reward, done, info = env.step(
            action
        )

        info["confidence"] = judge_result.get(
            "confidence",
            "",
        )
        info["classification"] = judge_result["classification"]
        info["judge_diagnostics"] = judge_result

        info["reason"] = judge_result.get(
            "reason",
            "",
        )

        info["exact_supporting_papers"] = len(
            judge_result.get(
                "exact_supporting_papers",
                [],
            )
        )

        info["related_papers"] = len(
            judge_result.get(
                "related_papers",
                [],
            )
        )

        return info
