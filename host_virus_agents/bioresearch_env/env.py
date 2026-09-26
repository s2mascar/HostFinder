from copy import deepcopy

from evidence_agent import analyze_paper, get_pair_taxonomy

from literature_search import (
    merge_papers,
)

from search_agent import (
    search_query_sources,
    candidate_score,
)

from taxonomy_aliases import (
    get_host_aliases,
)

from bioresearch_env.actions import (
    SEARCH,
    ANALYZE_PAPER,
    SUBMIT,
    validate_action,
)

from bioresearch_env.rewards import (
    action_cost,
    terminal_reward,
    MAX_STEPS_PENALTY,
)

from bioresearch_env.biological_context_agent import (
    run_biological_context_agent,
)


class BioResearchEnv:
    """
    Evidence-grounded host-virus research environment.

    Episode phases:
        1. SEARCH
        2. ANALYZE
        3. SUBMIT

    All search results are pooled before any paper is analyzed. They are
    then deduplicated, globally ranked, and truncated to max_papers.

    expected_status is hidden from observations and never shown to the
    agent.
    """

    def __init__(
        self,
        task,
        max_steps=15,
        max_papers=8,
    ):
        self.task = deepcopy(task)

        self.host = task["host"]
        self.virus = task["virus"]
        self._expected_status = task["expected_status"]

        self.max_steps = max_steps
        self.max_papers = max_papers

        self.host_aliases = []
        self.biological_context = {}

        self.raw_candidate_papers = []
        self.papers = {}
        self.analyzed_results = {}

        self.search_metadata = {}
        self.history = []

        self.steps = 0
        self.search_count = 0
        self.analysis_count = 0
        self.total_reward = 0.0

        self.phase = "SEARCH"
        self.search_finalized = False
        self.done = False
        self.predicted_status = None

        self.last_action_result = None
        self.final_info = {}

    # ========================================================
    # RESET
    # ========================================================

    def reset(self):
        self.host_aliases = get_host_aliases(
            self.host
        )

        self.biological_context = (
            run_biological_context_agent(
                self.host,
                self.virus,
            )
        )

        self.raw_candidate_papers = []
        self.papers = {}
        self.analyzed_results = {}

        self.search_metadata = {
            "queries_attempted": [],
            "queries_planned": 0,
            "queries_executed": 0,
            "source_successes": 0,
            "source_failures": [],
            "candidate_papers": 0,
            "retrieval_complete": True,
            "host_aliases": self.host_aliases,
        }
        self.search_metadata["taxonomy_resolution"] = get_pair_taxonomy(self.host, self.virus)
        self.search_metadata["taxonomy_resolved"] = self.search_metadata["taxonomy_resolution"]["resolved"]

        self.history = []

        self.steps = 0
        self.search_count = 0
        self.analysis_count = 0
        self.total_reward = 0.0

        self.phase = "SEARCH"
        self.search_finalized = False
        self.done = False
        self.predicted_status = None

        self.last_action_result = None
        self.final_info = {}

        return self._observation()

    # ========================================================
    # SEARCH
    # ========================================================

    def _search(self, action):
        if self.phase != "SEARCH":
            raise RuntimeError(
                "SEARCH is only allowed before paper analysis starts."
            )

        host_term = action["host_term"]
        virus_term = action["virus_term"]

        papers, statuses = search_query_sources(
            host_term,
            virus_term,
        )

        self.search_count += 1
        self.search_metadata["queries_executed"] += 1
        self.search_metadata["queries_attempted"].append({"host_term": host_term, "virus_term": virus_term,
                                                         "source_status": statuses, "retrieved_papers": papers})

        for status in statuses:
            if status.get("success"):
                self.search_metadata["source_successes"] += 1
            else:
                self.search_metadata["retrieval_complete"] = False
                self.search_metadata["source_failures"].append({
                    "host_term": host_term,
                    "virus_term": virus_term,
                    "source": status.get("source"),
                    "error": status.get("error"),
                })

        self.raw_candidate_papers.extend(
            papers
        )

        # Search results have changed, so global selection must be
        # recomputed before analysis begins.
        self.search_finalized = False

        result = {
            "action": SEARCH,
            "host_term": host_term,
            "virus_term": virus_term,
            "papers_returned": len(papers),
            "raw_candidate_papers": len(
                self.raw_candidate_papers
            ),
            "retrieved_papers": [
                self._raw_paper_summary(paper)
                for paper in papers[:10]
            ],
        }

        return result

    # ========================================================
    # GLOBAL PAPER SELECTION
    # ========================================================

    def finalize_searches(self):
        """
        Match the original HostFinder search behavior:

        pool all searches -> merge duplicates -> globally rank -> top N
        """

        if self.search_finalized:
            return

        unique_papers = merge_papers(
            self.raw_candidate_papers
        )

        unique_papers = sorted(
            unique_papers,
            key=lambda paper: candidate_score(
                paper,
                self.host_aliases,
                self.virus,
                host=self.host,
            ),
            reverse=True,
        )

        selected = unique_papers[
            :self.max_papers
        ]

        self.papers = {}

        for index, paper in enumerate(
            selected
        ):
            paper_id = f"P{index}"

            paper_copy = deepcopy(
                paper
            )

            paper_copy["paper_id"] = paper_id
            paper_copy["candidate_score"] = candidate_score(
                paper_copy,
                self.host_aliases,
                self.virus,
                host=self.host,
            )

            self.papers[paper_id] = paper_copy

        self.search_metadata["candidate_papers"] = len(
            self.papers
        )

        self.search_finalized = True
        self.phase = "ANALYZE"

    # ========================================================
    # ANALYZE PAPER
    # ========================================================

    def _analyze_paper(self, action):
        if not self.search_finalized:
            self.finalize_searches()

        paper_id = action["paper_id"]

        if paper_id not in self.papers:
            raise ValueError(
                f"Unknown paper_id: {paper_id}"
            )

        if paper_id in self.analyzed_results:
            return {
                "action": ANALYZE_PAPER,
                "paper_id": paper_id,
                "cached": True,
                "result": self.analyzed_results[paper_id],
            }

        paper = self.papers[paper_id]

        result = analyze_paper(
            self.host,
            self.host_aliases,
            self.virus,
            paper,
            biological_context=self.biological_context,
        )

        self.analyzed_results[paper_id] = result
        self.analysis_count += 1

        return {
            "action": ANALYZE_PAPER,
            "paper_id": paper_id,
            "cached": False,
            "result": result,
        }

    # ========================================================
    # SUBMIT
    # ========================================================

    def _submit(self, action):
        predicted_status = action["status"]

        supporting_paper_ids = (
            action.get(
                "supporting_paper_ids",
                [],
            )
            or []
        )

        self.predicted_status = predicted_status

        reward = terminal_reward(
            predicted_status=predicted_status,
            expected_status=self._expected_status,
            supporting_paper_ids=supporting_paper_ids,
            analyzed_results=self.analyzed_results,
        )

        self.done = True
        self.phase = "DONE"

        correct = (
            predicted_status
            == self._expected_status
        )

        self.final_info = {
            "predicted_status": predicted_status,
            "expected_status": self._expected_status,
            "correct": correct,
            "supporting_paper_ids": supporting_paper_ids,
            "steps": self.steps,
            "searches": self.search_count,
            "papers_retrieved": len(self.papers),
            "papers_analyzed": len(
                self.analyzed_results
            ),
        }

        return (
            {
                "action": SUBMIT,
                "predicted_status": predicted_status,
                "correct": correct,
            },
            reward,
        )

    # ========================================================
    # STEP
    # ========================================================

    def step(self, action):
        if self.done:
            raise RuntimeError(
                "Episode has already finished."
            )

        validate_action(action)

        self.steps += 1
        action_type = action["type"]

        reward = action_cost(
            action_type
        )

        if action_type == SEARCH:
            result = self._search(action)

        elif action_type == ANALYZE_PAPER:
            result = self._analyze_paper(action)

        elif action_type == SUBMIT:
            result, submit_reward = self._submit(
                action
            )
            reward += submit_reward

        else:
            raise ValueError(
                f"Unsupported action: {action_type}"
            )

        self.last_action_result = deepcopy(
            result
        )

        if (
            not self.done
            and self.steps >= self.max_steps
        ):
            self.done = True
            self.phase = "DONE"
            reward += MAX_STEPS_PENALTY

            self.final_info = {
                "predicted_status": None,
                "expected_status": self._expected_status,
                "correct": False,
                "supporting_paper_ids": [],
                "steps": self.steps,
                "searches": self.search_count,
                "papers_retrieved": len(self.papers),
                "papers_analyzed": len(
                    self.analyzed_results
                ),
                "terminated_reason": "MAX_STEPS",
            }

        self.total_reward += reward

        self.history.append({
            "step": self.steps,
            "action": deepcopy(action),
            "reward": reward,
            "result": deepcopy(result),
        })

        observation = self._observation()
        info = {}

        if self.done:
            info = deepcopy(
                self.final_info
            )
            info["episode_reward"] = self.total_reward

        return (
            observation,
            reward,
            self.done,
            info,
        )

    # ========================================================
    # OBSERVATIONS
    # ========================================================

    def _raw_paper_summary(self, paper):
        return {
            "title": paper.get("title"),
            "pmid": paper.get("pmid"),
            "pmcid": paper.get("pmcid"),
            "doi": paper.get("doi"),
            "sources": paper.get("sources", []),
        }

    def _paper_summaries(self):
        summaries = []

        for paper_id, paper in self.papers.items():
            summaries.append({
                "paper_id": paper_id,
                "title": paper.get("title"),
                "pmid": paper.get("pmid"),
                "pmcid": paper.get("pmcid"),
                "doi": paper.get("doi"),
                "sources": paper.get("sources", []),
                "candidate_score": paper.get(
                    "candidate_score"
                ),
                "analyzed": (
                    paper_id
                    in self.analyzed_results
                ),
            })

        return summaries

    def _safe_biological_context(self):
        return {
            key: deepcopy(value)
            for key, value in self.biological_context.items()
            if not key.startswith("_")
        }

    def _observation(self):
        return {
            "host": self.host,
            "virus": self.virus,
            "phase": self.phase,
            "steps_used": self.steps,
            "steps_remaining": max(
                0,
                self.max_steps - self.steps,
            ),
            "available_actions": (
                [SEARCH, ANALYZE_PAPER, SUBMIT]
                if self.phase == "SEARCH"
                else [ANALYZE_PAPER, SUBMIT]
                if self.phase == "ANALYZE"
                else []
            ),
            "biological_context": self._safe_biological_context(),
            "papers": self._paper_summaries(),
            "analyzed_results": deepcopy(
                self.analyzed_results
            ),
            "last_action_result": deepcopy(
                self.last_action_result
            ),
        }

    # ========================================================
    # PUBLIC HELPERS
    # ========================================================

    def get_search_metadata(self):
        return deepcopy(
            self.search_metadata
        )

    def get_analyzed_results(self):
        return deepcopy(
            self.analyzed_results
        )

    def get_paper_ids(self):
        if not self.search_finalized:
            self.finalize_searches()

        return list(
            self.papers.keys()
        )

    def get_biological_context(self):
        return deepcopy(
            self.biological_context
        )
