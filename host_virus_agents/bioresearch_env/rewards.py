SEARCH_COST = -0.05
ANALYZE_COST = -0.02

CORRECT_REWARD = 1.0
INCORRECT_REWARD = -1.0

GOOD_EVIDENCE_BONUS = 0.50
UNSUPPORTED_EVIDENCE_PENALTY = -0.50

MAX_STEPS_PENALTY = -1.0


def action_cost(action_type):
    if action_type == "SEARCH":
        return SEARCH_COST

    if action_type == "ANALYZE_PAPER":
        return ANALYZE_COST

    return 0.0


def evidence_reward(
    predicted_status,
    supporting_paper_ids,
    analyzed_results,
):
    supporting_paper_ids = (
        supporting_paper_ids
        or []
    )

    if predicted_status == "KNOWN":
        if not supporting_paper_ids:
            return UNSUPPORTED_EVIDENCE_PENALTY

        for paper_id in supporting_paper_ids:
            result = analyzed_results.get(
                paper_id,
                {},
            )

            if result.get("classification") == "EXACT_SUPPORT":
                return GOOD_EVIDENCE_BONUS

        return UNSUPPORTED_EVIDENCE_PENALTY

    if predicted_status == "POSSIBLY_KNOWN":
        if not supporting_paper_ids:
            return UNSUPPORTED_EVIDENCE_PENALTY

        for paper_id in supporting_paper_ids:
            result = analyzed_results.get(
                paper_id,
                {},
            )

            if result.get("classification") == "TARGET_HOST_RELATED":
                return GOOD_EVIDENCE_BONUS

        return UNSUPPORTED_EVIDENCE_PENALTY

    return 0.0


def terminal_reward(
    predicted_status,
    expected_status,
    supporting_paper_ids,
    analyzed_results,
):
    reward = 0.0

    if predicted_status == expected_status:
        reward += CORRECT_REWARD
    else:
        reward += INCORRECT_REWARD

    reward += evidence_reward(
        predicted_status,
        supporting_paper_ids,
        analyzed_results,
    )

    return reward
