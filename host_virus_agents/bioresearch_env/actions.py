SEARCH = "SEARCH"
ANALYZE_PAPER = "ANALYZE_PAPER"
SUBMIT = "SUBMIT"

VALID_ACTIONS = {
    SEARCH,
    ANALYZE_PAPER,
    SUBMIT,
}

VALID_STATUSES = {
    "KNOWN",
    "POSSIBLY_KNOWN",
    "NO_EVIDENCE_FOUND",
    "UNCLEAR",
}


def validate_action(action):
    if not isinstance(action, dict):
        raise TypeError("Action must be a dictionary.")

    action_type = action.get("type")

    if action_type not in VALID_ACTIONS:
        raise ValueError(
            f"Unknown action type: {action_type}"
        )

    if action_type == SEARCH:
        if not action.get("host_term"):
            raise ValueError(
                "SEARCH requires host_term."
            )

        if not action.get("virus_term"):
            raise ValueError(
                "SEARCH requires virus_term."
            )

    elif action_type == ANALYZE_PAPER:
        if not action.get("paper_id"):
            raise ValueError(
                "ANALYZE_PAPER requires paper_id."
            )

    elif action_type == SUBMIT:
        status = action.get("status")

        if status not in VALID_STATUSES:
            raise ValueError(
                f"Invalid submitted status: {status}"
            )

    return True
