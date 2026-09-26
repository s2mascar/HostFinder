"""Actual public entry points, with no model/network required."""
import unittest
from unittest.mock import patch


class IntegrationTests(unittest.TestCase):
    def test_public_classifier_and_judge(self):
        from evidence_agent import classify_relationship
        from judge_agent import judge_interaction
        text = "Example virus 1 was detected in Example animal."
        result = dict(study_host="Example animal", host_virus_name="Example virus 1",
                      host_virus_passage=text, host_virus_relationship_type="DETECTED_IN",
                      extraction_status="SUCCESS")
        self.assertEqual(classify_relationship(result, ["Example animal"], "Example virus 1", text), "EXACT_SUPPORT")
        decision = judge_interaction("Example animal", "Example virus 1", [result], {"taxonomy_resolved": True})
        self.assertEqual(decision["classification"], "KNOWN")

    def test_rescue_role_mismatch(self):
        from evidence_agent import _extraction_rescue_reason
        e = dict(study_host="Example animal", host_virus_name="Virus A",
                 comparison_source_virus_name="Virus B", comparison_virus_name="Virus C")
        self.assertEqual(_extraction_rescue_reason(e, "", "Virus C", ["Virus C"]), "COMPARISON_SOURCE_MISMATCH")

    def test_failed_model_reaches_judge_as_uncertainty(self):
        from evidence_agent import analyze_paper
        from judge_agent import judge_interaction
        with patch("evidence_agent.get_virus_taxonomy_context", return_value={"aliases": [], "resolved": True}), \
             patch("evidence_agent.generate_text", side_effect=RuntimeError("model unavailable")), \
             patch("evidence_agent.cleanup_gpu"):
            result = analyze_paper("Example animal", ["Example animal"], "Example virus 1", {"title": "Example", "abstract": "Example virus 1 was detected in Example animal."})
        self.assertEqual(result["extraction_status"], "FAILED")
        self.assertEqual(result["classification"], "UNCLEAR")
        self.assertEqual(judge_interaction("Example animal", "Example virus 1", [result], {"taxonomy_resolved": True})["classification"], "INSUFFICIENT_EVIDENCE")

    def test_flexible_name_regex(self):
        from evidence_agent import _flexible_name_matches
        self.assertTrue(_flexible_name_matches("Example-virus 1", "Example virus 1"))

    def test_unbound_positive_gets_targeted_rescue(self):
        from evidence_agent import analyze_paper
        primary = dict(study_host="Example animal", host_virus_name="Example virus 1",
                       host_virus_relationship_type="DETECTED_IN", host_virus_passage="Example virus 1")
        rescue = dict(primary, host_virus_passage="Example virus 1 was detected in Example animal.")
        text = "Example virus 1 is a name. Example virus 1 was detected in Example animal."
        with patch("evidence_agent.get_virus_taxonomy_context", return_value={"aliases": [], "resolved": True}), \
             patch("evidence_agent._run_extraction", return_value=primary), \
             patch("evidence_agent._run_rescue_extraction", return_value=rescue) as repair, \
             patch("evidence_agent.cleanup_gpu"):
            result = analyze_paper("Example animal", ["Example animal"], "Example virus 1", {"title": "Discovery", "abstract": text})
        self.assertEqual(result["classification"], "EXACT_SUPPORT")
        # The original fragment also occurs in a non-assertive sentence; the
        # verifier must not cherry-pick a later occurrence for that fragment.
        self.assertEqual(repair.call_count, 1)

    def test_reversed_comparison_requires_model_reextraction(self):
        from evidence_agent import analyze_paper
        text = "Other virus was detected in Example animal. Other virus is similar to Example virus 1."
        primary = dict(study_host="Example animal", host_virus_name="Example virus 1",
                       host_virus_relationship_type="DETECTED_IN", host_virus_passage=text,
                       comparison_source_virus_name="Example virus 1", comparison_virus_name="Other virus",
                       comparison_relationship_passage="Other virus is similar to Example virus 1.")
        rescue = dict(primary, host_virus_name="Other virus", comparison_source_virus_name="Other virus", comparison_virus_name="Example virus 1")
        with patch("evidence_agent.get_virus_taxonomy_context", return_value={"aliases": [], "resolved": True}), \
             patch("evidence_agent._run_extraction", return_value=primary), \
             patch("evidence_agent._run_rescue_extraction", return_value=rescue) as repair, \
             patch("evidence_agent.cleanup_gpu"):
            result = analyze_paper("Example animal", ["Example animal"], "Example virus 1", {"title": "Discovery", "abstract": text})
        self.assertEqual(repair.call_count, 1)
        self.assertEqual(result["classification"], "TARGET_HOST_RELATED")
        self.assertEqual(result["extraction_trace"][0]["output"]["host_virus_name"], "Example virus 1")
        self.assertEqual(result["extraction_trace"][1]["output"]["host_virus_name"], "Other virus")

    def test_full_source_failed_extraction_can_be_independently_rejected(self):
        from evidence_agent import analyze_paper
        text = "Other virus was detected in Different animal."
        with patch("evidence_agent.get_virus_taxonomy_context", return_value={"aliases": [], "resolved": True}), \
             patch("evidence_agent.generate_text", side_effect=RuntimeError("model failed")) as model, \
             patch("evidence_agent.cleanup_gpu"):
            result = analyze_paper("Example animal", ["Example animal"], "Example virus 1", {"title": "Discovery", "full_text": text})
        self.assertEqual(result["extraction_status"], "FAILED")
        self.assertEqual(result["evidence_state"], "IRRELEVANT_OR_REJECTED")
        self.assertEqual(model.call_count, 1)

    def test_judge_summary_prefers_exact_over_other_host(self):
        from evidence_agent import classify_relationship
        from judge_agent import judge_interaction
        papers = []
        for host in ("Different animal", "Example animal"):
            text = f"Example virus 1 was detected in {host}."
            e = dict(study_host=host, host_virus_name="Example virus 1", host_virus_passage=text)
            classify_relationship(e, ["Example animal"], "Example virus 1", text)
            papers.append(e)
        result = judge_interaction("Example animal", "Example virus 1", papers, {"taxonomy_resolved": True})
        self.assertEqual(result["strongest_evidence"], papers[1]["structured_evidence"]["supporting_text"])

    def test_judge_paper_list_does_not_trust_label_alone(self):
        from judge_agent import judge_interaction
        result = judge_interaction("Example animal", "Example virus 1", [{"classification": "EXACT_SUPPORT"}], {"taxonomy_resolved": True})
        self.assertEqual(result["exact_supporting_papers"], [])


if __name__ == "__main__":
    unittest.main()
