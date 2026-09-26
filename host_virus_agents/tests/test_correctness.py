"""Synthetic assertion tests; names are unrelated to benchmark examples."""
import unittest

from evidence_semantics import classify, decide, entity_equal


class BindingTests(unittest.TestCase):
    def check_edge(self, text, expected, **fields):
        extraction = dict(study_host="Example animal", host_virus_name="Example virus 1",
                          host_virus_passage=text, study_host_passage=text,
                          host_virus_relationship_type="DETECTED_IN",
                          extraction_status="SUCCESS")
        extraction.update(fields)
        result = classify(extraction, ["Example animal"], "Example virus 1", text)
        self.assertEqual(result, expected, extraction)
        return extraction

    def test_exact_detection(self):
        self.check_edge("Example virus 1 was detected in Example animal.", "EXACT_SUPPORT")

    def test_other_host(self):
        self.check_edge("Example virus 1 was detected in Different animal.", "VIRUS_OTHER_HOST",
                        study_host="Different animal")

    def test_other_virus(self):
        self.check_edge("Different virus was detected in Example animal.", "NO_SUPPORT",
                        host_virus_name="Different virus")

    def test_comparison(self):
        self.check_edge("Different virus was detected in Example animal and was similar to Example virus 1.", "UNCLEAR")

    def test_similarity_only(self):
        self.check_edge("Example animal protein is similar to Example virus 1 protein.", "NO_SUPPORT",
                        host_virus_relationship_type="SEQUENCE_SIMILARITY")

    def test_environment(self):
        e = self.check_edge("Example virus 1 was detected in wastewater near Example animal.", "NO_SUPPORT")
        self.assertEqual(e["structured_evidence"]["evidence_scope"], "ENVIRONMENTAL_ASSOCIATION")

    def test_experimental(self):
        e = self.check_edge("Example animal was experimentally infected with Example virus 1.", "EXACT_SUPPORT")
        self.assertEqual(e["structured_evidence"]["natural_vs_experimental"], "EXPERIMENTAL")

    def test_natural(self):
        e = self.check_edge("Example animal was naturally infected with Example virus 1.", "EXACT_SUPPORT")
        self.assertEqual(e["structured_evidence"]["natural_vs_experimental"], "NATURAL")

    def test_negation(self):
        self.check_edge("Example virus 1 was not detected in Example animal.", "NO_SUPPORT")

    def test_exposure_not_infection(self):
        self.check_edge("Example animal was inoculated with Example virus 1.", "UNCLEAR")

    def test_antigen_not_virus(self):
        self.check_edge("Example virus 1 protein was detected in Example animal.", "NO_SUPPORT")

    def test_background(self):
        self.check_edge("Previous studies reported Example virus 1 was detected in Example animal.", "NO_SUPPORT")

    def test_rank(self):
        self.check_edge("Example virus 1 was detected in Example animal.", "NO_SUPPORT", host_scope="GENUS")

    def test_numbered_virus_not_alias(self):
        self.assertFalse(entity_equal("Example virus 1", "Example virus 16"))

    def test_alias(self):
        self.assertTrue(entity_equal("Old designation", "Example virus 1", ["Old designation"]))

    def test_shortened_host_not_inferred(self):
        self.check_edge("Example virus 1 was detected in E. animal.", "UNCLEAR")

    def test_failed_extraction(self):
        self.check_edge("Example virus 1 was detected in Example animal.", "UNCLEAR", extraction_status="FAILED")

    def test_source_conflict(self):
        self.check_edge("Example virus 1 was detected in Example animal.", "UNCLEAR",
                        comparison_source_virus_name="Different virus", comparison_virus_name="Another virus")

    def test_related_chain(self):
        text = "Different virus was detected in Example animal. Different virus is similar to Example virus 1."
        e = self.check_edge(text, "TARGET_HOST_RELATED", host_virus_name="Different virus",
                           comparison_source_virus_name="Different virus", comparison_virus_name="Example virus 1",
                           comparison_relationship_type="SEQUENCE_SIMILARITY",
                           comparison_relationship_passage="Different virus is similar to Example virus 1.")
        self.assertTrue(e["comparison_edge"]["verified"])

    def test_specimen_linkage(self):
        text = "Specimen batch-Q7 was collected from Example animal. Example virus 1 was detected in specimen batch-Q7."
        self.check_edge(text, "EXACT_SUPPORT", shared_specimen="batch-Q7",
                        study_host_passage="Specimen batch-Q7 was collected from Example animal.",
                        host_virus_passage="Example virus 1 was detected in specimen batch-Q7.")

    def test_unlinked_discovery(self):
        self.check_edge("Example animal was studied in a virome survey. Example virus 1 genome is 1400 nucleotides long.", "UNCLEAR")

    def test_explicit_local_acronym(self):
        self.check_edge("Example virus 1 (EV1) was examined. EV1 was detected in Example animal.", "EXACT_SUPPORT")

    def test_alias_abbreviation_ambiguous(self):
        e = dict(study_host="Example animal", host_virus_name="Example virus 1",
                 host_virus_passage="Example virus 1 was detected in E. animal.")
        self.assertEqual(classify(e, ["Example animal", "E. animal"], "Example virus 1", e["host_virus_passage"]), "UNCLEAR")

    def test_hypothesis_is_not_observation(self):
        self.check_edge("We asked whether Example virus 1 was detected in Example animal.", "UNCLEAR")

    def test_fecal_source_is_not_host_proof(self):
        self.check_edge("Example virus 1 was detected in Example animal feces.", "NO_SUPPORT")

    def test_cell_line_is_not_whole_host(self):
        self.check_edge("Example virus 1 was detected in Example animal cell cultures.", "NO_SUPPORT")

    def test_host_name_inside_virus_name(self):
        self.check_edge("Example virus 1 was detected in Example animal virus preparations.", "UNCLEAR")

    def test_failed_comparison_does_not_create_known(self):
        self.check_edge("Different virus was detected in Example animal, unlike Example virus 1.", "UNCLEAR")


class AggregationTests(unittest.TestCase):
    complete = dict(retrieval_complete=True, source_successes=3, source_failures=[],
                    queries_planned=1, queries_executed=1, taxonomy_resolved=True)

    def test_search_failure(self):
        self.assertEqual(decide([], dict(self.complete, source_failures=["timeout"]))["classification"], "INSUFFICIENT_EVIDENCE")

    def test_zero_is_not_novel(self):
        self.assertEqual(decide([], self.complete)["classification"], "INSUFFICIENT_EVIDENCE")

    def test_unresolved_taxonomy(self):
        self.assertEqual(decide([], dict(self.complete, taxonomy_resolved=False))["classification"], "INSUFFICIENT_EVIDENCE")

    def test_label_alone_cannot_support(self):
        self.assertNotEqual(decide([{"classification": "EXACT_SUPPORT"}], self.complete)["classification"], "KNOWN")

    def test_exact_dominates_background_and_failed_candidates(self):
        text = "Example virus 1 was detected in Example animal."
        e = dict(study_host="Example animal", host_virus_name="Example virus 1", host_virus_passage=text)
        classify(e, ["Example animal"], "Example virus 1", text)
        self.assertEqual(decide([{"extraction_status": "FAILED"}, e], self.complete)["classification"], "KNOWN")

    def test_related_never_exact(self):
        text = "Other virus was detected in Example animal. Other virus is similar to Example virus 1."
        e = dict(study_host="Example animal", host_virus_name="Other virus", host_virus_passage=text,
                 comparison_source_virus_name="Other virus", comparison_virus_name="Example virus 1",
                 comparison_relationship_passage="Other virus is similar to Example virus 1.")
        classify(e, ["Example animal"], "Example virus 1", text)
        self.assertEqual(decide([e], self.complete)["literature_status"], "POSSIBLY_KNOWN")

    def test_partial_success_is_incomplete(self):
        self.assertEqual(decide([], dict(self.complete, queries_planned=3, queries_executed=2))["literature_status"], "UNCLEAR")


if __name__ == "__main__":
    unittest.main()
