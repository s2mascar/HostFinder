"""Generic regression fixtures for over-abstention; no benchmark entities."""
import unittest
from evidence_semantics import classify, decide

COMPLETE = dict(taxonomy_resolved=True, retrieval_complete=True, queries_planned=1,
                queries_executed=1, source_successes=3, source_failures=[])


def extract(text, **kwargs):
    e = dict(study_host="Example animal", host_virus_name="Example virus 1",
             host_virus_relationship_type="DETECTED_IN", host_virus_passage=text,
             evidence_source_complete=True, extraction_status="SUCCESS")
    e.update(kwargs)
    classify(e, ["Example animal"], "Example virus 1", text)
    return e


class MaterialityTests(unittest.TestCase):
    def test_material_failed_extraction(self):
        e = extract("Example virus 1 was recovered from Example animal.", extraction_status="FAILED")
        self.assertEqual(e["evidence_state"], "MATERIAL_UNRESOLVED")
        self.assertEqual(decide([e], COMPLETE)["literature_status"], "UNCLEAR")

    def test_irrelevant_failed_extraction(self):
        e = extract("Other virus was isolated from Different animal.", extraction_status="FAILED")
        self.assertEqual(e["evidence_state"], "IRRELEVANT_OR_REJECTED")
        self.assertEqual(decide([e], COMPLETE)["literature_status"], "NO_EVIDENCE_FOUND")

    def test_absent_target_in_incomplete_source_not_exclusion(self):
        e = extract("Other virus was isolated from Different animal.", extraction_status="FAILED", evidence_source_complete=False)
        self.assertEqual(e["evidence_state"], "MATERIAL_UNRESOLVED")

    def test_receptor_assay_failure(self):
        e = extract("Example virus 1 receptor expression was measured in Example animal cell lines.", extraction_status="FAILED")
        self.assertEqual(e["evidence_state"], "IRRELEVANT_OR_REJECTED")

    def test_background_failure_does_not_hide_direct_evidence(self):
        text = "Example virus 1 receptor expression was measured in Example animal cell lines. Example virus 1 was recovered from Example animal."
        e = extract(text, extraction_status="FAILED")
        self.assertEqual(e["evidence_state"], "MATERIAL_UNRESOLVED")

    def test_materiality_label_cannot_bypass_verification(self):
        e = {"extraction_status": "FAILED", "evidence_state": "IRRELEVANT_OR_REJECTED"}
        self.assertEqual(decide([e], COMPLETE)["literature_status"], "UNCLEAR")

    def test_adjacent_explicit_specimens(self):
        e = extract("We collected specimens from Example animal. Example virus 1 was recovered from these specimens.")
        self.assertEqual(e["classification"], "EXACT_SUPPORT")
        self.assertEqual(e["structured_evidence"]["source_context"], "ADJACENT_SPECIMEN_LINK")

    def test_adjacent_different_sample_arm(self):
        e = extract("We collected specimens from Example animal. Example virus 1 was recovered from other specimens.")
        self.assertNotEqual(e["classification"], "EXACT_SUPPORT")

    def test_named_discovery(self):
        e = extract("A novel virus was recovered from Example animal. This virus was named Example virus 1.")
        self.assertEqual(e["classification"], "EXACT_SUPPORT")

    def test_named_comparison_not_discovery(self):
        e = extract("A novel virus was recovered from Example animal. This virus was similar to Example virus 1.")
        self.assertNotEqual(e["classification"], "EXACT_SUPPORT")

    def test_tested_positive(self):
        self.assertEqual(extract("Example animal tested positive for Example virus 1.")["classification"], "EXACT_SUPPORT")

    def test_ambiguous_acronym(self):
        text = "Example virus 1 (EV1) and Different virus (EV1) were compared. EV1 was recovered from Example animal."
        self.assertNotEqual(extract(text)["classification"], "EXACT_SUPPORT")

    def test_immunization(self):
        self.assertNotEqual(extract("Example animal was immunized with Example virus 1.")["classification"], "EXACT_SUPPORT")

    def test_vector_delivery(self):
        self.assertNotEqual(extract("Vector-mediated expression of Example virus 1 protein was detected in Example animal.")["classification"], "EXACT_SUPPORT")

    def test_taxonomic_inventory(self):
        self.assertNotEqual(extract("Taxonomy list: Example animal; Example virus 1; Different virus.")["classification"], "EXACT_SUPPORT")

    def test_population_detection_list(self):
        self.assertEqual(extract("We found several viruses in Example animal, including Example virus 1 (20%) and Different virus (10%).")["classification"], "EXACT_SUPPORT")

    def test_explicit_host_definition(self):
        self.assertEqual(extract("Example virus 1, the coronavirus of the animal (Example animal), causes disease.")["classification"], "EXACT_SUPPORT")

    def test_passaging_is_experimental(self):
        e = extract("We passaged Example virus 1 in Example animal.")
        self.assertEqual(e["classification"], "EXACT_SUPPORT")
        self.assertEqual(e["structured_evidence"]["natural_vs_experimental"], "EXPERIMENTAL")

    def test_separately_quoted_adjacent_units(self):
        text = "We collected specimens from Example animal. Example virus 1 was recovered from these specimens."
        e = extract(text, study_host_passage="We collected specimens from Example animal.",
                    host_virus_passage="Example virus 1 was recovered from these specimens.")
        self.assertEqual(e["classification"], "EXACT_SUPPORT")

    def test_source_offsets_are_literal(self):
        text = "Introductory sentence.  Example virus 1 was found in Example animal."
        e = extract(text)
        edge = e["structured_evidence"]
        self.assertEqual(text[edge["source_start"]:edge["source_end"]], edge["supporting_text"])

    def test_extracted_acronym_binds_defined_target(self):
        text = "Example virus 1 (EV1) was examined. EV1 was detected in Example animal."
        self.assertEqual(extract(text, host_virus_name="EV1")["classification"], "EXACT_SUPPORT")

    def test_comparison_source_local_acronym(self):
        text = "Other virus (OV1) was detected in Example animal. OV1 is similar to Example virus 1."
        e = extract(text, host_virus_name="Other virus", comparison_source_virus_name="OV1",
                    comparison_virus_name="Example virus 1", comparison_relationship_passage="OV1 is similar to Example virus 1.")
        self.assertEqual(e["classification"], "TARGET_HOST_RELATED")

    def test_repeated_paper_does_not_inflate_confidence(self):
        e = extract("Example virus 1 was detected in Example animal.")
        self.assertEqual(decide([e], COMPLETE)["confidence"], decide([e, e], COMPLETE)["confidence"])

    def test_wrong_background_selection_remains_material(self):
        text = "Previous studies discussed Example virus 1. Example virus 1 was detected in Example animal."
        e = extract(text, host_virus_passage="Previous studies discussed Example virus 1.")
        self.assertEqual(e["evidence_state"], "MATERIAL_UNRESOLVED")

    def test_other_host_does_not_hide_second_target_edge(self):
        text = "Example virus 1 was detected in Different animal. Example virus 1 was detected in Example animal."
        e = extract(text, study_host="Different animal", host_virus_passage="Example virus 1 was detected in Different animal.")
        self.assertEqual(e["evidence_state"], "MATERIAL_UNRESOLVED")

    def test_related_cannot_hide_material_exact_candidate(self):
        text = "Other virus was detected in Example animal. Other virus is similar to Example virus 1."
        e = extract(text, host_virus_name="Other virus", comparison_source_virus_name="Other virus",
                    comparison_virus_name="Example virus 1", comparison_relationship_passage="Other virus is similar to Example virus 1.")
        failed = extract("Example virus 1 was detected in Example animal.", extraction_status="FAILED")
        self.assertEqual(decide([e, failed], COMPLETE)["literature_status"], "UNCLEAR")

    def test_grounded_taxonomy_mention_state(self):
        text = "Taxonomy list: Example animal; Example virus 1; Different virus."
        e = extract(text, host_virus_relationship_type="MENTION_ONLY")
        self.assertEqual(e["evidence_state"], "MENTION_ONLY")


if __name__ == "__main__":
    unittest.main()
