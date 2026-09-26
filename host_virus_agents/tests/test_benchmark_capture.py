import csv
import io
import json
from contextlib import redirect_stdout
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import benchmark_capture
from evaluate_pairs import run_benchmark


class CaptureTests(unittest.TestCase):
    def setUp(self):
        # Offline tests must not hash real staged model weights on Nibi.
        fixture_model = patch.dict("os.environ", {"MODEL_PATH": "tests/_absent_model_fixture"})
        fixture_model.start()
        self.addCleanup(fixture_model.stop)

    def test_duplicate_rows_and_failure_retained(self):
        with tempfile.TemporaryDirectory(dir=".") as directory:
            source = Path(directory) / "input.csv"
            output = Path(directory) / "output.csv"
            source.write_text("host,virus,expected_status\nExample animal,Example virus,KNOWN\nExample animal,Example virus,KNOWN\n", encoding="utf-8")
            with patch("evaluate_pairs.run_judge_agent", side_effect=[dict(literature_status="KNOWN", classification="KNOWN"), RuntimeError("timeout")]), redirect_stdout(io.StringIO()):
                run_benchmark(source, output)
            with output.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[1]["classification"], "INSUFFICIENT_EVIDENCE")
            records = [json.loads(line) for line in Path(str(output) + ".diagnostics.jsonl").read_text(encoding="utf-8").splitlines()]
            self.assertEqual([r["kind"] for r in records], ["manifest", "pair", "pair"])
            self.assertEqual(records[-1]["diagnostics"]["error"], "timeout")
            with self.assertRaises(FileExistsError):
                benchmark_capture.begin(source, output)

    def test_environment_failure_is_captured(self):
        from evaluate_env import run_environment_benchmark
        from validate_live_benchmark import validate
        with tempfile.TemporaryDirectory(dir=".") as directory:
            source, output = Path(directory) / "input.csv", Path(directory) / "env.csv"
            source.write_text("host,virus,expected_status\nExample animal,Example virus,KNOWN\n", encoding="utf-8")
            with patch("evaluate_env.BioResearchEnv", side_effect=RuntimeError("retrieval failed")), \
                 patch("evaluate_env.BaselineResearchAgent"), redirect_stdout(io.StringIO()):
                run_environment_benchmark(source, output)
            result = validate(str(output) + ".diagnostics.jsonl", source)
            self.assertEqual(result["correct"], 0)
            self.assertFalse(result["automated_checks_passed"])
            self.assertEqual(result["cases"][0]["predicted"], "UNCLEAR")


if __name__ == "__main__":
    unittest.main()
