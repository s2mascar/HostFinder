"""Exercise network-independent function bodies against literal XML fixtures."""
import ast
from pathlib import Path
import re
import unittest
from unittest.mock import Mock
import xml.etree.ElementTree as ET


def load_functions(path, names, namespace):
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    module = ast.Module(body=[n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names], type_ignores=[])
    exec(compile(module, path, "exec"), namespace)
    return namespace


class TaxonomyTests(unittest.TestCase):
    def test_virus_search_does_not_choose_unmatched_first_hit(self):
        ns = load_functions("bioresearch_env/virus_taxonomy_agent.py", {"_fetch_by_name", "_parse_taxon", "_all_name_texts"},
                            dict(ET=ET, normalize_name=lambda x: x.lower(), EMAIL="", TOOL=""))
        response = Mock(text="<TaxaSet><Taxon><TaxId>12</TaxId><ScientificName>Different virus</ScientificName><Rank>species</Rank></Taxon></TaxaSet>")
        response.json.return_value = {"esearchresult": {"idlist": ["12"]}}
        ns["ncbi_get"] = Mock(return_value=response)
        self.assertFalse(ns["_fetch_by_name"]("Example virus")["resolved"])

    def test_host_aliases_do_not_choose_unmatched_first_hit(self):
        ns = load_functions("taxonomy_aliases.py", {"fetch_taxonomy_aliases", "scientific_abbreviation"},
                            dict(ET=ET, normalize_name=lambda x: x.lower(), EMAIL="", TOOL=""))
        response = Mock(text="<TaxaSet><Taxon><ScientificName>Different animal</ScientificName></Taxon></TaxaSet>")
        response.json.return_value = {"esearchresult": {"idlist": ["12"]}}
        ns["ncbi_get"] = Mock(return_value=response)
        self.assertEqual(ns["fetch_taxonomy_aliases"]("Example animal"), ["Example animal"])


if __name__ == "__main__":
    unittest.main()
