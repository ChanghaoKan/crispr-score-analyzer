"""The current app cites the 1.1.0 Zenodo deposit, not the April 2026 snapshot."""

import ast
from pathlib import Path


def test_citation_uses_the_1_1_0_zenodo_deposit():
    names = {"build_citations"}
    constants = {
        "TOOL_VERSION", "TOOL_RELEASED", "TOOL_AUTHORS", "TOOL_YEAR", "DATA_VERSION",
        "ZENODO_CONCEPT_DOI", "ZENODO_VERSION_DOI", "ZENODO_VERSION", "ZENODO_RELEASED",
        "ZENODO_PREVIOUS_DOI", "ZENODO_PREVIOUS_VERSION", "ZENODO_PREVIOUS_RELEASED",
        "GITHUB_URL",
    }
    tree = ast.parse((Path(__file__).resolve().parents[1] / "app.py").read_text())
    nodes = [node for node in tree.body
             if (isinstance(node, ast.FunctionDef) and node.name in names)
             or (isinstance(node, ast.Assign)
                 and any(isinstance(target, ast.Name) and target.id in constants
                         for target in node.targets))]
    namespace = {}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), "app.py", "exec"), namespace)
    current, previous, bibtex = namespace["build_citations"]()
    assert current.startswith("Kan, C. (2026). CRISPR Score Analyzer (1.1.0)")
    assert "10.5281/zenodo.24424595" in current
    assert "1.0.0" in previous
    assert "10.5281/zenodo.19607603" in previous
    assert "version = {1.1.0}" in bibtex
    assert "10.5281/zenodo.24424595" in bibtex
    assert "19607603" not in bibtex
