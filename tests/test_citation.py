"""The running app must not be labeled as the Zenodo 1.0.0 deposit."""

import ast
from pathlib import Path


def test_running_version_is_separate_from_the_zenodo_deposit():
    names = {"build_citations"}
    constants = {
        "TOOL_VERSION", "TOOL_RELEASED", "TOOL_AUTHORS", "TOOL_YEAR", "DATA_VERSION",
        "ZENODO_CONCEPT_DOI", "ZENODO_VERSION_DOI", "ZENODO_VERSION", "ZENODO_RELEASED",
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
    running, archived, bibtex = namespace["build_citations"]()
    assert "1.1.0" in running
    assert "DepMap Public 26Q3" in running
    assert "1.0.0" in archived
    assert "10.5281/zenodo.19607603" in archived
    assert "version = {1.0.0}" in bibtex
    assert "10.5281/zenodo.19607603" in bibtex
    assert "1.1.0" not in bibtex
