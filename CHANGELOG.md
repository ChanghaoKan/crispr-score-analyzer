# Changelog

All notable changes to **CRISPR Score Analyzer** are documented in this file.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Installable Python library (`crispr_score_analyzer`) via PEP 621 `pyproject.toml` (`pip install -e .` without requiring Streamlit; `pip install -e ".[app,dev]"` for the full web UI and test suite).
- Streamlit-free plotting module (`crispr_score_analyzer.plotting`) exposing CNS/publication Nature NPG and 7-stop ColorBrewer `RdBu` Plotly figure builders for scripts and Jupyter notebooks.
- Command-line interface (`crispr-score-analyzer` / `python -m crispr_score_analyzer`) supporting `info`, `rank`, `codep` (single-gene genome-wide scan & multi-gene matrix), and `app`, with vector export and provenance ZIP generation.
- Deterministic `filter_gene_suggestions` helper for exact/prefix/HGNC-alias gene lookup.
- Developer governance and documentation: `CONTRIBUTING.md`, `CHANGELOG.md`, `docs/api.md`, and JOSS manuscript draft (`paper/paper.md`, `paper/paper.bib`).
- Package API and CLI test suites (`tests/test_package_api.py`, `tests/test_cli.py`).

### Changed
- Upgraded visual design across Light (`#f5f5f7`) and deep-charcoal Dark (`#0f0f11`) themes with colorblind-safe Nature Publishing Group (NPG) categorical hues and 7-stop ColorBrewer `RdBu` diverging scales.
- Replaced fragile 18,434-option combobox dropdowns on Co-dependency `Query gene` and the `Search` popover with native `st.text_input` fields plus deterministic prefix/alias suggestion pills, eliminating dropdown hijacking and IME digit-`1` collisions when typing numbered symbols such as `TEAD1`.
- Preserved backward-compatible root shims (`analysis_core.py`, `codependency.py`, `figure_export.py`) and root `app.py` / `requirements.txt` for Streamlit Community Cloud.

### Fixed
- Fixed typing numbered gene symbols (e.g., `TEAD1`) where combobox dropdown selection or IME interactions could substitute the first dropdown item instead of the trailing digit `1`.
- Fixed the Co-dependency `Clear` button not resetting `codep_query` cleanly when uncommitted combobox text was present.
- Fixed dark-mode tab header and caption contrast on Streamlit 1.67+.

## [1.2.0] - 2026-11-23

Zenodo software archive: [10.5281/zenodo.24458535](https://doi.org/10.5281/zenodo.24458535) (MIT source; excludes bundled DepMap matrix via `.gitattributes`).

### Added
- Standalone **Co-dependency** module with chunked `float64` genome-wide Pearson/Spearman single-gene scans (including cohort-scoped single-lineage scans), 2–25 gene correlation matrices, and integrated lineage-colored cell-line scatter plots.
- **Lineage preference** diverging lollipop / heatmap deep-dive tab and **Top cell lines** paired sensitive/resistant model lookup tab.
- Offline HGNC alias resolution (`data/gene_aliases.json`, 40,402 unambiguous mappings) and shareable URL query parameters.
- Explicit Broad Institute DepMap Terms and Conditions notice (`data/DEPMAP_TERMS.md`) and `.gitattributes` exclusion of `data/depmap_chronos.npz` from MIT release archives.

## [1.1.0] - 2026-11-21

Zenodo software archive: [10.5281/zenodo.24424595](https://doi.org/10.5281/zenodo.24424595) (MIT).

### Added
- Updated bundled reference dataset to DepMap Public 26Q3 Chronos Gene Effect (1,361 cell lines × 18,434 genes across 30 Oncotree lineages) with SHA-256/MD5 provenance (`data/depmap_chronos.provenance.json`).
- Strict `SYMBOL (Entrez ID)` schema detection, custom CSV column diagnostics (`diagnostics.csv`), deterministic white-paper PDF/SVG/300-DPI PNG figure exports, and on-demand analysis provenance ZIP bundles.

## [1.0.0] - 2026-04-16

Zenodo software archive: [10.5281/zenodo.19607603](https://doi.org/10.5281/zenodo.19607603) (CC-BY-4.0).

### Added
- Initial public Streamlit web application for DepMap CRISPR gene essentiality ranking, cancer-type boxplot comparison, and multi-layer gene-set annotation.
