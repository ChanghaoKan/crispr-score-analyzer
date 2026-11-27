# Contributing to CRISPR Score Analyzer

Contributions, bug reports, documentation improvements, and scientific workflow feedback are welcome.

## Getting help and reporting issues

- **Bug reports & feature requests:** Open an issue on [GitHub Issues](https://github.com/ChanghaoKan/crispr-score-analyzer/issues).
- **Direct contact:** Changhao Kan (`kanchanghao@szbl.ac.cn`), Institute of Molecular Physiology, Shenzhen Bay Laboratory.

When reporting a bug, please include:
1. Python and dependency versions (`crispr-score-analyzer info` or `pip list`).
2. Steps to reproduce (and a small synthetic CSV snippet when reporting a custom-data issue).
3. Expected vs. observed behavior (or a screenshot for UI issues).

## Development setup

Clone the repository and install the package in editable mode with both the Streamlit UI and test dependencies:

```bash
git clone https://github.com/ChanghaoKan/crispr-score-analyzer.git
cd crispr-score-analyzer
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[app,dev]"
```

Run the full automated test suite (runs offline using deterministic synthetic fixtures; does not download DepMap files):

```bash
python -m pytest tests -v
```

To skip the Kaleido Chromium/binary image rendering test in restricted sandboxes:

```bash
python -m pytest tests -q -k "not actual_svg_is_vector_pdf_is_valid_and_png_has_dpi"
```

Launch the local web interface:

```bash
crispr-score-analyzer app
# or directly:
python -m streamlit run app.py
```

## Architecture guidelines

The repository separates reusable scientific computation from interactive UI state:

- `crispr_score_analyzer/analysis_core.py`: CSV/NPZ loaders, strict schema validation, HGNC alias resolution, gene ranking, lineage preference, correlation matrix, and provenance ZIP bundling.
- `crispr_score_analyzer/codependency.py`: Chunked `float64` pairwise-complete Pearson/Spearman genome-wide co-dependency scanning and cell-line scatter builders.
- `crispr_score_analyzer/plotting.py`: Streamlit-free Plotly builders using colorblind-safe Nature Publishing Group (NPG) and 7-stop ColorBrewer `RdBu` palettes.
- `crispr_score_analyzer/figure_export.py`: Deterministic white-paper PDF/SVG/300-DPI PNG figure export and fingerprinting.
- `crispr_score_analyzer/cli.py`: Reproducible command-line entrypoint (`crispr-score-analyzer`).
- `app.py`: Streamlit interactive web frontend (kept thin so scientific routines remain importable without Streamlit).
- Root `analysis_core.py`, `codependency.py`, `figure_export.py`: Thin backward-compatible re-export shims.

## Pull request checklist

1. Keep analytical functions in `crispr_score_analyzer/` free of `streamlit` imports.
2. Add or update unit/integration tests under `tests/` for any behavioral change or bug fix.
3. Verify `python -m pytest tests -q` passes locally before opening a pull request.
4. Keep software code and third-party data licensing strictly separate: all contributed code must be compatible with the [MIT License](License), and pull requests must **not** bundle proprietary or non-commercial third-party datasets into Python distribution artifacts (`pyproject.toml` excludes `data/depmap_chronos.npz`).
5. Preserve clear scientific framing: outputs are descriptive summaries and must not claim formal hypothesis-testing significance where none is computed.
