# Python Library & CLI Guide

`crispr-score-analyzer` separates its computational core, publication plotting builders, vector export helpers, and command-line interface from the Streamlit web frontend (`app.py`). Every analytical workflow can be run headless in Python scripts, Jupyter notebooks, or shell pipelines using either custom CSV score matrices or the optional DepMap Public 26Q3 matrix.

## Installation

Install the Streamlit-free core Python library and CLI:

```bash
git clone https://github.com/ChanghaoKan/crispr-score-analyzer.git
cd crispr-score-analyzer
python -m pip install -e .
```

Install with the interactive Streamlit web frontend and test dependencies:

```bash
python -m pip install -e ".[app,dev]"
```

> **Licensing note:** The MIT-licensed Python package ships the CC0/public-domain HGNC symbol alias table (`40,402` unambiguous mappings) and **never** bundles `data/depmap_chronos.npz`, which is governed separately by Broad Institute DepMap non-commercial terms (`data/DEPMAP_TERMS.md`).

---

## Self-contained Python API quickstart (synthetic CSV)

The example below runs end-to-end without downloading DepMap data:

```python
import io
import crispr_score_analyzer as csa

csv_bytes = b"""ModelID,cell_line_display_name,lineage,MYC (4609),PTEN (5728),SOX10 (6663),BRAF (673),MAPK1 (5594)
ACH-000001,SK-MEL-1,Skin,-1.95,0.28,-1.82,-1.65,-1.40
ACH-000002,SK-MEL-2,Skin,-2.10,0.31,-1.60,-1.48,-1.25
ACH-000003,A375,Skin,-1.88,0.22,-1.44,-1.52,-1.31
ACH-000004,A549,Lung,-1.90,0.35,0.04,-0.08,-0.12
ACH-000005,H460,Lung,-2.05,0.25,-0.02,-0.15,-0.10
ACH-000006,H1299,Lung,-1.75,0.40,0.01,-0.05,-0.07
"""

# 1. Load and validate the score matrix
df = csa.read_score_csv(csv_bytes)
rankings, diagnostics = csa.analyze_dataset(df)
print(rankings[["gene", "mean_score", "rank", "percentile", "n_valid"]])

# 2. Resolve user gene symbols (including HGNC aliases such as HER2 -> ERBB2)
aliases = csa.load_gene_aliases()
parsed = csa.parse_gene_text("sox10, BRAF (673), MAPK1")
matched, missing, resolved = csa.resolve_gene_list(rankings, parsed["genes"], aliases)

# 3. Compute lineage preference and top/bottom cell lines
lineage_col = csa.find_lineage_column(df)
plot_data = csa.get_lineage_frame(df, matched, rankings, lineage_col)
pref = csa.lineage_preference(plot_data)
models = csa.top_cell_lines(plot_data, matched, top_n=2)

# 4. Run a genome-wide co-dependency scan (Pearson r or Spearman rho)
scan = csa.genomewide_codependency(df, "BRAF", rankings, method="pearson", min_shared=3)
print(scan[["gene", "correlation", "n_shared"]])

# 5. Build CNS-palette Plotly figures and export white-paper vector SVG/PDF/300-DPI PNG
rank_fig = csa.create_rank_plot(rankings, matched, n_cell_lines=len(df), theme="light")
svg_bytes = csa.render_figure_bytes(rank_fig, "svg", width=640, height=480)

# 6. Package a self-contained reproducibility ZIP bundle
zip_bytes = csa.make_analysis_bundle(
    rankings,
    matched,
    plot_data,
    metadata={"source": "synthetic_example", "lineage_preference": pref.to_dict(orient="records")},
    diagnostics=diagnostics,
)
```

---

## Command-line interface (`crispr-score-analyzer`)

Every subcommand works on custom CSV matrices via `--input path/to/matrix.csv` or defaults to the repository's `data/depmap_chronos.npz` when run inside a repository checkout.

### Inspect package and dataset metadata

```bash
crispr-score-analyzer info --input my_scores.csv
```

### Rank genes (pan-cancer or cohort-scoped) with exports

```bash
crispr-score-analyzer rank \
  --input my_scores.csv \
  --genes "KRAS,SOX10,BRAF,HER2" \
  --cohort Skin \
  --output-dir results/skin_rank \
  --plot-format svg \
  --bundle
```

Outputs written to `results/skin_rank/`:
- `selected_results.csv`, `full_rankings.csv`, `plot_data.csv`, `diagnostics.csv`, `lineage_preference.csv`
- `rank_plot.svg` (or `.pdf` / `.png`)
- `analysis_bundle.zip` (when `--bundle` is passed)

For custom CSV files with plain gene headers (without DepMap `SYMBOL (Entrez ID)` suffixes), pass `--gene-columns "MYC,PTEN,SOX10,BRAF"`.

### Co-dependency scan (single gene) or correlation matrix (multi-gene)

```bash
# Single-gene genome-wide Pearson or Spearman scan restricted to Skin cell lines
crispr-score-analyzer codep \
  --input my_scores.csv \
  --query BRAF \
  --cohort Skin \
  --method spearman \
  --top 10 \
  --output results/braf_skin_codep.csv \
  --plot results/braf_skin_codep.svg

# Multi-gene pairwise correlation matrix + heatmap
crispr-score-analyzer codep \
  --input my_scores.csv \
  --genes "TP53,CHEK2,MDM2,MDM4,PPM1D,USP7" \
  --method pearson \
  --output results/tp53_corr.csv \
  --plot results/tp53_corr.svg
```

### Launch the interactive web interface

```bash
crispr-score-analyzer app --port 8501
```

---

## Public module reference

| Module | Key exports | Streamlit required? |
|---|---|---|
| `crispr_score_analyzer.analysis_core` | `read_score_csv`, `read_builtin_matrix`, `analyze_dataset`, `parse_gene_text`, `parse_gene_upload`, `load_gene_aliases`, `resolve_gene_list`, `filter_gene_suggestions`, `lineage_preference`, `top_cell_lines`, `gene_correlation_matrix`, `make_analysis_bundle` | No |
| `crispr_score_analyzer.codependency` | `genomewide_codependency`, `gene_pair_scatter_frame`, `create_codependency_scan_plot`, `create_gene_pair_scatter` | No |
| `crispr_score_analyzer.plotting` | `create_rank_plot`, `create_multilayer_rank_plot`, `create_lineage_boxplot`, `create_preference_plot`, `create_preference_heatmap`, `create_correlation_heatmap`, `PLOT_COLORS`, `PLOT_THEMES`, `preference_colorscale` | No |
| `crispr_score_analyzer.figure_export` | `prepare_export_figure`, `render_figure_bytes`, `export_fingerprint`, `MAX_PNG_PIXELS` | No |
| `crispr_score_analyzer.cli` | `main`, `build_parser` (`crispr-score-analyzer` console script) | Only for `app` subcommand |
