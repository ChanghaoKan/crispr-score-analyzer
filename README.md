# CRISPR Score Analyzer

Explore gene dependency in DepMap CRISPR screens through a browser. Locate genes in a cohort-wide ranking, compare cancer lineages, and export figures with their underlying data and analysis settings.

[Live app](https://crispr-score-analyzer-szbl-denglab.streamlit.app/) · [Report an issue](https://github.com/ChanghaoKan/crispr-score-analyzer/issues) · [Citation](CITATION.cff)

This README describes the source in this checkout. The hosted app may lag behind these changes until its deployment is updated.

## Features

- Instant first-screen gene-ranking render on the built-in DepMap matrix, unified **Examples** menu and **Search** autocomplete, single-control co-dependency search, cross-view gene transfer, shareable URL parameters, exception-only input match feedback, and a quiet unsaved-edit indicator (**Update chart**) whenever edits differ from the committed result.
- Gene dependency rankings within all available cell lines, selected cancer lineages, or optional primary-disease subtypes.
- Cancer-type comparisons with valid sample counts, median or alphabetical ordering, stacked or 2-column grids, single-axis multi-gene comparison, and optional individual cell-line points.
- Two-layer gene-set annotation with configurable colors, reference genes (with alias resolution), and compact or wide canvas options.
- Standalone **Co-dependency** module supporting both **genome-wide single-gene scans** across all ranked genes in the active sidebar scope (including single-lineage cohorts) and **multi-gene correlation matrices** (2–25 genes), each paired with interactive lineage-colored cell-line scatter plots, Pearson/Spearman switching, and figure/CSV exports.
- Deep-dive tabs below ranking, comparison, and annotation results: **Lineage preference** (lollipop with Top/Bottom 6 focus or compact diverging heatmap, cohort highlighting, median difference, `% < -0.5`, and figure export) and **Top cell lines** (most/least dependent cell-line models with CSV export).
- English and Chinese interfaces with instant top-of-sidebar language and light/dark theme controls, compact hero statistics, and scrollable result tables with localized headers.
- Custom score-matrix uploads with explicit gene-column mapping and downloadable validation records.
- PDF, SVG, and PNG figures; result tables; and a ZIP containing rankings, observations, `lineage_preference.csv`, diagnostics, and provenance.

## Quick start

Open the [live app](https://crispr-score-analyzer-szbl-denglab.streamlit.app/) to use the deployed version without a local installation.

For local use, install Python **3.10 or newer**, then run:

```bash
git clone https://github.com/ChanghaoKan/crispr-score-analyzer.git
cd crispr-score-analyzer
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

On Windows, create the environment with `python -m venv .venv` and activate it using `.venv\Scripts\Activate.ps1` in PowerShell. Open the local URL printed by Streamlit, normally `http://localhost:8501`.

The built-in matrix is stored in the repository as `data/depmap_chronos.npz`. Opening the app reads that file from disk. It does not download the score matrix from Hugging Face. The application uses Streamlit, pandas, NumPy, Plotly, Pillow, and Kaleido; dependency constraints are in [requirements.txt](requirements.txt). Streamlit 1.64 or newer is required for on-demand CSV generation. Kaleido is pinned to `0.2.1`, which includes its own rendering engine.

## Usage

1. Switch **Language / 语言** and **Theme** directly at the top of the sidebar, then select **Built-in** or **Custom CSV** under **Data Source**. The upload control appears for custom data; map plain-symbol gene columns explicitly when using a custom matrix.
2. Choose **Cancer types** (and optional **Disease / subtype** when present in the dataset) under **Analysis scope**. An empty selection uses all available lineages. Reference genes, labels, point size, and canvas aspect ratio are under **Figure settings** next to a ranking chart, so those display changes do not reload the dataset.
3. Choose **Gene ranking**, **Cancer-type comparison**, **Gene-set annotation**, or **Co-dependency**. On first load of the built-in matrix, **Gene ranking** immediately renders an iconic cancer dependency panel; edit the gene list, toggle **Upload file**, or pick curated sets from **Examples** and **Search**, then click **Run analysis** (or **Update chart** after editing).
4. Review the run's genes, cohort, and UTC time in the status strip above the plot, click **Lineages** / **Ranking** / **Co-dependency** to move genes seamlessly across views, or click **Share** to write `?view=...&genes=...&cohort=...` into the browser address bar. A compact results table appears right below the primary plot beside **Export figure** and **Download results** controls, followed by deep-dive tabs.

Gene-list and per-view analysis edits remain drafts until Run/Update is clicked. Input matching follows exception-only feedback: when every symbol matches cleanly, the canvas stays uncluttered; previous/alias symbols, duplicates, and unmatched names are reported immediately before and after submission. Alias resolution uses HGNC previous symbols and aliases that point to exactly one approved symbol in the bundled matrix (`data/gene_aliases.json`). When an edited draft differs from the committed result, the status indicator turns amber (`Unsaved edit`) and the primary button morphs to **Update chart**. Each view retains its last successful result within the current session; an unmatched submission leaves that result available. Running one view does not submit another view's inputs. Changing the dataset, gene-column mapping, or sidebar cohort clears previous results and requires another Run. Presentation settings, such as reference genes and label visibility, update the saved result's display directly and persist across tab switches. The background ranking curve is reused when only those display settings change.

The **Lineage preference** tab compares each cancer type against each gene's full-dataset average (independent of the sidebar cohort filter, with active sidebar lineages outlined/starred). Switch between per-gene **Lollipop** charts (with a **Top & bottom 6** vs **All types** toggle) and a compact diverging **Heatmap** (selected automatically when more than six genes are shown). Hover displays mean difference, median difference, fraction of cell lines below `-0.5`, and valid sample counts. Types with fewer than five cell lines are left off the figure, while all rows are saved in `lineage_preference.csv` and `metadata.json` inside the analysis ZIP. All preference values are descriptive and are not a test of selectivity or significance.

The **Top cell lines** tab lists the most dependent (lowest score) and least dependent (highest score) cell lines per gene in the active scope for wet-lab model selection.

### Gene ranking

Locate matched genes within the selected cohort's mean-score ranking. Hover over a point for the gene, mean score, rank, rank percentile, and valid sample count. Reference A and Reference B are configurable visual anchors (supporting HGNC alias resolution and optional clearing).

Ranking charts default to a centered 4:3 layout (640 × 480 px) and also support a wide 16:10 canvas (840 × 500 px), with transparent gene labels and leader lines that keep the underlying score curve visible. Zoom, pan, and reset-axes controls remain available in the chart toolbar. Charts shrink to fit narrower screens. Gene-set annotation uses the same presentation, and figure exports default to these dimensions.

### Cancer-type comparison

Choose genes, then open **Comparison options** to optionally restrict the displayed lineages within the sidebar cohort, choose median-score or alphabetical ordering, switch between **By cancer type** and **Compare genes together** (single shared axis), or lay out multiple genes in **1 column** or **2 columns**. Each group shows its valid `n`. Enable **Show all cell lines** to display every observation; point hover includes the cell-line name and identifier when available.

If more than eight genes match, choose up to eight explicitly before generating the plot. The complete submitted input remains recorded in the analysis metadata. The accompanying gene ranking uses the sidebar cohort; group summaries and plotted observations use the boxplot's additional lineage selection.

### Gene-set annotation

Enter a background gene set and a highlight set, choose their colors under **Annotation colors**, and click **Run analysis**. Both lists receive matching feedback. Highlight genes take precedence when the sets overlap. This view uses the same cohort ranking as the gene-ranking view; it does not run a gene-set enrichment test.

### Co-dependency

Explore pairwise CRISPR dependency correlations across the cell lines in the active sidebar scope using either **Pearson r** or **Spearman ρ**:

- **Genome-wide scan (1 gene):** Correlates a single query gene against every other ranked gene (~18,433 genes on the bundled matrix) using exact pairwise-complete observations. Because the scan respects the sidebar cancer-type and subtype filters, you can discover lineage-restricted co-dependencies within a specific cancer type as well as pan-cancer partners. Displays a diverging lollipop of top positive and negative partners, an interactive lineage-colored **cell-line scatter plot** with OLS trendline for any selected partner, top-50 positive/negative tables, one-click transfer of top partners to a correlation matrix or other views, and full-genome CSV / vector figure exports.
- **Multi-gene matrix (2–25 genes):** Computes a pairwise correlation heatmap for a custom gene panel alongside an interactive Gene X vs Gene Y cell-line scatter plot, with figure and CSV exports.

All co-dependency correlations are descriptive and may reflect shared pathways, lineage specificity, or common-essentiality structure rather than direct genetic interaction.

## Data and input formats

### Built-in data

The default is **DepMap Public 26Q3 Chronos Gene Effect**, bundled with the app as `data/depmap_chronos.npz`. `data/depmap_chronos.provenance.json` records the release, the source-file checksums, and the matrix checksum.

| Setting | Value |
|---|---|
| Served file | `data/depmap_chronos.npz` (88 MB) |
| Matrix | 1,361 models × 18,434 genes |
| Scores | float32 Chronos gene effect; missing observations are NaN |
| Metadata columns, in order | `depmap_id`, `cell_line_display_name`, `lineage` |
| Gene columns | `SYMBOL (Entrez ID)`, taken from `CRISPRGeneEffect.csv` |
| Lineage | `OncotreeLineage` from the same release's `Model.csv` |

The three metadata columns are not gene scores. `lineage` is the cancer-type scope. Uploaded files do not inherit these column names. The number of ranked genes can vary after column validation or cohort-specific missingness.

The app reads the bundled matrix from disk. It does not download a score matrix when a page opens. On this machine, reading the file and ranking every gene takes about 1.2 seconds. The first visit after Community Cloud has slept still waits for the platform to start the container.

Rebuild the bundle from an official portal download of the same release:

```bash
python scripts/build_builtin_matrix.py \
  --effect data/raw/CRISPRGeneEffect.csv \
  --model data/raw/Model.csv \
  --release "DepMap Public 26Q3"
```

`data/raw/` is gitignored. DepMap's terms cover the bundled matrix; the MIT license covers the application code. Cite the release used for a study: DepMap, Broad (2026). DepMap Public 26Q3. https://depmap.org/portal/

### Custom score matrices

Upload a UTF-8 CSV with one row per cell line and numeric gene scores in columns, for example:

```csv
ModelID,cell_line_display_name,lineage,MYC (4609),PTEN (5728)
ACH-000001,Example A,Liver,-1.95,0.12
ACH-000002,Example B,Lung,-2.13,0.08
```

- Automatic detection requires `SYMBOL (numeric Entrez ID)`, such as `MYC (4609)`. Numeric metadata is not automatically treated as a gene.
- For plain-symbol headers such as `MYC`, open **Gene column mapping** in the sidebar, enable **Select columns manually**, and select the score columns. Other descriptive headers must first be renamed to valid symbols.
- Include a cell-line identifier such as `ModelID`, `depmap_id`, or `cell_line`. Include `cell_line_display_name` or `CellLineName` for readable point labels. Without an identifier, row indices are used.
- Include `lineage`, `OncotreeLineage`, or another recognized lineage column for cohort filtering and boxplots. Optional `OncotreePrimaryDisease` / `primary_disease` columns also enable subtype filtering. Missing labels become `Unknown`; ranking remains available when lineage metadata is absent.
- All columns sharing an ambiguous gene symbol, ignoring case, are excluded. Repeated recognized metadata headers are rejected. The analysis ZIP includes column inclusion/exclusion reasons in `diagnostics.csv`.
- Ensure each cell line occurs once. The app does not deduplicate repeated observations, so repeated rows would contribute repeatedly to means and group counts.

### Data validation

The app validates columns and scores before ranking, and reports input errors where they occur. Cell-line identifiers, names and lineage metadata describe samples and are not part of the gene ranking. Full validation records are available as `diagnostics.csv` inside **Prepare analysis ZIP** under **Details & reproducibility**. If no gene columns can be analyzed, the error message provides a separate **Download validation report** button.

The `position` field in the full CSV is the column's **zero-based position in the original table**, not its gene rank. For the bundled matrix, the three metadata columns occupy positions 0–2, so the first gene starts at position 3.

### Gene lists

Gene-list input is separate from score-matrix upload. Use gene symbols such as `MYC`; pasted DepMap-style headers such as `MYC (4609)` are also automatically stripped of their numeric Entrez suffix.

Paste symbols separated by whitespace, newlines, tabs, commas, or semicolons; Chinese commas and semicolons are also accepted. Matching is case-insensitive. Duplicate inputs are removed while preserving their first occurrence, and unmatched genes are shown for correction.

The ranking and boxplot views also accept UTF-8 `.txt`, `.csv`, and `.tsv` lists. CSV/TSV uploads read the first column. A recognized header such as `gene`, `gene_symbol`, or `symbol` is optional; a headerless file keeps its first gene. Gene-set annotation uses the two text inputs.

## Analysis and interpretation

Each gene's score is the arithmetic mean of its **finite numeric observations** in the selected cohort. Blank, nonnumeric, and infinite values are excluded from the mean and counted as missing. Diagnostics also report nonnumeric and infinite values separately.

Genes are retained regardless of variance or mean-score range if they have at least one finite value. Genes without finite values cannot receive a rank. Check `n_valid`, especially for small cohorts or incomplete custom matrices.

- **Rank:** ascending mean score, starting at 1. Tied means retain the dataset's column order.
- **Rank percentile:** `100 × rank / number of ranked genes`. Lower values indicate stronger mean dependency within that ranking.
- **Reference cutoff:** the `−0.5` line is a descriptive mean-score screen, not a classification of common essential genes.
- **Missingness:** `n_missing` counts all excluded observations; `missing_fraction = n_missing / cohort_rows`.

Ranks depend on both the selected cell lines and the analyzable gene columns. A custom subset does not produce a genome-wide ranking. Boxplot differences and dependency ranks are descriptive; they do not establish cancer-type selectivity, statistical significance, mechanism, or a therapeutic window.

## Exports

Open **Export figure** beside the results heading (or inside the Lineage preference / Co-dependency tabs), click **PDF**, **PNG**, or **SVG**, then download the resulting file. Use **Figure size** within this panel to change dimensions. Exports use a white page background, white marker halos, and transparent text labels regardless of the interface theme. Cached bytes are tied to the figure content, format, and dimensions, so a changed result does not reuse an earlier figure.

| Format | Behavior |
|---|---|
| PDF / SVG | Full vector export for all supported analysis views; interactive WebGL ranking points are converted to vector traces. |
| PNG | Raster export at 3× the selected width and height, with embedded 300-DPI metadata. |

For example, an export size of 1,000 × 600 produces a 3,000 × 1,800-pixel PNG. DPI metadata does not increase pixel resolution. PNG output is capped at **24 million pixels**; reduce dimensions or choose PDF/SVG for larger figures.

The visible results table shows gene symbols, mean scores, ranks, rank percentiles, and valid sample counts. Browser labels follow the selected language; means display three decimal places and percentiles display a percent sign. **Download results** exports CSV with canonical field names and full numeric precision, independent of these display formats.

Open **Details & reproducibility** for missing-value counts, cancer-type summaries, plot-data CSV, and analysis metadata. Click **Prepare analysis ZIP** to generate the complete bundle, then download it. Bundles are generated only on request; changes to the analysis or recorded display/export settings require a new bundle. Pending input edits do not change exported results until Run is clicked:

| ZIP member | Contents |
|---|---|
| `selected_results.csv` | Selected matched genes with raw column names, symbols, means, ranks, percentiles, and valid/missing counts. |
| `full_rankings.csv` | Complete ranking for the current dataset, column mapping, and sidebar cohort. |
| `plot_data.csv` | Finite observations with `gene`, `lineage`, `crispr_score`, `cell_line`, and `cell_line_name`; boxplot lineage restrictions are applied. |
| `lineage_preference.csv` | Per-lineage mean/median differences, all-lineage averages, `% < -0.5` dependent fraction, and valid counts (when lineage metadata is available). |
| `metadata.json` | Data source, release, pinned revision, file SHA-256, cohort, column mapping, input matching, parameters, timestamp, software versions, ranking method, and display/export settings. |
| `diagnostics.csv` | Column inclusion/exclusion status and reason, with valid, missing, invalid, and infinite counts for recognized candidate columns. |

For gene-ranking and gene-set annotation views, observation data covers the selected genes; `full_rankings.csv` also records the background ranking. The ZIP does not contain figure files or a copy of the entire original score matrix. Keep the source matrix and downloaded figures alongside the bundle when archiving an analysis.

## Tests and development

Install application dependencies as above, then install the test runner and run the full suite from the repository root:

```bash
python -m pip install pytest
python -m pytest tests -q
```

The suite covers parsing, ranking and missingness, draft and submitted-result state, cohort changes, localized table formatting without loss of CSV precision, on-demand bundle generation, export cache invalidation, citation wording, and actual SVG/PDF/PNG rendering. Tests use small deterministic fixtures. They do not fetch DepMap data. Rendering tests need permission to launch Kaleido's bundled renderer.

For a check that excludes the actual image-rendering test:

```bash
python -m pytest tests -q -k "not actual_svg_is_vector_pdf_is_valid_and_png_has_dpi"
```

GitHub Actions runs the full suite on Python 3.12. Automated checks do not replace visual review of the deployed app or validation of a custom dataset.

## Citation

Cite CRISPR Score Analyzer 1.1.0, deposited 2026-11-21: [10.5281/zenodo.24424595](https://doi.org/10.5281/zenodo.24424595). The concept DOI [10.5281/zenodo.19607602](https://doi.org/10.5281/zenodo.19607602) resolves to that version. The Zenodo record and the code in this repository are both MIT. Record DepMap Public 26Q3 with an analysis from the current app. The source archive on Zenodo does not include the bundled matrix.

The 2026-04-16 snapshot remains version 1.0.0: [10.5281/zenodo.19607603](https://doi.org/10.5281/zenodo.19607603). That record's Zenodo license is CC-BY-4.0.

Also cite the underlying DepMap release and relevant methodology. The bundled matrix is DepMap Public 26Q3. The repository records this foundational reference:

> DepMap, Broad (2026). DepMap Public 26Q3. https://depmap.org/portal/

> Tsherniak, A., et al. (2017). Defining a Cancer Dependency Map. *Cell*, 170, 564–576. [doi:10.1016/j.cell.2017.06.010](https://doi.org/10.1016/j.cell.2017.06.010).

## Contributing and contact

Report bugs and request features through [GitHub Issues](https://github.com/ChanghaoKan/crispr-score-analyzer/issues). Include reproduction steps, relevant package versions, and a small synthetic example when reporting a data problem. Pull requests should describe the behavior changed and the checks performed.

Maintainer: [Changhao Kan](https://github.com/ChanghaoKan).

## License and acknowledgements

The application code is licensed under the [MIT License](License). DepMap and other source datasets retain their own terms; consult the [DepMap Portal](https://depmap.org/portal/) before redistributing data.

Thanks to the DepMap Consortium for the dependency data.

### Development assistance

[Claude (Anthropic)](https://www.anthropic.com/claude) and [Codex (OpenAI)](https://openai.com/codex/) provided AI development assistance. Codex assisted with code and interface improvements, testing, and documentation. The maintainer is responsible for reviewing these contributions and validating the software's scientific use.
