# CRISPR Score Analyzer

Explore gene dependency in DepMap CRISPR screens through a browser. Locate genes in a cohort-wide ranking, compare cancer lineages, and export figures with their underlying data and analysis settings.

[Live app](https://crispr-score-analyzer-szbl-denglab.streamlit.app/) · [Dataset](https://huggingface.co/datasets/ChanghaoKan/crispr-depmap) · [Report an issue](https://github.com/ChanghaoKan/crispr-score-analyzer/issues) · [Citation](CITATION.cff)

This README describes the source in this checkout. The hosted app may lag behind these changes until its deployment is updated.

## Features

- Explicit **Run analysis / Run · 生成图表** controls, example and clear actions, input matching feedback, and a record of the last submitted analysis in each view.
- Gene dependency rankings within all available cell lines or selected cancer lineages.
- Cancer-type comparisons with valid sample counts, median or alphabetical ordering, and optional individual cell-line points.
- Two-layer gene-set annotation with configurable colors and reference genes.
- English and Chinese interfaces, light and dark themes, a card-based layout with a consistent product title, compact navigation, and readable result tables with localized headers.
- Custom score-matrix uploads with explicit gene-column mapping and column diagnostics.
- PDF, SVG, and PNG figures; result tables; and a ZIP containing rankings, observations, and provenance.

The experimental gene–drug correlation code is retained in the repository but disabled in the public interface.

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

The initial built-in dataset download requires an internet connection. The application uses Streamlit, pandas, NumPy, Plotly, Pillow, and Kaleido; dependency constraints are in [requirements.txt](requirements.txt). Kaleido is pinned to `0.2.1`, which includes its own rendering engine.

## Usage

1. Select **Built-in** or **Custom CSV** under **Data Source** in the sidebar. The upload control appears for custom data; map plain-symbol gene columns explicitly when using a custom matrix.
2. Choose **Cancer types** under **Analysis scope**. An empty selection uses all available lineages. Figure settings follow the cohort controls; language and theme are under **Preferences**.
3. Choose **Gene ranking**, **Cancer-type comparison**, or **Gene-set annotation**, enter gene symbols, and click **Run analysis**. **Load example** and **Clear** help edit the input.
4. Review the run's genes, cohort, and UTC time above the plot. A compact results table appears below the plot, beside **Export figure** and **Download results** controls.

Gene-list and per-view analysis edits are drafts until Run is clicked. After leaving an edited input, the app indicates when it differs from the displayed analysis. Each view retains its last successful result within the current session; an unmatched submission leaves that result available. Running one view does not submit another view's inputs. Changing the dataset, gene-column mapping, or sidebar cohort clears previous results and requires another Run. Presentation settings, such as reference genes and label visibility, update the saved result's display directly.

### Gene ranking

Locate matched genes within the selected cohort's mean-score ranking. Hover over a point for the gene, mean score, rank, rank percentile, and valid sample count. Reference A and Reference B are configurable visual anchors.

Ranking charts use a compact, centered 4:3 layout (640 × 480 px), with transparent gene labels and leader lines that keep the underlying score curve visible. Charts shrink to fit narrower screens. Gene-set annotation uses the same presentation, and figure exports default to these dimensions.

### Cancer-type comparison

Choose genes, then open **Comparison options** to optionally restrict the displayed lineages within the sidebar cohort and select median-score or alphabetical ordering. Each group shows its valid `n`. Enable **Show all cell lines** to display every observation; point hover includes the cell-line name and identifier when available.

If more than eight genes match, choose up to eight explicitly before generating the plot. The complete submitted input remains recorded in the analysis metadata. The accompanying gene ranking uses the sidebar cohort; group summaries and plotted observations use the boxplot's additional lineage selection.

### Gene-set annotation

Enter a background gene set and a highlight set, choose their colors under **Annotation colors**, and click **Run analysis**. Both lists receive matching feedback. Highlight genes take precedence when the sets overlap. This view uses the same cohort ranking as the gene-ranking view; it does not run a gene-set enrichment test.

## Data and input formats

### Built-in data

The default is the curated **DepMap Public 25Q3 Chronos Gene Effect** matrix hosted in [ChanghaoKan/crispr-depmap](https://huggingface.co/datasets/ChanghaoKan/crispr-depmap).

| Setting | Value |
|---|---|
| File | `CRISPR_(DepMap_Public_25Q3+Score,_Chronos)_subsetted.csv` |
| Pinned Hugging Face revision | `8400bf566411ed1df4e0784fdb9e97fdbaa371fd` |
| Metadata columns, in order | `depmap_id`, `cell_line_display_name`, `lineage_1`, `lineage_2`, `lineage_3`, `lineage_6`, `lineage_4` |
| Gene columns | Plain gene symbols following these seven metadata columns |

The app explicitly maps this known schema only for the built-in file. It uses `lineage_1` for its cancer-type scope. Uploaded files do not inherit the built-in column assumptions. The number of ranked genes can vary after column validation or cohort-specific missingness.

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
- Include `lineage`, `OncotreeLineage`, or another recognized lineage column for cohort filtering and boxplots. Missing labels become `Unknown`; ranking remains available when lineage metadata is absent.
- All columns sharing an ambiguous gene symbol, ignoring case, are excluded. Repeated recognized metadata headers are rejected. Inspect column diagnostics for each exclusion reason.
- Ensure each cell line occurs once. The app does not deduplicate repeated observations, so repeated rows would contribute repeatedly to means and group counts.

### Gene lists

Gene-list input is separate from score-matrix upload. Use gene symbols, such as `MYC`, rather than matrix headers such as `MYC (4609)`.

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

Open **Export figure** beside the results heading, click **PDF**, **PNG**, or **SVG**, then download the resulting file. Use **Figure size** within this panel to change dimensions. Exports use a white page background and transparent text labels regardless of the interface theme. Cached bytes are tied to the figure content, format, and dimensions, so a changed result does not reuse an earlier figure.

| Format | Behavior |
|---|---|
| PDF / SVG | Full vector export for the three supported analysis views; interactive WebGL ranking points are converted to vector traces. |
| PNG | Raster export at 3× the selected width and height, with embedded 300-DPI metadata. |

For example, an export size of 1,000 × 600 produces a 3,000 × 1,800-pixel PNG. DPI metadata does not increase pixel resolution. PNG output is capped at **24 million pixels**; reduce dimensions or choose PDF/SVG for larger figures.

The visible results table shows gene symbols, mean scores, ranks, rank percentiles, and valid sample counts. Browser labels follow the selected language; means display three decimal places and percentiles display a percent sign. **Download results** exports CSV with canonical field names and full numeric precision, independent of these display formats.

Open **Details & reproducibility** for missing-value counts, cancer-type summaries, plot-data CSV, and analysis metadata. Click **Prepare analysis ZIP** to generate the complete bundle, then download it. Bundles are generated only on request; changes to the analysis or recorded display/export settings require a new bundle. Pending input edits do not change exported results until Run is clicked:

| ZIP member | Contents |
|---|---|
| `selected_results.csv` | Selected matched genes with raw column names, symbols, means, ranks, percentiles, and valid/missing counts. |
| `full_rankings.csv` | Complete ranking for the current dataset, column mapping, and sidebar cohort. |
| `plot_data.csv` | Finite observations with `gene`, `lineage`, `crispr_score`, `cell_line`, and `cell_line_name`; boxplot lineage restrictions are applied. |
| `metadata.json` | Data source, release, pinned revision, file SHA-256, cohort, column mapping, input matching, parameters, timestamp, software versions, ranking method, and display/export settings. |
| `diagnostics.csv` | Column inclusion/exclusion status and reason, with valid, missing, invalid, and infinite counts for recognized candidate columns. |

For gene-ranking and gene-set annotation views, observation data covers the selected genes; `full_rankings.csv` also records the background ranking. The ZIP does not contain figure files or a copy of the entire original score matrix. Keep the source matrix and downloaded figures alongside the bundle when archiving an analysis.

## Tests and development

Install application dependencies as above, then install the test runner and run the full suite from the repository root:

```bash
python -m pip install pytest
python -m pytest tests -q
```

The suite covers parsing, ranking and missingness, draft and submitted-result state, cohort changes, localized table formatting without loss of CSV precision, on-demand bundle generation, export cache invalidation, and actual SVG/PDF/PNG rendering. Tests use small deterministic fixtures and mock the public dataset download; they do not fetch DepMap data. Rendering tests need permission to launch Kaleido's bundled renderer.

For a check that excludes the actual image-rendering test:

```bash
python -m pytest tests -q -k "not actual_svg_is_vector_pdf_is_valid_and_png_has_dpi"
```

GitHub Actions runs the full suite on Python 3.12. Automated checks do not replace visual review of the deployed app or validation of a custom dataset.

## Citation

See [CITATION.cff](CITATION.cff) for the recorded software citation. It identifies the project concept DOI as [10.5281/zenodo.19607602](https://doi.org/10.5281/zenodo.19607602). Cite the archived software version used for a study, and record the commit when using changes that are not part of that release.

Also cite the underlying DepMap release and relevant methodology. The repository records this foundational reference:

> Tsherniak, A., et al. (2017). Defining a Cancer Dependency Map. *Cell*, 170, 564–576. [doi:10.1016/j.cell.2017.06.010](https://doi.org/10.1016/j.cell.2017.06.010).

## Contributing and contact

Report bugs and request features through [GitHub Issues](https://github.com/ChanghaoKan/crispr-score-analyzer/issues). Include reproduction steps, relevant package versions, and a small synthetic example when reporting a data problem. Pull requests should describe the behavior changed and the checks performed.

Maintainer: [Changhao Kan](https://github.com/ChanghaoKan), Deng Lab, Shenzhen Bay Laboratory.

## License and acknowledgements

The application code is licensed under the [MIT License](License). DepMap and other source datasets retain their own terms; consult the [DepMap Portal](https://depmap.org/portal/) before redistributing data.

Thanks to the DepMap Consortium for the dependency data and Deng Lab for institutional support.
