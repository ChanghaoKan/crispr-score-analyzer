# CRISPR Score Analyzer: an interactive, cohort-scoped browser platform for DepMap gene essentiality and co-dependency analysis

**Changhao Kan¹,\***

¹ Shenzhen Bay Laboratory, Shenzhen, China  
\* Correspondence: `kanchanghao@szbl.ac.cn`  
ORCID: [0009-0000-4014-5923](https://orcid.org/0009-0000-4014-5923)

**Category:** Bioinformatics (Application Note)

---

## Abstract

**Summary:** Genome-scale CRISPR-Cas9 loss-of-function screens in the Cancer Dependency Map (DepMap) are widely used to prioritize candidate cancer targets, evaluate lineage-restricted vulnerabilities, and inspect functional gene relationships. Moving beyond single-gene portal lookups—such as placing multi-gene candidate panels on a cohort-wide essentiality curve, separating pan-cancer essentiality from tissue-specific dependency, pairing sensitive and resistant cell lines for experimental validation, running single-lineage genome-wide co-dependency scans, and exporting vector figures with analysis provenance—often requires ad-hoc scripting. We present **CRISPR Score Analyzer (v1.2.0)**, an open-source Python/Streamlit browser platform bundling the DepMap Public 26Q3 Chronos gene effect matrix (1,361 cancer cell lines × 18,434 genes across 30 Oncotree lineages) and 40,402 unambiguous HGNC alias mappings, while also supporting validated custom CSV uploads. Across four workflows—*Gene ranking*, *Cancer-type comparison*, *Gene-set annotation*, and *Co-dependency*—users can interrogate pan-cancer or lineage-restricted cohorts, run chunked vectorized Pearson or Spearman scans in seconds using NumPy and pandas without requiring SciPy, and export white-paper PDF/SVG/300-DPI PNG figures alongside SHA-256 ZIP provenance bundles. Benchmark use cases demonstrate recapitulating lineage-selective dependencies (`SOX10` in skin, `ESR1` in ER-positive breast models), pan-cancer `MDM2`–`TP53` antagonism (`r = −0.702`), and strengthened intra-lineage `BRAF`–MAPK co-dependencies within skin cancer cell lines (`MAPK1` `r = +0.798`).

**Availability and Implementation:** Source code for v1.2.0 is freely available under the MIT license at https://github.com/ChanghaoKan/crispr-score-analyzer (release tag `v1.2.0`: https://github.com/ChanghaoKan/crispr-score-analyzer/releases/tag/v1.2.0; Zenodo concept archive DOI: `10.5281/zenodo.19607602`). A hosted web instance is accessible at https://crispr-score-analyzer-szbl-denglab.streamlit.app/. Because MIT release source archives exclude the bundled DepMap matrix under DepMap portal terms, exact reproduction of the analyses and figures requires cloning tag `v1.2.0`, verifying `data/depmap_chronos.npz` against the SHA-256 checksum in `data/depmap_chronos.provenance.json` (`ddf274da2aa5f89cf630bc8150dcd3a2c722466bd18ff6fb6d18fe161f8a0c67`, or rebuilding via `scripts/build_builtin_matrix.py`), and running `python preprint/generate_figures.py`.

**Keywords:** DepMap; CRISPR-Cas9 screen; Chronos gene effect; gene essentiality; co-dependency; cancer cell lines; reproducibility

---

## 1. Introduction

Pooled genome-scale CRISPR-Cas9 loss-of-function screens across human cancer cell lines—led by the Broad Institute Cancer Dependency Map (DepMap) (Tsherniak *et al.*, 2017; Meyers *et al.*, 2017; Dempster *et al.*, 2021)—provide a standard reference for oncology target prioritization, pathway mapping, and experimental model selection. Using copy-number- and kinetics-corrected Chronos gene effect scores (Dempster *et al.*, 2021), researchers routinely assess whether candidate genes from transcriptomic, proteomic, or genetic screens represent broad pan-cancer dependencies or context-specific vulnerabilities.

Although the DepMap Portal (https://depmap.org/portal/) and companion tools such as shinyDepMap (Shimada *et al.*, 2021) and co-essentiality atlases (Wainberg *et al.*, 2021) provide extensive access to dependency resources, everyday workflows still encounter four recurring bottlenecks:

1. **Multi-gene batch context:** Positioning a candidate panel of 5–50 genes simultaneously along a cohort-wide essentiality curve—or overlaying a prioritized subset onto a pathway background—typically requires repetitive single-gene lookups or ad-hoc R/Python scripting.
2. **Bridging lineage summaries to wet-lab model selection:** Pan-cancer averages frequently dilute lineage- or subtype-restricted vulnerabilities, while selecting concrete sensitive and resistant cell lines within a chosen lineage scope requires cross-referencing multiple tables.
3. **Direct cohort-scoped genome-wide co-dependency scanning:** Default portal co-dependency tables report pan-cancer correlations where inter-lineage expression variance can confound pathway co-essentiality; running a live single-lineage or subtype-restricted genome-wide scan in one browser step typically requires multi-step custom explorer setup.
4. **Reproducible publication export:** Producing vector-editable figures (PDF/SVG) paired with exact numeric tables, data checksums, input-matching diagnostics, and software environment metadata remains fragmented across tools.

To streamline these exploratory workflows without requiring programming expertise, we developed **CRISPR Score Analyzer**, an open-source interactive web application designed around explicit input-state isolation, transparent descriptive statistics, vector figure export, and self-contained provenance archiving.

---

## 2. Implementation and Core Capabilities

### 2.1 Bundled Reference Data, Custom Uploads, and Symbol Resolution

CRISPR Score Analyzer is implemented in Python (Streamlit, pandas, NumPy, Plotly, Kaleido, and Pillow) and provides English and Chinese interface localization. To ensure deterministic startup without remote matrix downloads, the repository bundles a compact `float32` binary matrix (`data/depmap_chronos.npz`, 88 MB; maximum absolute `float32` rounding deviation `1.19 × 10⁻⁷` across the first 100 genes) of the **DepMap Public 26Q3 Chronos Gene Effect** dataset comprising **1,361 cancer cell lines** and **18,434 genes** across **30 Oncotree lineages**, paired with SHA-256 and MD5 provenance metadata (`data/depmap_chronos.provenance.json`). The bundled matrix stores `OncotreeLineage` metadata; optional subtype filtering (`OncotreePrimaryDisease` / `primary_disease`) activates automatically when present in custom uploaded tables.

Researchers can upload custom UTF-8 CRISPR score CSV matrices. To guard against silent misinterpretation of numeric metadata columns, automatic gene detection strictly requires DepMap-style `SYMBOL (Entrez ID)` headers, whereas plain-symbol tables require explicit manual score-column mapping in the sidebar. Ambiguous duplicate gene symbols and repeated metadata columns are excluded automatically and logged in a downloadable `diagnostics.csv` validation report.

Gene-list inputs accept arbitrary ASCII or full-width delimiters (newlines, spaces, tabs, commas, and semicolons), `.txt/.csv/.tsv` file uploads, and pasted DepMap headers (e.g., `MYC (4609)`). An offline HGNC symbol map (`data/gene_aliases.json`, **40,402** unambiguous previous-symbol and alias mappings) (Seal *et al.*, 2023) resolves legacy names (such as `HER2 → ERBB2`) automatically. Symbol matching reports resolved aliases, deduplicated entries, and unmatched queries explicitly before and after submission while leaving exact matches unflagged.

### 2.2 Modular Analytical Workflows

The interface organizes analyses into four top-level modules governed by a shared sidebar cohort scope (**Cancer types** and optional **Disease / subtype** when present in custom datasets; **Table 1**):

1. **Gene ranking:** Ranks all valid genes within the active cohort by ascending mean finite Chronos gene effect (`rank`, `percentile = 100 × rank / N`, `n_valid`, and missingness counts). Interactive scatter curves render background genes via WebGL (`Scattergl`) for responsiveness while converting all traces to vector paths during PDF/SVG export. Two configurable reference anchors (default `MYC` and `PTEN`, supporting alias resolution) and a dashed `−0.5` screening reference line provide visual orientation, accompanied by collision-avoided transparent gene labels with leader lines.
2. **Cancer-type comparison:** Generates horizontal per-lineage boxplots displaying valid sample counts (`n`) per group, sorted by median score or alphabetically, with optional 2-column grids, single-axis multi-gene comparison, and individual cell-line point overlays annotated with display name and `ModelID`. Submissions matching more than eight genes prompt explicit subset selection while preserving the complete submitted input in metadata.
3. **Gene-set annotation:** Plots a background gene set beneath a prioritized highlight layer on the cohort ranking curve with configurable colors, enabling pathway-level positioning without conflating visualization with formal enrichment testing.
4. **Co-dependency:** Supports two complementary modes with instant switching between **Pearson *r*** and **Spearman *ρ***:
   - *Genome-wide scan (1 gene):* Correlates a single query gene (selected via a searchable dropdown accepting approved symbols and HGNC aliases) against every other ranked gene in the active sidebar scope using exact pairwise-complete finite observations in `float64` precision across 2,048-gene column chunks. Pairs must meet a minimum shared observation threshold (`min_shared = min(10, max(3, n_valid_query))`, i.e., `n_shared ≥ 10` for cohorts of at least 10 lines) and non-zero variance; Spearman *ρ* resolves ties via average ranks using NumPy and pandas without requiring SciPy. Outputs include a diverging top-partner lollipop chart, top-50 positive and negative tables (displaying `n_shared` alongside correlation), one-click transfer of top partners to a correlation matrix, and full-genome CSV export.
   - *Multi-gene matrix (2–25 genes):* Computes an interactive pairwise correlation heatmap alongside an integrated **lineage-colored cell-line scatter plot** with ordinary least-squares trendline and both Pearson *r* and Spearman *ρ* statistics for any chosen gene pair.

**Table 1. Workflow comparison across CRISPR Score Analyzer (v1.2.0) and existing DepMap web interfaces (capabilities verified November 2026).**

| Workflow dimension | DepMap Portal | shinyDepMap | CRISPR Score Analyzer (v1.2.0) |
|---|---|---|---|
| Multi-gene placement on cohort ranking curve | Single-gene default view | Efficacy vs selectivity scatter | Batch ranking + two-layer gene-set annotation |
| Cohort-scoped genome-wide co-dependency scan | Pan-cancer default; subset scans via multi-step Custom Analysis | Pan-cancer co-dependency clusters | Single-step live scoped scan (Pearson *r* & Spearman *ρ*) |
| Pairwise multi-gene matrix & lineage scatter | Custom Analysis / Data Explorer small multiples | Cluster view | 2–25 gene heatmap + integrated lineage scatter |
| Lineage preference vs all-lineage mean | Per-gene lineage & subtype plots | Not provided | Diverging Top/Bottom-6 lollipop & multi-gene heatmap |
| Paired sensitive / resistant model lookup | Table filtering | Top sensitive cell lines | Paired most/least dependent tables + CSV export |
| Custom CRISPR score matrix upload & diagnostics | Single-feature custom vector (not full score matrix) | Not supported | Full UTF-8 CSV matrix + strict schema + `diagnostics.csv` |
| Offline batch HGNC alias resolution | Online single-gene search | Not provided | 40,402 unambiguous offline mappings (batch & single) |
| Vector export & analysis provenance archive | Separate chart and CSV downloads | Separate PNG and CSV | White-paper PDF/SVG/300-DPI PNG + SHA-256 ZIP bundle |

### 2.3 Deep-Dive Views, State Integrity, Performance, and Verification

Beneath ranking, comparison, and annotation results, two contextual deep-dive tabs provide immediate follow-up without switching modules: (i) **Lineage preference**, which compares every cancer lineage (with `n_valid ≥ 5`) against each gene's full-dataset mean (independent of the sidebar filter, outlining active sidebar lineages) in either Top/Bottom-6 lollipop or diverging heatmap layouts alongside median difference and the fraction of lines scoring `< −0.5`; and (ii) **Top cell lines**, listing the most dependent and least dependent models in the active scope for experimental selection.

To prevent mismatch between displayed outputs and edited inputs, gene-list modifications remain uncommitted drafts until explicitly re-executed, and changes to the dataset, column mapping, or sidebar cohort invalidate prior results and exports. All figure exports enforce white paper backgrounds, marker halos, and transparent label boxes, and cache exported bytes by exact figure fingerprint and dimensions (with PNG output capped at 24 million pixels). An on-demand **Complete analysis ZIP** packages `selected_results.csv`, `full_rankings.csv`, `plot_data.csv`, `lineage_preference.csv`, `diagnostics.csv`, and `metadata.json` (recording data SHA-256, release version, Zenodo DOI, cohort parameters, input matching, and Python package versions).

On a 2-vCPU Linux container (`Python 3.12.16`, `NumPy 2.5.5`, `pandas 2.3.3`, `Streamlit 1.67.0`, `Plotly 6.9.0`, `Kaleido 0.2.1`), median wall-clock times across three repeats on the bundled 88 MB DepMap Public 26Q3 matrix are **1.01 s** for matrix loading (`0.97–1.33 s`), **0.41 s** for ranking all 18,434 genes (`0.38–0.73 s`), and **0.82 s** (Pearson *r*, `0.66–0.85 s`) or **6.21 s** (Spearman *ρ*, `6.12–6.77 s`) for a pan-cancer genome-wide co-dependency scan across 18,433 partner genes. Correctness is verified by **90 automated tests** covering CSV parsing, missing/infinite handling, exact pairwise correlation parity, UI state isolation, citation integrity, and vector/DPI figure exports via continuous integration on GitHub Actions.

---

## 3. Benchmark Use Cases

### 3.1 Separating Pan-Cancer Essentiality from Lineage-Selective Vulnerability

Upon initial launch on DepMap Public 26Q3 (`n = 1,361` cell lines; `1,712` genes with pan-cancer mean `< −0.5`, range `−5.250` to `+0.398`), CRISPR Score Analyzer renders a default five-gene panel (`KRAS`, `SOX10`, `BRAF`, `ESR1`, `TP53`) anchored by `MYC` (mean `−1.971`, rank `187`, `1.01%`) and `PTEN` (mean `+0.377`, rank `18,428`, `99.97%`) (**Figure 1A**). Across the pan-cancer cohort, `KRAS` sits below the `−0.5` screening reference (mean `−0.729`, rank `1,263`, `6.85%`), whereas tumor suppressor `TP53` lies at the positive extreme (mean `+0.381`, rank `18,430`, `99.98%`), consistent with prior DepMap reports of knockout-induced proliferation advantage in *TP53*-wild-type backgrounds (Tsherniak *et al.*, 2017; Dempster *et al.*, 2021), although somatic mutation calls are not bundled in the default matrix. By contrast, two established lineage-restricted dependencies appear modest or near-neutral in pan-cancer averages: `SOX10` (mean `−0.347`, rank `2,225`, `12.07%`) and `ESR1` (mean `−0.065`, rank `7,347`, `39.86%`).

Inspecting the **Lineage preference** tab separates mean shifts from sample-size and distributional skew (**Figure 1C, D**, noting gene-specific x-axis scales). Relative to its all-lineage average (`−0.347`), `SOX10` exhibits a strong, consistent skin-lineage shift (`mean = −1.487`, `Δ = −1.140`, `median Δ = −1.368`, `78.7% < −0.5`, `n = 75`). In ocular models (`Eye`, `n = 15`), the mean shift (`Δ = −0.535`, `mean = −0.882`) is substantially larger than the median shift (`median Δ = −0.174`), reflecting strong dependency concentrated in a uveal melanoma subset (e.g., `Omm2.5` `−2.752`, `UPMM3` `−2.113`) alongside near-neutral retinoblastoma and retinal pigment epithelial lines. For `ESR1`, `Breast` (`n = 67`) shows a modest lineage-mean shift (`mean = −0.285`, `Δ = −0.220`, `median Δ = −0.061`, `17.9% < −0.5`), followed closely in Figure 1D by a small-sample `Ampulla of Vater` group (`n = 9`, `Δ = −0.174`, `median Δ = −0.059`). Because the bundled matrix stores only `OncotreeLineage` (omitting breast receptor/PAM50 subtype columns), lineage-level means compress ER-positive luminal dependency into a modest average—directly illustrating why paired inspection of **Top cell lines** is essential: that tab immediately surfaces canonical ER-positive models (`HCC1428` [`ACH-000352`, `−1.912`], `T-47D` [`ACH-000147`, `−1.845`], `MDA-MB-415` [`ACH-000876`, `−1.636`]) alongside resistant negative controls, while custom CSV uploads with `OncotreePrimaryDisease` enable subtype-scoped filtering. Complementarily, **Gene-set annotation** overlays mitotic kinases `PLK1` (mean `−2.691`, rank `38`, `0.21%`) and `AURKA` (mean `−1.236`, rank `626`, `3.40%`) onto a seven-gene cell-cycle CDK/cyclin background (`CDK1` `−2.588` through `CCNE1` `−0.194`; **Figure 1B**).

![Figure 1](figures/Figure1.png)

**Figure 1. Multi-gene Chronos gene effect positioning and lineage preference across DepMap Public 26Q3 (`n = 1,361` cell lines, `18,434` genes).** **(A)** Pan-cancer mean Chronos gene effect ranking of five benchmark cancer genes (`KRAS`, `SOX10`, `BRAF`, `ESR1`, `TP53`) relative to `MYC` and `PTEN` reference anchors and the `−0.5` screening reference line. **(B)** Two-layer gene-set annotation highlighting `PLK1` and `AURKA` over seven cell-cycle CDK/cyclin background genes (`CDK1`, `CDK2`, `CDK4`, `CDK6`, `CCNB1`, `CCND1`, `CCNE1`; listed in panel legend). **(C, D)** Diverging lineage preference lollipops (top and bottom 4 lineages with `n ≥ 5`, sample sizes shown on y-axis; note gene-specific x-axis scales) showing mean Chronos difference relative to each gene's all-lineage mean for `SOX10` (**C**, Skin and Eye shifts) and `ESR1` (**D**, Breast alongside small-sample Ampulla of Vater `n = 9`).

### 3.2 Benchmarking Pan-Cancer and Single-Lineage Co-Dependency Exploration

To benchmark exploratory co-dependency workflows in the **Co-dependency** module, we first performed a pan-cancer genome-wide Pearson scan of `MDM2` across all `18,433` other ranked genes (`0.82 s`; **Figure 2A**). The diverging lollipop recapitulates the established `TP53`–`MDM2` regulatory axis: top positive co-dependencies include `PPM1D` (`r = +0.664`), `MDM4` (`r = +0.614`), `USP7` (`r = +0.522`), `DDX31` (`r = +0.466`), and `UBE2D3` (`r = +0.398`, all `n_shared = 1,361`), whereas top anti-correlated partners include `TP53` (`r = −0.702`), `TP53BP1` (`r = −0.589`), `USP28` (`r = −0.531`), `CDKN1A` (`r = −0.530`), `CHEK2` (`r = −0.491`), and `ATM` (`r = −0.461`, all `n_shared = 1,361`). Surfacing `n_shared` directly on the scan output flags `TERF1` (`r = +0.395`, rank 6 positive) as having partial coverage (`n_shared = 228` vs `1,361`), preventing sparse-coverage correlations from being treated equivalently to full-panel hits. Transferring the core six-gene regulatory panel (`TP53`, `CHEK2`, `MDM2`, `MDM4`, `PPM1D`, `USP7`) to multi-gene matrix mode resolves the two-block antagonistic structure (**Figure 2B**), while the integrated cell-line scatter plot of `MDM2` vs `TP53` (`r = −0.702`, `ρ = −0.704`, `n = 1,361`; **Figure 2C**) confirms that the inverse relationship spans CNS/Brain, Lung, Lymphoid, Skin, Bowel, and Breast lineages rather than reflecting a single-tissue artifact.

Conversely, pan-cancer correlations of tissue-enriched oncogenes combine pathway co-essentiality with inter-lineage variance. In a pan-cancer scan of `BRAF`, core MAPK cascade components exhibit moderate-to-strong correlations (`MAPK1` `r = +0.644`, `MAP2K1` `r = +0.391`, `MAP2K2` `r = +0.271`) alongside melanoma lineage markers (`SOX10` `r = +0.544`, `DUSP4` `r = +0.474`; `SOX5` `r = +0.208`). Filtering the sidebar scope exclusively to **Skin** (`n = 75` cell lines; `17,785` partner genes meeting `n_shared ≥ 10`, with `648` sparser genes excluded) and re-running the genome-wide scan yields top seven positive correlations for `MAPK1` (`r = +0.798`), `SOX10` (`+0.633`), `SOX5` (`+0.562`), `MAP2K1` (`+0.552`), `MAP2K2` (`+0.550`), `UBQLN1` (`+0.526`), and `AKIRIN1` (`+0.496`), alongside negative tails led by `STYXL1` (`−0.534`) (**Figure 2D**). Thus, restricting to Skin strengthens canonical RAF–MEK–ERK co-dependencies substantially—and an offline partial-correlation check regressing out `SOX10` within Skin still retains `MAPK1` (`r = +0.665`), `MAP2K1` (`+0.467`), and `MAP2K2` (`+0.398`)—yet the top-ranked list is not a pure MAPK module. The persistence of `SOX10` and `SOX5` reflects Oncotree `Skin` heterogeneity across *SOX10*-high cutaneous melanoma and *SOX10*-low non-melanoma or dedifferentiated skin models (e.g., `HSC-5`, `MCC26`, `RVH421_DAB_R`), while at `n = 75` across `17,785` genes, correlations near `|r| ≈ 0.5` (such as `UBQLN1`, `AKIRIN1`, or `STYXL1`) sit close to genome-wide multiple-testing bounds and should be treated strictly as exploratory candidates rather than confirmed biological interactions.

![Figure 2](figures/Figure2.png)

**Figure 2. Pan-cancer and single-lineage genome-wide co-dependency benchmarks.** **(A)** Pan-cancer genome-wide Pearson *r* scan of `MDM2` across `18,433` genes (`n = 1,361` cell lines), showing the top 7 positive (`PPM1D`, `MDM4`, `USP7`, `DDX31`, `UBE2D3`, `TERF1 [n = 228]`, `PPM1G`) and top 7 negative (`TP53`, `TP53BP1`, `USP28`, `CDKN1A`, `CHEK2`, `ATM`, `UBE2K`) partners. **(B)** Pairwise Pearson correlation matrix for the six-gene `TP53`–`MDM2` regulatory axis. **(C)** Lineage-colored cell-line scatter plot of `MDM2` vs `TP53` Chronos gene effect (`r = −0.702`, Spearman `ρ = −0.704`, `n = 1,361`) with linear trendline. **(D)** Cohort-scoped genome-wide `BRAF` scan restricted to `Skin` cell lines (`n = 75`, `17,785` genes with `n_shared ≥ 10`), showing strengthened MAPK cascade correlations (`MAPK1` `+0.80`, `MAP2K1` `+0.55`, `MAP2K2` `+0.55`) co-ranking with lineage factors (`SOX10` `+0.63`, `SOX5` `+0.56`) and exploratory hits near `|r| ≈ 0.5`.

---

## 4. Scope, Limitations, and Scientific Interpretation

CRISPR Score Analyzer is intentionally designed for transparent exploratory visualization rather than formal statistical hypothesis testing. Accordingly, all summary outputs carry explicit interpretive boundaries:

- **Descriptive metrics only:** Cohort gene ranks, lineage preference mean/median differences, and co-dependency correlation coefficients are descriptive summaries. They do not represent formal statistical tests of selectivity, differential dependency, multiple-testing-corrected significance, or synthetic lethality, and should be followed by appropriate statistical modeling or experimental validation.
- **Reference threshold:** The `−0.5` line serves solely as a visual screening reference on Chronos gene effect scores, not a formal classifier of common-essential genes.
- **Confounding and multiplicity in co-dependency scans:** High correlation between two genes can reflect physical complex membership, shared signaling pathways, residual subtype or cell-state heterogeneity, common-essentiality variance, or small-sample multiplicity when scanning ~18,000 genes in a single lineage. Inspecting shared sample counts (`n_shared`), Spearman *ρ* alongside Pearson *r*, and lineage-colored cell-line scatter plots is essential before inferring direct genetic interaction.
- **Bundled metadata and custom upload assumptions:** The bundled DepMap Public 26Q3 matrix includes `OncotreeLineage` but omits mutation, expression, and PAM50/receptor subtype annotations; users requiring subtype-level stratification should supply `OncotreePrimaryDisease` or custom metadata columns in uploaded CSVs. Uploaded score matrices are assumed to contain one row per unique cell line; duplicate rows are not silently deduplicated.

---

## Acknowledgements and Disclosures

**Data Source:** Reference dependency data are derived from DepMap Public 26Q3 Chronos Gene Effect scores provided by the Broad Institute DepMap Portal (https://depmap.org/portal/).

**AI Assistance Disclosure:** Large language model assistants (Claude, Codex, and related coding agents) were used to assist with software engineering, automated testing, and manuscript editing. All scientific specifications, biological case studies, verification, and final manuscript decisions were reviewed and approved by the author.

**Funding:** This work received no specific external funding.

**Competing Interests:** The author declares no competing financial or non-financial interests.

**Data and Code Availability:** The MIT-licensed source code (v1.2.0), automated test suite, figure reproduction script (`preprint/generate_figures.py`), and matrix build scripts are hosted at https://github.com/ChanghaoKan/crispr-score-analyzer (release tag `v1.2.0`: https://github.com/ChanghaoKan/crispr-score-analyzer/releases/tag/v1.2.0; Zenodo concept DOI: `10.5281/zenodo.19607602`). The hosted web tool is accessible at https://crispr-score-analyzer-szbl-denglab.streamlit.app/. Because release source archives exclude `data/depmap_chronos.npz` under DepMap portal terms, exact reproduction requires cloning tag `v1.2.0` and verifying the bundled matrix against `data/depmap_chronos.provenance.json` (`matrix_sha256`: `ddf274da2aa5f89cf630bc8150dcd3a2c722466bd18ff6fb6d18fe161f8a0c67`).

---

## References

- Dempster, J. M., Boyle, I., Vazquez, F., Root, D. E., Boehm, J. S., Hahn, W. C., Tsherniak, A., & McFarland, J. M. (2021). Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. *Genome Biology*, 22(1), 343. https://doi.org/10.1186/s13059-021-02540-7
- DepMap, Broad. (2026). DepMap Public 26Q3 Chronos Gene Effect dataset. Broad Institute. https://depmap.org/portal/
- Kan, C. (2026). CRISPR Score Analyzer (v1.2.0) [Software]. GitHub / Zenodo. https://doi.org/10.5281/zenodo.19607602 (Release tag `v1.2.0`: https://github.com/ChanghaoKan/crispr-score-analyzer/releases/tag/v1.2.0)
- Meyers, R. M., Bryan, J. G., McFarland, J. M., Weir, B. A., Sizemore, A. E., Xu, H., *et al.* (2017). Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. *Nature Genetics*, 49(12), 1779–1784. https://doi.org/10.1038/ng.3984
- Seal, R. L., Braschi, B., Gray, K., Jones, T. E. M., Tweedie, S., Haim-Vilmovsky, L., & Bruford, E. A. (2023). Genenames.org: the HGNC resources in 2023. *Nucleic Acids Research*, 51(D1), D1003–D1009. https://doi.org/10.1093/nar/gkac888
- Shimada, K., Bachman, J. A., Muhlich, J. L., & Mitchison, T. J. (2021). shinyDepMap, a tool to identify targetable cancer genes and their functional connections from Cancer Dependency Map data. *eLife*, 10, e57116. https://doi.org/10.7554/eLife.57116
- Tsherniak, A., Vazquez, F., Montgomery, P. G., Weir, B. A., Kryukov, G., Cowley, G. S., *et al.* (2017). Defining a Cancer Dependency Map. *Cell*, 170(3), 564–576.e16. https://doi.org/10.1016/j.cell.2017.06.010
- Wainberg, M., Kamber, R. A., Balsubramani, A., Meyers, R. M., Sinnott-Armstrong, N., Hornburg, D., *et al.* (2021). A genome-wide atlas of co-essential modules assigns function to uncharacterized genes. *Nature Genetics*, 53(5), 638–649. https://doi.org/10.1038/s41588-021-00840-z
