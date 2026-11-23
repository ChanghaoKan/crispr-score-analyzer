# CRISPR Score Analyzer: an interactive, cohort-scoped browser platform for DepMap gene essentiality and co-dependency analysis

**Changhao Kan¹,\***

¹ Shenzhen Bay Laboratory, Shenzhen, China  
\* Correspondence: `kanchanghao@szbl.ac.cn`  
ORCID: [0009-0000-4014-5923](https://orcid.org/0009-0000-4014-5923)

**Category:** Bioinformatics (Application Note)

---

## Abstract

**Summary:** Genome-scale CRISPR-Cas9 loss-of-function screens in the Cancer Dependency Map (DepMap) are widely used to prioritize candidate cancer targets, identify lineage-restricted vulnerabilities, and map functional gene networks. In day-to-day research, however, experimental and computational biologists routinely encounter workflow friction when moving beyond single-gene portal lookups—such as locating multi-gene candidate lists on a cohort-wide essentiality curve, separating pan-cancer essentiality from tissue-specific dependency, selecting sensitive and resistant cell-line models for wet-lab validation, performing single-lineage genome-wide co-dependency scans, and exporting publication-ready vector figures alongside reproducible analysis provenance. We present **CRISPR Score Analyzer (v1.2.0)**, an open-source, bilingual (English/Chinese) browser platform built on Python and Streamlit. The application ships directly with the DepMap Public 26Q3 Chronos gene effect matrix (1,361 cancer cell lines × 18,434 genes across 30 Oncotree lineages) together with 40,402 unambiguous HGNC alias mappings, and also supports validated custom CRISPR score CSV uploads. Across four integrated workflows—*Gene ranking*, *Cancer-type comparison*, *Gene-set annotation*, and *Co-dependency*—users can interrogate pan-cancer or lineage-restricted cohorts, inspect lineage preference and top sensitive/resistant cell lines, run chunked vectorized genome-wide Pearson or Spearman co-dependency scans in seconds without external statistical dependencies, and export full-vector PDF/SVG figures, 300-DPI PNGs, and self-contained ZIP provenance bundles. We demonstrate the platform's utility through two biological use cases: distinguishing lineage-selective dependencies (`SOX10` in skin and `ESR1` in breast cancer) from pan-cancer averages, and uncovering both pan-cancer `MDM2`–`TP53` antagonism (`r = -0.702`) and sharpened intra-lineage `BRAF`–MAPK co-dependencies within skin cancer cell lines (`MAPK1` `r = +0.798`, `MAP2K1` `r = +0.552`, `MAP2K2` `r = +0.550`).

**Availability and Implementation:** Source code for v1.2.0 is freely available under the MIT license at https://github.com/ChanghaoKan/crispr-score-analyzer (release tag `v1.2.0`: https://github.com/ChanghaoKan/crispr-score-analyzer/releases/tag/v1.2.0; Zenodo concept DOI: `10.5281/zenodo.19607602`). A live web instance is hosted at https://crispr-score-analyzer-szbl-denglab.streamlit.app/. The bundled DepMap Public 26Q3 matrix is governed by DepMap portal terms.

**Keywords:** DepMap; CRISPR-Cas9 screen; gene essentiality; co-dependency; cancer cell lines; interactive visualization; reproducibility

---

## 1. Introduction

Pooled genome-scale CRISPR-Cas9 loss-of-function screens across hundreds of human cancer cell lines—spearheaded by the Broad Institute Cancer Dependency Map (DepMap) (Tsherniak *et al.*, 2017; Meyers *et al.*, 2017; Dempster *et al.*, 2021)—have become a foundational resource for oncology target prioritization, pathway annotation, and experimental model selection. Using copy-number- and kinetics-corrected Chronos gene effect scores (Dempster *et al.*, 2021), researchers routinely evaluate whether candidate genes from transcriptomic, proteomic, or genetic screens represent broad pan-cancer dependencies or context-specific vulnerabilities.

Although the DepMap Portal (https://depmap.org/portal/) and companion web tools such as shinyDepMap (Shimada *et al.*, 2021) and co-essentiality atlases (Wainberg *et al.*, 2021) provide rich access to dependency datasets, everyday research workflows still expose recurring practical gaps:

1. **Multi-gene batch context:** Positioning a candidate panel of 5–50 genes simultaneously along a cohort-wide essentiality curve—or overlaying a prioritized subset onto a broader pathway background—typically requires repetitive single-gene queries or ad-hoc R/Python scripting.
2. **Bridging pan-cancer summaries to wet-lab model selection:** Pan-cancer average scores frequently obscure lineage-restricted vulnerabilities, while selecting concrete positive (sensitive) and negative (resistant) cell lines within a chosen lineage scope requires cross-referencing multiple views.
3. **Cohort-scoped genome-wide co-dependency discovery:** Standard portal co-dependency lookups compute correlations across the full pan-cancer panel, where shared lineage expression can confound pathway co-essentiality; interactively scanning a query gene against all other genes *within a single cancer lineage or disease subtype* is rarely available in one browser workflow.
4. **Publication export and reproducibility:** Generating clean, vector-editable figures (PDF/SVG) paired with exact numeric tables, data checksums, input-matching records, and software environment metadata remains fragmented across tools.

To address these workflow bottlenecks without requiring programming expertise, we developed **CRISPR Score Analyzer**, an open-source interactive web application designed with explicit draft-vs-committed state isolation, transparent descriptive statistics, bilingual English/Chinese support, and end-to-end provenance export.

---

## 2. Implementation and Core Capabilities

### 2.1 Bundled Reference Data, Custom Uploads, and Symbol Resolution

CRISPR Score Analyzer is implemented in Python (Streamlit, pandas, NumPy, Plotly, Kaleido, and Pillow). To eliminate network latency and external download fragility on startup, the repository bundles a compact `float32` binary representation (`data/depmap_chronos.npz`, 88 MB; maximum absolute `float32` rounding deviation `1.19 × 10⁻⁷` across the first 100 genes) of the **DepMap Public 26Q3 Chronos Gene Effect** matrix comprising **1,361 cancer cell lines** and **18,434 genes** across **30 Oncotree lineages**, paired with SHA-256 and MD5 provenance metadata (`data/depmap_chronos.provenance.json`).

Researchers can alternatively upload custom UTF-8 CRISPR score CSV matrices. To guard against silent misinterpretation of numeric metadata columns, automatic gene detection strictly requires DepMap-style `SYMBOL (Entrez ID)` headers, whereas plain-symbol tables require explicit manual score-column mapping in the sidebar. Duplicate ambiguous gene symbols and repeated metadata columns are automatically excluded and documented in a downloadable `diagnostics.csv` validation report.

Gene-list inputs accept arbitrary delimiters (newlines, spaces, tabs, English/Chinese commas and semicolons), `.txt/.csv/.tsv` file uploads, and pasted DepMap headers (e.g., `MYC (4609)`). An offline HGNC symbol map (`data/gene_aliases.json`, **40,402** unambiguous previous-symbol and alias mappings) (Seal *et al.*, 2023) resolves legacy names (such as `HER2 → ERBB2`) automatically. Matching follows an exception-only interface principle: clean matches remain uncluttered, whereas resolved aliases, removed duplicates, and unmatched entries are surfaced immediately before and after submission.

### 2.2 Modular Analytical Workflows

The interface organizes analyses into four top-level modules governed by a shared sidebar scope (**Cancer types** and optional **Disease / subtype** when present in custom datasets; **Table 1**):

1. **Gene ranking:** Ranks all valid genes within the active cohort by ascending mean finite Chronos score (`rank`, `percentile = 100 × rank / N`, `n_valid`, and missingness counts). Interactive scatter curves render background genes via WebGL (`Scattergl`) for fluid responsiveness while converting all traces to crisp vector paths during PDF/SVG export. Two configurable reference anchors (default `MYC` and `PTEN`, supporting alias resolution) and a dashed `−0.5` screening reference line provide immediate visual orientation, accompanied by collision-avoided transparent gene labels with leader lines.
2. **Cancer-type comparison:** Generates horizontal per-lineage boxplots displaying valid sample counts (`n`) per group, sorted by median score or alphabetically, with optional 2-column grids, single-axis multi-gene comparison, and individual cell-line point overlays annotated with display name and `ModelID`. Submissions matching more than eight genes prompt explicit subset selection while preserving the complete submitted input in metadata.
3. **Gene-set annotation:** Plots a background gene set beneath a prioritized highlight layer on the cohort ranking curve with configurable colors, enabling intuitive pathway-level positioning without conflating visualization with formal enrichment testing.
4. **Co-dependency:** Supports two complementary modes with instant switching between **Pearson *r*** and **Spearman *ρ***:
   - *Genome-wide scan (1 gene):* Correlates a single query gene (selected via a unified searchable dropdown accepting approved symbols and HGNC aliases) against every other ranked gene in the active sidebar scope using exact pairwise-complete finite observations. Chunked vectorized NumPy/pandas routines compute both Pearson and Spearman correlations without requiring SciPy, rendering a diverging top-partner lollipop chart, top-50 positive and negative tables, one-click transfer of top partners to a correlation matrix, and full-genome CSV export.
   - *Multi-gene matrix (2–25 genes):* Computes an interactive pairwise correlation heatmap alongside an integrated **lineage-colored cell-line scatter plot** with ordinary least-squares trendline and both Pearson *r* and Spearman *ρ* statistics for any chosen gene pair.

**Table 1. Comparison of workflow capabilities across CRISPR Score Analyzer and existing DepMap web interfaces.**

| Capability | DepMap Portal | shinyDepMap | CRISPR Score Analyzer |
|---|:---:|:---:|:---:|
| Multi-gene batch placement on cohort ranking curve | Limited (single-gene default) | No (efficacy/selectivity scatter) | **Yes** (ranking + two-layer annotation) |
| Cohort-scoped / single-lineage genome-wide co-dependency scan | Pan-cancer precomputed | Pan-cancer clusters | **Yes** (live Pearson *r* & Spearman *ρ*) |
| Pairwise multi-gene correlation matrix + lineage scatter | Custom explorer steps | Cluster view | **Yes** (2–25 genes, integrated scatter) |
| Lineage preference lollipop/heatmap vs all-lineage mean | Partial (per-gene violin) | No | **Yes** (Top/Bottom 6 lollipop & heatmap) |
| Wet-lab top sensitive / resistant cell-line lookup | Table filtering | Top sensitive list | **Yes** (paired sensitive/resistant + CSV) |
| Custom CRISPR score CSV upload with column diagnostics | No | No | **Yes** (strict schema + `diagnostics.csv`) |
| Offline HGNC alias resolution (40,402 symbols) | Built-in search | No | **Yes** (batch + single-gene + references) |
| Vector PDF/SVG + 300-DPI PNG + SHA-256 ZIP bundle | PNG/CSV separately | PNG/CSV | **Yes** (unified white-paper export + ZIP) |
| Bilingual interface (English / Chinese) & Light/Dark themes | English only | English only | **Yes** (instant switch, shareable URLs) |

### 2.3 Deep-Dive Views, State Integrity, Performance, and Verification

Beneath ranking, comparison, and annotation results, two contextual deep-dive tabs provide immediate follow-up without switching modules: (i) **Lineage preference**, which compares every cancer lineage (with `n_valid ≥ 5`) against each gene's own full-dataset mean (independent of the sidebar filter, outlining active sidebar lineages) in either Top/Bottom-6 lollipop or diverging heatmap layouts alongside median difference and the fraction of lines scoring `< -0.5`; and (ii) **Top cell lines**, listing the most dependent and least dependent models in the active scope for experimental selection.

To prevent accidental mismatch between figures and inputs, gene-list edits remain isolated drafts until explicitly committed: pending edits switch the status dot to amber (`Unsaved edit`) and morph the primary button to **Update chart**, while changing the dataset, column mapping, or sidebar cohort invalidates prior results until re-run. All figure exports enforce white paper backgrounds, crisp marker halos, and transparent label boxes regardless of interface theme, and cache exported bytes by exact figure fingerprint and dimensions. An on-demand **Complete analysis ZIP** packages `selected_results.csv`, `full_rankings.csv`, `plot_data.csv`, `lineage_preference.csv`, `diagnostics.csv`, and `metadata.json` (recording data SHA-256, release version, Zenodo DOI, cohort parameters, input matching, and Python package versions).

On a standard Linux CPU container (Python 3.12), loading the bundled 88 MB DepMap Public 26Q3 matrix takes **2.12 s**, ranking all 18,434 genes takes **0.92 s**, and a pan-cancer genome-wide co-dependency scan across 18,433 partner genes completes in **1.18 s** (Pearson *r*) or **6.29 s** (Spearman *ρ*). Correctness is verified by **90 automated tests** covering CSV parsing, missing/infinite handling, exact pairwise correlation parity, Streamlit UI state transitions, citation integrity, and vector/DPI figure exports via continuous integration on GitHub Actions.

---

## 3. Demonstration and Biological Use Cases

### 3.1 Separating Pan-Cancer Essentiality from Lineage-Selective Vulnerability

Upon initial launch on DepMap Public 26Q3 (`n = 1,361` cell lines; `1,712` genes with pan-cancer mean `< -0.5`, range `−5.250` to `+0.398`), CRISPR Score Analyzer renders a curated five-gene benchmark panel (`KRAS`, `SOX10`, `BRAF`, `ESR1`, `TP53`) anchored by `MYC` (mean `−1.971`, rank `187`, `1.01%`) and `PTEN` (mean `+0.377`, rank `18,428`, `99.97%`) (**Figure 1A**). Across the pan-cancer cohort, `KRAS` sits below the `−0.5` screening reference (mean `−0.729`, rank `1,263`, `6.85%`), whereas tumor suppressor `TP53` lies at the positive extreme (mean `+0.381`, rank `18,430`, `99.98%`), consistent with knockout-induced proliferation advantage in *TP53*-wild-type lines. Critically, two classic lineage drivers appear modest or near-neutral in the pan-cancer mean curve: `SOX10` (mean `−0.347`, rank `2,225`, `12.07%`) and `ESR1` (mean `−0.065`, rank `7,347`, `39.86%`).

Inspecting the **Lineage preference** tab immediately resolves this masking effect (**Figure 1C, D**). Relative to its all-lineage average (`−0.347`), `SOX10` exhibits a pronounced skin-lineage shift (`mean = −1.487`, `Δ = −1.140`, median `Δ = −1.368`, with `78.7%` of `n = 75` skin lines `< −0.5`), followed by ocular models (`mean = −0.882`, `Δ = −0.535`, `n = 15`), matching its neural-crest master-regulator role in cutaneous and uveal melanoma. Similarly, `ESR1` shifts selectively toward breast cancer (`mean = −0.285`, `Δ = −0.220`, `n = 67`, `17.9% < −0.5`, reflecting the ER-positive luminal fraction), while the **Top cell lines** tab surfaces canonical ER-positive experimental models—`HCC1428` (`ACH-000352`, `−1.912`), `T-47D` (`ACH-000147`, `−1.845`), and `MDA-MB-415` (`ACH-000876`, `−1.636`)—paired with resistant negative controls for assay design. Complementarily, **Gene-set annotation** overlays mitotic kinases `PLK1` (mean `−2.691`, rank `38`, `0.21%`) and `AURKA` (mean `−1.236`, rank `626`, `3.40%`) onto a seven-gene CDK/cyclin background (`CDK1` `−2.588` through `CCNE1` `−0.194`) in a single vector figure (**Figure 1B**).

![Figure 1](figures/Figure1.png)

**Figure 1. Multi-gene essentiality positioning and lineage preference across DepMap Public 26Q3 (`n = 1,361` cell lines, `18,434` genes).** **(A)** Pan-cancer mean Chronos score ranking locating five canonical cancer genes (`KRAS`, `SOX10`, `BRAF`, `ESR1`, `TP53`) relative to `MYC` and `PTEN` reference anchors and the `−0.5` screening line. **(B)** Two-layer gene-set annotation highlighting mitotic kinases (`PLK1`, `AURKA`) over a cell-cycle CDK/cyclin background (`CDK1`, `CDK2`, `CDK4`, `CDK6`, `CCNB1`, `CCND1`, `CCNE1`). **(C, D)** Diverging lineage preference lollipops (top and bottom 4 lineages with `n ≥ 5`) showing mean Chronos difference relative to each gene's all-lineage average for `SOX10` (**C**, strong Skin/Eye selectivity) and `ESR1` (**D**, Breast enrichment).

### 3.2 Pan-Cancer and Single-Lineage Co-Dependency Discovery

To illustrate functional network discovery in the **Co-dependency** module, we first performed a pan-cancer genome-wide Pearson scan of `MDM2` across all `18,433` other ranked genes (`1.18 s`; **Figure 2A**). The diverging lollipop cleanly recovers the symmetric `p53–MDM2` regulatory axis: top positive co-dependencies include `PPM1D` (`r = +0.664`), `MDM4` (`r = +0.614`), `USP7` (`r = +0.522`), `DDX31` (`r = +0.466`), and `UBE2D3` (`r = +0.398`), whereas top anti-correlated partners include `TP53` (`r = −0.702`), `TP53BP1` (`r = −0.589`), `USP28` (`r = −0.531`), `CDKN1A` (`r = −0.530`), `CHEK2` (`r = −0.491`), and `ATM` (`r = −0.461`). Transferring the core six-gene regulatory panel (`TP53`, `CHEK2`, `MDM2`, `MDM4`, `PPM1D`, `USP7`) to multi-gene matrix mode displays the crisp two-block antagonistic structure (**Figure 2B**), while the integrated cell-line scatter plot of `MDM2` vs `TP53` (`r = −0.702`, `ρ = −0.704`, `n = 1,361`; **Figure 2C**) verifies that the inverse relationship holds broadly across CNS/Brain, Lung, Lymphoid, Skin, Bowel, and Breast lineages rather than reflecting a single-tissue artifact.

Conversely, pan-cancer correlations of tissue-enriched oncogenes frequently mix pathway interactions with inter-lineage variance. In a pan-cancer scan of `BRAF`, core MAPK cascade components exhibit moderate-to-strong correlations (`MAPK1` `r = +0.644`, `MAP2K1` `r = +0.391`, `MAP2K2` `r = +0.271`) alongside melanoma lineage markers (`SOX10` `r = +0.544`, `DUSP4` `r = +0.474`). Filtering the sidebar scope exclusively to **Skin** (`n = 75` cell lines; `17,785` genes passing finite-observation checks) and re-running the genome-wide scan sharpens canonical RAF–MEK–ERK pathway co-dependency dramatically (**Figure 2D**): `MAPK1` rises to `r = +0.798`, `MAP2K1` increases from `+0.391` to `+0.552`, and `MAP2K2` more than doubles from `+0.271` to `+0.550`—directly demonstrating the utility of cohort-scoped co-dependency scanning for dissecting intra-lineage signaling dependencies.

![Figure 2](figures/Figure2.png)

**Figure 2. Pan-cancer and single-lineage genome-wide co-dependency discovery.** **(A)** Pan-cancer genome-wide Pearson *r* scan of `MDM2` across `18,433` genes (`n = 1,361` cell lines), showing the top 7 positive (`PPM1D`, `MDM4`, `USP7`, `DDX31`, `UBE2D3`, `TERF1`, `PPM1G`) and top 7 negative (`TP53`, `TP53BP1`, `USP28`, `CDKN1A`, `CHEK2`, `ATM`, `UBE2K`) partners. **(B)** Pairwise Pearson correlation matrix for the six-gene `p53–MDM2` regulatory axis. **(C)** Lineage-colored cell-line scatter plot of `MDM2` vs `TP53` (`r = −0.702`, Spearman `ρ = −0.704`, `n = 1,361`) with linear trendline, confirming pan-lineage antagonism. **(D)** Cohort-scoped genome-wide `BRAF` scan restricted to `Skin` cell lines (`n = 75`, `17,785` genes), demonstrating sharpened intra-lineage MAPK cascade co-dependencies (`MAPK1` `+0.80`, `MAP2K1` `+0.55`, `MAP2K2` `+0.55`).

---

## 4. Scope, Limitations, and Scientific Interpretation

CRISPR Score Analyzer is intentionally designed for transparent exploratory visualization rather than black-box statistical inference. Accordingly, all summary outputs carry explicit interpretive boundaries:

- **Descriptive metrics only:** Cohort gene ranks, lineage preference mean/median differences, and co-dependency correlation coefficients are descriptive summaries. They do not represent formal statistical tests of selectivity, differential dependency, or synthetic lethality, and should be followed by appropriate statistical modeling or experimental validation.
- **Reference threshold:** The `−0.5` line serves solely as a visual screening reference on Chronos gene effect scores, not a formal classifier of common-essential genes.
- **Confounding in co-dependency scans:** High correlation between two genes can reflect physical complex membership, shared signaling pathways, residual subtype composition, or common-essentiality variance. Inspecting shared sample counts (`n_shared`), Spearman *ρ* alongside Pearson *r*, and lineage-colored cell-line scatter plots is essential before inferring direct genetic interaction.
- **Custom upload assumptions:** Uploaded score matrices are assumed to contain one row per unique cell line; duplicate rows are not silently deduplicated.

---

## Acknowledgements and Disclosures

**Data Source:** Reference dependency data are derived from DepMap Public 26Q3 Chronos Gene Effect scores provided by the Broad Institute DepMap Portal (https://depmap.org/portal/).

**AI Assistance Disclosure:** Large language model coding assistants (Claude, Codex, and agentic development tooling) were used to assist with software engineering, interface refinement, automated test expansion, and manuscript editing. All scientific requirements, biological case studies, code verification, and final manuscript decisions were reviewed and approved by the author.

**Funding:** This work received no specific external funding.

**Competing Interests:** The author declares no competing financial or non-financial interests.

**Data and Code Availability:** The MIT-licensed source code (v1.2.0), automated test suite, and build scripts are hosted at https://github.com/ChanghaoKan/crispr-score-analyzer (release tag `v1.2.0`: https://github.com/ChanghaoKan/crispr-score-analyzer/releases/tag/v1.2.0; Zenodo concept DOI: `10.5281/zenodo.19607602`). The hosted web tool is accessible at https://crispr-score-analyzer-szbl-denglab.streamlit.app/.

---

## References

1. Dempster, J. M., Boyle, I., Vazquez, F., Root, D. E., Boehm, J. S., Hahn, W. C., Tsherniak, A., & McFarland, J. M. (2021). Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. *Genome Biology*, 22(1), 343. https://doi.org/10.1186/s13059-021-02540-7
2. DepMap, Broad. (2026). DepMap Public 26Q3 Chronos Gene Effect dataset. Broad Institute. https://depmap.org/portal/
3. Kan, C. (2026). CRISPR Score Analyzer (v1.2.0) [Software]. GitHub / Zenodo. https://doi.org/10.5281/zenodo.19607602 (Release tag `v1.2.0`: https://github.com/ChanghaoKan/crispr-score-analyzer/releases/tag/v1.2.0)
4. Meyers, R. M., Bryan, J. G., McFarland, J. M., Weir, B. A., Sizemore, A. E., Xu, H., *et al.* (2017). Computational correction of copy number effect improves specificity of CRISPR-Cas9 essentiality screens in cancer cells. *Nature Genetics*, 49(12), 1779–1784. https://doi.org/10.1038/ng.3984
5. Seal, R. L., Braschi, B., Gray, K., Jones, T. E. M., Tweedie, S., Haim-Vilmovsky, L., & Bruford, E. A. (2023). Genenames.org: the HGNC resources in 2023. *Nucleic Acids Research*, 51(D1), D1003–D1009. https://doi.org/10.1093/nar/gkac888
6. Shimada, K., Bachman, J. A., Muhlich, J. L., & Mitchison, T. J. (2021). shinyDepMap, a tool to identify targetable cancer genes and their functional connections from Cancer Dependency Map data. *eLife*, 10, e57116. https://doi.org/10.7554/eLife.57116
7. Tsherniak, A., Vazquez, F., Montgomery, P. G., Weir, B. A., Kryukov, G., Cowley, G. S., *et al.* (2017). Defining a Cancer Dependency Map. *Cell*, 170(3), 564–576.e16. https://doi.org/10.1016/j.cell.2017.06.010
8. Wainberg, M., Kamber, R. A., Balsubramani, A., Meyers, R. M., Sinnott-Armstrong, N., Hornburg, D., *et al.* (2021). A genome-wide atlas of co-essential modules assigns function to uncharacterized genes. *Nature Genetics*, 53(5), 638–649. https://doi.org/10.1038/s41588-021-00840-z
