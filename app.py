"""
CRISPR Score Analyzer - Public Gene Essentiality Analysis Tool
A multilingual, theme-aware interactive platform for DepMap CRISPR data.

Author: Deng Lab
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import io
import re
import copy
import hashlib
from html import escape
import json
from datetime import datetime, timezone
from importlib.metadata import version
from analysis_core import (analyze_dataset, parse_gene_text, parse_gene_upload,
                           match_genes, find_lineage_column, get_lineage_frame,
                           lineage_summary, make_analysis_bundle, read_score_csv)
from figure_export import export_fingerprint, render_figure_bytes, MAX_PNG_PIXELS

# =============================================================================
# 页面配置
# =============================================================================
st.set_page_config(
    page_title="CRISPR Score Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="auto"
)

# =============================================================================
# 数据源配置 (Hugging Face)
# =============================================================================
HF_REPO_ID = "ChanghaoKan/crispr-depmap"
HF_FILENAME = "CRISPR_(DepMap_Public_25Q3+Score,_Chronos)_subsetted.csv"
HF_REVISION = "8400bf566411ed1df4e0784fdb9e97fdbaa371fd"
USE_HUGGINGFACE = True
DATA_VERSION = "DepMap Public 25Q3"
SCORE_TYPE = "Chronos Gene Effect"
ESSENTIALITY_THRESHOLD = -0.5

# 暂时不在公开网页中展示基因 × 药物分析。底层分析代码仍保留，便于后续恢复。
ENABLE_GENE_DRUG_UI = False

# ---- Gene × Drug correlation module (DepMap 26Q1 + GDSC2) ----
# 同一个 HF 数据集仓库内放置以下三个文件即可（行/列键已对齐到 DepMap ModelID）：
#   GDSC2: index=ModelID(ACH-xxxxxx), columns=CompoundID(DPC-xxxxxx), values=AUC
#   CRISPR(26Q1): index=ModelID, columns="GENE (entrez)", values=Gene Effect(Chronos)
#   Compounds: 列含 CompoundID/CompoundName/GeneSymbolOfTargets/TargetOrMechanism
HF_GDSC_FILENAME = "GDSC2_AUC_Matrix.csv"
HF_CRISPR26Q1_FILENAME = "CRISPRGeneEffect_26Q1.csv"
HF_COMPOUNDS_FILENAME = "PortalCompounds.csv"
HF_MODEL_FILENAME = "Model.csv"  # 提供 ModelID -> OncotreeLineage/PrimaryDisease

# 肿瘤类型分层：粗(lineage) / 细(primary disease)
LINEAGE_COL_COARSE = "OncotreeLineage"
LINEAGE_COL_FINE = "OncotreePrimaryDisease"
MIN_N_PER_GROUP = 10  # 分层相关每组最少细胞系数

# 默认演示用的基因 / 药物关键词（用户可改）
DEFAULT_CORR_GENE = "DUSP6"
DEFAULT_CORR_DRUG_QUERY = "PARP"

# =============================================================================
# Citation / DOI 配置
# =============================================================================
ZENODO_DOI = "10.5281/zenodo.19607603"
TOOL_VERSION = "v1.0"
TOOL_AUTHORS = "Kan, C."
TOOL_YEAR = "2026"
GITHUB_URL = "https://github.com/ChanghaoKan/crispr-score-analyzer"

# =============================================================================
# 国际化
# =============================================================================
TRANSLATIONS = {
    'en': {
        'app_title': 'CRISPR Score Analyzer',
        'app_subtitle': 'Explore gene dependencies across cancer cell lines.',
        'hero_kicker': 'DEPMAP · CRISPR SCREENING',
        'sidebar_settings': 'Settings',
        'language': 'Language',
        'theme': 'Theme',
        'light': 'Light',
        'dark': 'Dark',
        'data_source': 'Data Source',
        'upload_custom': 'Upload custom data (optional)',
        'upload_csv': 'Upload CRISPR Score CSV',
        'reference_genes': 'Reference Genes',
        'essential': 'Reference A',
        'nonessential': 'Reference B',
        'display_settings': 'Display',
        'show_labels': 'Show gene labels',
        'point_size': 'Point size',
        'export_size': 'Export Size',
        'export_height': 'Image height (px)',
        'cell_lines': 'Cell lines',
        'gene_count': 'Genes',
        'essential_genes': 'Mean score < {threshold}',
        'score_range': 'Score range',
        'resources': 'Citation & data source',
        'custom_dataset': 'Custom uploaded dataset',
        'score_guide': 'Lower score = stronger gene dependency',
        'essential_status': 'Mean-score screen',
        'essential_yes': 'Mean below cutoff',
        'essential_no': 'Mean not below cutoff',
        'first_eight_only': 'Showing the first 8 genes to keep the figure readable.',
        'lineage_missing': 'No lineage column was found in this dataset.',
        'tab1': 'Gene ranking',
        'tab2': 'Cancer-type comparison',
        'tab3': 'Gene-set annotation',
        'tab4': '🔗 Gene × Drug',
        'gene_ranking_title': 'Gene dependency ranking',
        'gene_ranking_desc': 'Locate your genes of interest in genome-wide CRISPR screens.',
        'input_target_genes': 'Input target genes',
        'input_method': 'Input method',
        'input_direct': 'Direct input',
        'input_file': 'Upload gene list',
        'gene_list': 'Gene list',
        'gene_list_help': 'One gene per line, or comma-separated',
        'matched': 'Matched',
        'not_found': 'Not found',
        'export_title': 'Export high-quality figure',
        'gene_details': 'Gene details',
        'boxplot_title': 'CRISPR Score Distribution by Cancer Type',
        'multilayer_title': 'Gene-set annotation',
        'bg_gene_set': 'Background gene set',
        'hl_gene_set': 'Highlight genes',
        'bg_color': 'Background color',
        'hl_color': 'Highlight color',
        'download_pdf': '📄 PDF (vector)',
        'download_png': '🖼️ PNG (300 DPI)',
        'download_svg': '✏️ SVG (editable)',
        'download_csv': '📊 Download data (CSV)',
        'download_hint': '💡 Click a button to generate the file. Recommend **PDF** for publications. Use **SVG** for Illustrator/Inkscape editing.',
        'no_data_warn': '⚠️ Please configure a data source or upload data',
        'loading_upload': 'Loading uploaded data...',
        'loading_hf': 'Loading the reference dataset...',
        'loaded': '✅ Loaded',
        'acknowledgements': 'Acknowledgements',
        'data_from': 'Data Source',
        'dev_with': 'Development Assistance',
        'ai_dev': 'AI-assisted development',
        'citation': '📚 How to Cite',
        'cite_this_tool': 'Cite this tool',
        'copy_bibtex': 'Copy BibTeX',
        'doi_pending': 'DOI pending — please cite by URL until release',
        # ---- Tab4: Gene x Drug correlation ----
        'corr_title': 'Gene Dependency × Drug Sensitivity Correlation',
        'corr_desc': 'Test whether cell lines more dependent on a gene (CRISPR) are also more sensitive to a drug (GDSC2 AUC), aligned by DepMap ModelID.',
        'corr_gene_label': 'Gene of interest (CRISPR dependency)',
        'corr_gene_help': 'e.g. DUSP6 — matched against the 26Q1 CRISPR gene-effect matrix.',
        'corr_drug_search': 'Search drug (name or target)',
        'corr_drug_search_help': 'e.g. PARP (matches all PARP inhibitors via target), or comma-separated names: olaparib, talazoparib, niraparib.',
        'corr_drug_select': 'Select compound',
        'corr_gdsc2_only': 'Only compounds with GDSC2 data',
        'corr_gdsc2_only_help': 'PortalCompounds lists many compounds from various sources; only some have GDSC2 AUC values. Keep this on to show only analyzable drugs.',
        'corr_no_drug': 'No compound matches your search.',
        'corr_run': 'Run correlation',
        'corr_loading': 'Loading correlation datasets (GDSC2 + CRISPR 26Q1 + Compounds)...',
        'corr_gene_not_found': 'Gene not found in CRISPR matrix.',
        'corr_drug_not_found': 'Compound not found in GDSC2 matrix.',
        'corr_too_few': 'Too few overlapping cell lines for a reliable estimate (n < 10).',
        'corr_result_dir_pos': 'Positive ρ: cell lines more dependent on this gene (lower Gene Effect) tend to be MORE sensitive to the drug (lower AUC) — consistent with a shared-vulnerability hypothesis.',
        'corr_result_dir_neg': 'Negative ρ: cell lines more dependent on this gene tend to be LESS sensitive to the drug — opposite to the shared-vulnerability hypothesis.',
        'corr_result_dir_ns': 'Correlation is not statistically significant (p ≥ 0.05); no clear association in this pan-cancer set.',
        'corr_axis_x': 'CRISPR Gene Effect (lower = more dependent)',
        'corr_axis_y': 'GDSC2 AUC (lower = more sensitive)',
        'corr_stat_n': 'Overlapping cell lines',
        'corr_stat_rho': "Spearman ρ",
        'corr_stat_p': 'p-value',
        'corr_caveat': '⚠️ Pan-cancer correlation only. Lineage and mutation background (e.g. BRCA/HR status for PARP inhibitors) are NOT controlled here — a significant ρ may reflect confounding. Interpret as association, not mechanism.',
        'corr_lineage_gran': 'Cancer-type granularity',
        'corr_gran_coarse': 'Lineage (broad)',
        'corr_gran_fine': 'Primary disease (fine)',
        'corr_restrict_lineage': 'Restrict to one cancer type',
        'corr_all_lineages': 'All cancer types',
        'corr_no_model': 'Model.csv not found in the dataset repo — lineage stratification is disabled. Upload Model.csv to enable per-cancer-type analysis.',
        'corr_strat_title': 'Per-lineage correlation (confounding check)',
        'corr_strat_desc': 'Spearman ρ computed within each cancer type (n ≥ 10). If the overall correlation holds within lineages, it is less likely to be driven purely by cancer-type composition. Red points: BH-FDR < 0.05.',
        'corr_strat_insufficient': 'Not enough cancer types with ≥10 cell lines for stratified analysis.',
        'corr_forest_sub': 'Point size ∝ n; red = BH-FDR < 0.05',
        'corr_tbl_lineage': 'Cancer type',
    },
    'zh': {
        'app_title': 'CRISPR 基因依赖分析',
        'app_subtitle': '探索癌症细胞系中的基因依赖性。',
        'hero_kicker': 'DEPMAP · CRISPR 筛选',
        'sidebar_settings': '设置',
        'language': '语言',
        'theme': '主题',
        'light': '浅色',
        'dark': '深色',
        'data_source': '数据来源',
        'upload_custom': '上传自定义数据（可选）',
        'upload_csv': '上传 CRISPR Score CSV',
        'reference_genes': '参考基因',
        'essential': '参考基因 A',
        'nonessential': '参考基因 B',
        'display_settings': '显示设置',
        'show_labels': '显示基因名标签',
        'point_size': '点大小',
        'export_size': '导出尺寸',
        'export_height': '图片高度 (px)',
        'cell_lines': '细胞系',
        'gene_count': '基因数',
        'essential_genes': '平均分 < {threshold}',
        'score_range': '分数范围',
        'resources': '引用与数据来源',
        'custom_dataset': '自定义上传数据',
        'score_guide': '分数越低，基因依赖越强',
        'essential_status': '平均分筛选',
        'essential_yes': '平均分低于阈值',
        'essential_no': '平均分未低于阈值',
        'first_eight_only': '为保证图表清晰，仅展示前 8 个基因。',
        'lineage_missing': '该数据中未找到 lineage 列。',
        'tab1': '基因排名',
        'tab2': '癌种比较',
        'tab3': '基因集标注',
        'tab4': '🔗 基因×药物',
        'gene_ranking_title': '基因依赖排名',
        'gene_ranking_desc': '在全基因组 CRISPR 筛选数据中定位您关注的基因。',
        'input_target_genes': '输入目标基因',
        'input_method': '输入方式',
        'input_direct': '直接输入',
        'input_file': '上传基因列表',
        'gene_list': '基因列表',
        'gene_list_help': '每行一个基因名，或用逗号分隔',
        'matched': '匹配成功',
        'not_found': '未找到',
        'export_title': '导出高质量图片',
        'gene_details': '基因详细信息',
        'boxplot_title': '按癌症类型的 CRISPR Score 分布',
        'multilayer_title': '基因集标注',
        'bg_gene_set': '背景基因集',
        'hl_gene_set': '高亮基因',
        'bg_color': '背景颜色',
        'hl_color': '高亮颜色',
        'download_pdf': '📄 PDF (矢量)',
        'download_png': '🖼️ PNG (300 DPI)',
        'download_svg': '✏️ SVG (可编辑)',
        'download_csv': '📊 下载数据表 (CSV)',
        'download_hint': '💡 点击按钮生成文件。论文投稿推荐 **PDF**（矢量图）。需要编辑选 **SVG**（可在 Illustrator/Inkscape 中修改）。',
        'no_data_warn': '⚠️ 请配置数据源或上传数据文件',
        'loading_upload': '正在加载上传的数据...',
        'loading_hf': '正在加载参考数据集...',
        'loaded': '✅ 已加载',
        'acknowledgements': '致谢',
        'data_from': '数据来源',
        'dev_with': '开发协助',
        'ai_dev': 'AI 辅助开发',
        'citation': '📚 引用方式',
        'cite_this_tool': '引用本工具',
        'copy_bibtex': '复制 BibTeX',
        'doi_pending': 'DOI 申请中 — 正式发布前请用网址引用',
        # ---- Tab4: 基因×药物相关 ----
        'corr_title': '基因依赖性 × 药物敏感性相关分析',
        'corr_desc': '检验：对某基因依赖性越强（CRISPR）的细胞系，是否也越敏感于某药物（GDSC2 AUC）。按 DepMap ModelID 对齐。',
        'corr_gene_label': '目标基因（CRISPR 依赖性）',
        'corr_gene_help': '如 DUSP6，在 26Q1 CRISPR 基因效应矩阵中匹配。',
        'corr_drug_search': '搜索药物（名称或靶点）',
        'corr_drug_search_help': '如 PARP（按靶点匹配所有 PARP 抑制剂），或逗号分隔多个药名：olaparib, talazoparib, niraparib。',
        'corr_drug_select': '选择化合物',
        'corr_gdsc2_only': '仅显示 GDSC2 有数据的化合物',
        'corr_gdsc2_only_help': 'PortalCompounds 收录了多来源的化合物，但只有部分在 GDSC2 里有 AUC 数据。保持勾选可只显示能分析的药物。',
        'corr_no_drug': '没有匹配的化合物。',
        'corr_run': '运行相关分析',
        'corr_loading': '加载相关分析数据集（GDSC2 + CRISPR 26Q1 + 化合物表）...',
        'corr_gene_not_found': '在 CRISPR 矩阵中未找到该基因。',
        'corr_drug_not_found': '在 GDSC2 矩阵中未找到该化合物。',
        'corr_too_few': '重叠细胞系太少，结果不可靠（n < 10）。',
        'corr_result_dir_pos': 'ρ 为正：对该基因依赖越强（Gene Effect 越低）的细胞系，往往对该药更敏感（AUC 越低）——与「共享脆弱性」假设方向一致。',
        'corr_result_dir_neg': 'ρ 为负：对该基因依赖越强的细胞系，往往对该药更不敏感——与「共享脆弱性」假设方向相反。',
        'corr_result_dir_ns': '相关不显著（p ≥ 0.05）；在该泛癌集合中无明确关联。',
        'corr_axis_x': 'CRISPR Gene Effect（越低=依赖越强）',
        'corr_axis_y': 'GDSC2 AUC（越低=越敏感）',
        'corr_stat_n': '重叠细胞系数',
        'corr_stat_rho': "Spearman ρ",
        'corr_stat_p': 'p 值',
        'corr_caveat': '⚠️ 仅为泛癌相关。此处未控制谱系与突变背景（如 PARP 抑制剂的 BRCA/HR 状态）——显著的 ρ 可能来自混杂。结论应表述为关联，而非机制。',
        'corr_lineage_gran': '肿瘤类型粒度',
        'corr_gran_coarse': '谱系（粗）',
        'corr_gran_fine': '原发病种（细）',
        'corr_restrict_lineage': '仅分析某一肿瘤类型',
        'corr_all_lineages': '全部肿瘤类型',
        'corr_no_model': '数据集仓库中未找到 Model.csv —— 肿瘤分层功能已禁用。上传 Model.csv 即可启用按癌种分析。',
        'corr_strat_title': '分层相关（混杂检查）',
        'corr_strat_desc': '在每个肿瘤类型内部分别计算 Spearman ρ（n ≥ 10）。若整体相关在各谱系内部依然成立，则更不可能纯粹由癌种构成驱动。红点：BH-FDR < 0.05。',
        'corr_strat_insufficient': '细胞系数 ≥10 的肿瘤类型不足，无法做分层分析。',
        'corr_forest_sub': '点大小 ∝ n；红色 = BH-FDR < 0.05',
        'corr_tbl_lineage': '肿瘤类型',
    }
}


def t(key: str) -> str:
    lang = st.session_state.get('lang', 'en')
    return TRANSLATIONS[lang].get(key, key)


# =============================================================================
# Session State
# =============================================================================
if st.session_state.get('lang') not in TRANSLATIONS:
    st.session_state.lang = 'en'
if st.session_state.get('theme') not in {'light', 'dark'}:
    st.session_state.theme = 'light'


# =============================================================================
# 主题
# =============================================================================
THEMES = {
    'light': {
        'bg': '#f5f7fa',
        'bg_secondary': '#edf2f7',
        'bg_card': '#ffffff',
        'text': '#172231',
        'text_muted': '#607080',
        'border': '#dce4ec',
        'accent': '#245f96',
        'accent_hover': '#1f568d',
        'success': '#1f8a67',
        'danger': '#b54552',
        'shadow': '0 1px 3px rgba(31, 54, 78, 0.04)',
        'shadow_hover': '0 12px 30px rgba(31, 54, 78, 0.10)',
        'plot_bg': '#ffffff',
        'plot_text': '#172231',
        'plot_axis': '#2f3d4b',
        'plot_reference': 'rgba(47,61,75,0.32)',
        'plot_grid': 'rgba(47,61,75,0.06)',
        'plot_scatter_bg': 'rgba(126,139,151,0.23)',
    },
    'dark': {
        'bg': '#0f1720',
        'bg_secondary': '#131e29',
        'bg_card': '#182430',
        'text': '#edf3f8',
        'text_muted': '#9eafbf',
        'border': '#2a3a49',
        'accent': '#76a9da',
        'accent_hover': '#95bee5',
        'success': '#55c39a',
        'danger': '#ee8490',
        'shadow': '0 10px 28px rgba(0,0,0,0.22)',
        'shadow_hover': '0 14px 34px rgba(0,0,0,0.30)',
        'plot_bg': '#182430',
        'plot_text': '#edf3f8',
        'plot_axis': '#b6c3cf',
        'plot_reference': 'rgba(182,195,207,0.32)',
        'plot_grid': 'rgba(255,255,255,0.08)',
        'plot_scatter_bg': 'rgba(170,183,195,0.22)',
    }
}


def get_theme():
    return THEMES.get(st.session_state.get('theme'), THEMES['light'])


def inject_css():
    th = get_theme()
    st.markdown(f"""
    <style>
        html, body, [class*="css"] {{
            font-family: Inter, -apple-system, BlinkMacSystemFont, "Segoe UI",
                         "PingFang SC", "Microsoft YaHei", sans-serif;
        }}
        .stApp {{ background: {th['bg']}; color: {th['text']}; }}
        [data-testid="stMainBlockContainer"] {{
            max-width: 1280px; padding: 1.25rem 2rem 2rem;
        }}
        [data-testid="stVerticalBlock"] {{ gap: 0.85rem; }}
        section[data-testid="stSidebar"] {{
            background: {th['bg_secondary']}; border-right: 1px solid {th['border']};
        }}
        [data-testid="stSidebarUserContent"] {{ padding-top: 1.25rem; }}
        h1, h2, h3, h4, h5, label {{ color: {th['text']} !important; }}
        h2 {{ font-size: 1.15rem !important; }}
        h3 {{ font-size: 1.1rem !important; padding-bottom: 0.25rem !important; }}
        p, .stCaption {{ line-height: 1.5; }}
        [data-testid="stCaptionContainer"] {{ color: {th['text_muted']}; }}
        .hero-shell {{
            display: flex; align-items: center; justify-content: space-between;
            gap: 1rem; padding: 0.35rem 0 0.7rem;
        }}
        .main-header {{
            font-size: clamp(1.35rem, 2.3vw, 1.85rem) !important;
            font-weight: 720; letter-spacing: -0.045em; line-height: 1.2;
            margin: 0 !important; padding: 0 !important;
        }}
        .sub-header {{ margin: 0.35rem 0 0; color: {th['text_muted']}; font-size: 0.9rem; }}
        .hero-version {{
            flex-shrink: 0; border: 1px solid {th['border']}; border-radius: 6px;
            color: {th['accent']}; background: {th['bg_card']};
            font-size: 0.75rem; font-weight: 650; padding: 0.35rem 0.6rem;
        }}
        .context-strip {{
            display: flex; flex-wrap: wrap; align-items: center; gap: 0.65rem 1.25rem;
            padding: 0.8rem 1rem; border: 1px solid {th['border']};
            background: {th['bg_card']}; border-radius: 8px;
            color: {th['text_muted']}; font-size: 0.82rem;
        }}
        .context-strip strong {{ color: {th['text']}; font-size: 1rem; margin-right: 0.25rem; }}
        .context-scope {{ margin-left: auto; overflow-wrap: anywhere; }}
        .status-strip {{
            padding: 0.65rem 0.85rem; border-left: 3px solid {th['success']};
            background: {th['bg_card']}; border-radius: 5px;
            color: {th['text_muted']}; font-size: 0.8rem; overflow-wrap: anywhere;
        }}
        .status-strip strong {{ color: {th['text']}; }}
        .gene-tag {{
            display: inline-block; color: {th['accent']}; background: {th['accent']}0d;
            border: 1px solid {th['accent']}24; border-radius: 5px;
            padding: 0.18rem 0.5rem; margin: 0.12rem 0.3rem 0.12rem 0;
            font-size: 0.82rem; font-weight: 600;
        }}
        [data-testid="stButtonGroup"] {{ margin: 0.1rem 0 0.2rem; }}
        [data-testid="stButtonGroup"] button {{ min-height: 2.45rem; font-weight: 600; }}
        [data-testid="stButtonGroup"] button[aria-checked="true"],
        [data-testid="stButtonGroup"] button[aria-pressed="true"] {{
            background: {th['accent']}14; color: {th['accent']}; border-color: {th['accent']};
        }}
        [data-testid="stVerticalBlockBorderWrapper"] {{
            border-color: {th['border']} !important; border-radius: 10px !important;
        }}
        [data-testid="stVerticalBlockBorderWrapper"] > div {{ background: {th['bg_card']}; }}
        [data-testid="stPlotlyChart"] {{
            background: {th['plot_bg']}; border: 1px solid {th['border']};
            border-radius: 10px; overflow: hidden;
        }}
        button {{ border-radius: 7px !important; }}
        .stButton button, .stDownloadButton button, [data-testid="stPopover"] button {{
            min-height: 2.4rem; font-weight: 600;
        }}
        button[kind="primary"] {{ background: {th['accent']}; border-color: {th['accent']}; color: {th['bg_card']}; }}
        button[kind="primary"]:hover {{ background: {th['accent_hover']}; border-color: {th['accent_hover']}; }}
        .stTextInput input, .stTextArea textarea {{
            background: {th['bg_card']} !important; color: {th['text']} !important;
            border-color: {th['border']} !important; line-height: 1.55;
        }}
        .stSelectbox div[data-baseweb="select"] > div,
        .stMultiSelect div[data-baseweb="select"] > div {{
            background: {th['bg_card']}; color: {th['text']}; border-color: {th['border']};
        }}
        .stSelectbox div[data-baseweb="select"] span {{ color: {th['text']}; }}
        div[data-baseweb="popover"], div[role="listbox"] {{ background: {th['bg_card']}; color: {th['text']}; }}
        [data-testid="stExpander"] summary {{ background: {th['bg_card']}; font-weight: 550; }}
        [data-testid="stDataFrame"] {{ border-radius: 8px; overflow: hidden; border: 1px solid {th['border']}; }}
        .footer-card {{ padding: 1rem; background: {th['bg_card']}; border: 1px solid {th['border']}; border-radius: 8px; }}
        .footer-card h4 {{ margin: 0 0 0.5rem; }}
        .footer-card p {{ color: {th['text_muted']}; margin: 0; font-size: 0.85rem; }}
        .footer-card a {{ color: {th['accent']}; }}
        .cite-box {{
            padding: 1rem; border-left: 3px solid {th['accent']}; background: {th['bg_card']};
            font-size: 0.85rem; overflow-wrap: anywhere;
        }}
        .doi-badge {{ font-size: 0.8rem; margin-left: 0.5rem; color: {th['accent']}; }}
        #MainMenu, footer {{ visibility: hidden; }}
        header[data-testid="stHeader"] {{ background: transparent; }}
        button:focus-visible, input:focus-visible, textarea:focus-visible {{
            outline: 3px solid {th['accent']}80 !important; outline-offset: 2px !important;
        }}
        @media (max-width: 900px) {{
            [data-testid="stMainBlockContainer"] {{ padding: 1rem 0.85rem 2rem; }}
            .hero-shell {{ align-items: flex-start; flex-wrap: wrap; gap: 0.55rem; }}
            .context-scope {{ width: 100%; margin-left: 0; }}
        }}
        @media (prefers-reduced-motion: reduce) {{
            *, *::before, *::after {{ scroll-behavior: auto !important; transition: none !important; }}
        }}
    </style>
    """, unsafe_allow_html=True)

# =============================================================================
# 数据加载
# =============================================================================
@st.cache_resource(show_spinner=False)
def download_from_huggingface(repo_id: str, filename: str, revision: str = HF_REVISION):
    try:
        from huggingface_hub import hf_hub_download
    except ImportError:
        return None, False, "huggingface_hub not installed"
    try:
        path = hf_hub_download(repo_id=repo_id, filename=filename, repo_type="dataset", revision=revision)
        digest = hashlib.sha256()
        with open(path, 'rb') as source:
            for chunk in iter(lambda: source.read(1024 * 1024), b''):
                digest.update(chunk)
        df = read_score_csv(path)
        df.attrs['sha256'] = digest.hexdigest()
        return df, True, None
    except Exception as e:
        return None, False, f"HF error: {str(e)}"


@st.cache_resource(show_spinner=False, max_entries=3)
def load_uploaded_data(file_hash: str, _file_content: bytes):
    """按内容摘要缓存上传数据；返回的数据框在应用中只读使用。"""
    if not file_hash:
        raise ValueError("Missing upload content hash")
    return read_score_csv(_file_content)


def extract_gene_name(col_name: str) -> str:
    match = re.match(r'^([A-Za-z0-9_.-]+)\s*\(', str(col_name))
    if match:
        return match.group(1)
    return str(col_name)


@st.cache_data(show_spinner=False, max_entries=4)
def compute_gene_rankings(df_hash: str, _df: pd.DataFrame, gene_columns=None):
    return analyze_dataset(_df, gene_columns=gene_columns)


# =============================================================================
# 基因 × 药物 相关分析模块
# =============================================================================
@st.cache_resource(show_spinner=False)
def load_corr_datasets(repo_id, gdsc_file, crispr_file, compounds_file, model_file=None):
    """加载并缓存相关分析数据集。返回 (gdsc, crispr, compounds, model, err)。
    gdsc/crispr 的第一列是 ModelID，设为 index。
    model 为 ModelID->lineage 的精简表；若文件缺失则为 None（不报错，仅禁用分层）。"""
    try:
        from huggingface_hub import hf_hub_download
    except ImportError:
        return None, None, None, None, "huggingface_hub not installed"
    try:
        gdsc = pd.read_csv(
            hf_hub_download(repo_id=repo_id, filename=gdsc_file, repo_type="dataset"),
            index_col=0)
        crispr = pd.read_csv(
            hf_hub_download(repo_id=repo_id, filename=crispr_file, repo_type="dataset"),
            index_col=0)
        compounds = pd.read_csv(
            hf_hub_download(repo_id=repo_id, filename=compounds_file, repo_type="dataset"))
        model = None
        if model_file:
            try:
                m = pd.read_csv(
                    hf_hub_download(repo_id=repo_id, filename=model_file,
                                    repo_type="dataset"))
                keep = ['ModelID'] + [c for c in [LINEAGE_COL_COARSE, LINEAGE_COL_FINE]
                                      if c in m.columns]
                if 'ModelID' in m.columns and len(keep) > 1:
                    model = m[keep].set_index('ModelID')
            except Exception:
                model = None  # Model.csv 缺失或格式不符：分层功能禁用，主分析照常
        return gdsc, crispr, compounds, model, None
    except Exception as e:
        return None, None, None, None, f"HF error: {str(e)}"


def find_crispr_gene_column(crispr_df, gene_name):
    """在 'GENE (entrez)' 形式的列里按基因名(忽略大小写)定位列。"""
    target = gene_name.strip().upper()
    for col in crispr_df.columns:
        if extract_gene_name(col).upper() == target:
            return col
    return None


# 药名别名：搜常用名也能命中 GDSC2 里的旧代号
COMPOUND_ALIASES = {
    'rucaparib': ['ag-014699', 'ag014699', 'pf-01367338'],
    'veliparib': ['abt 888', 'abt-888', 'abt888'],
}


def search_compounds(compounds_df, query, gdsc2_only=False, gdsc_drug_ids=None):
    """按 CompoundName / CompoundID / GeneSymbolOfTargets / TargetOrMechanism / Synonyms 模糊匹配。
    支持逗号或空格分隔的多关键词（任一命中即收录）。
    gdsc2_only=True 时，仅保留在 GDSC2 矩阵中真有列的化合物（需传 gdsc_drug_ids）。
    返回按药名排序、按 CompoundID 去重的 [(CompoundID, label)]。"""
    raw = query.strip().lower()
    if not raw:
        return []
    if ',' in raw:
        terms = [s.strip() for s in raw.split(',') if s.strip()]
    else:
        terms = [s for s in raw.split() if s]
    if not terms:
        return []
    # 展开别名：搜 rucaparib 也匹配 ag-014699 等
    expanded = set(terms)
    for term in terms:
        for canon, aliases in COMPOUND_ALIASES.items():
            if term == canon or term in aliases:
                expanded.add(canon)
                expanded.update(aliases)
    terms = list(expanded)

    id_col = 'CompoundID' if 'CompoundID' in compounds_df.columns else None
    if id_col is None:
        return []
    name_col = 'CompoundName' if 'CompoundName' in compounds_df.columns else None
    tgt_cols = [c for c in ['GeneSymbolOfTargets', 'TargetOrMechanism', 'Synonyms']
                if c in compounds_df.columns]
    has_sid = 'SampleIDs' in compounds_df.columns

    seen = set()
    results = []
    for _, row in compounds_df.iterrows():
        cid = str(row[id_col])
        if not cid.startswith('DPC') or cid in seen:
            continue
        # GDSC2-only 过滤：优先用矩阵真实列，否则退回 SampleIDs 标注
        if gdsc2_only:
            if gdsc_drug_ids is not None:
                if cid not in gdsc_drug_ids:
                    continue
            elif has_sid and 'GDSC2' not in str(row['SampleIDs']):
                continue
        hay = cid.lower()
        if name_col:
            hay += " " + str(row[name_col]).lower()
        for tc in tgt_cols:
            hay += " " + str(row[tc]).lower()
        if any(term in hay for term in terms):
            nm = str(row[name_col]) if name_col else cid
            tgt = (str(row['GeneSymbolOfTargets'])
                   if 'GeneSymbolOfTargets' in compounds_df.columns else '')
            label = f"{nm} ({cid})" + (f" · {tgt}" if tgt and tgt != 'nan' else "")
            results.append((cid, label, nm))
            seen.add(cid)
    results.sort(key=lambda x: x[2].lower())
    return [(cid, label) for cid, label, _ in results]


def _rankdata(a):
    """平均秩（处理并列），纯 numpy，等价于 scipy.stats.rankdata 默认行为。"""
    a = np.asarray(a, dtype=float)
    order = a.argsort()
    ranks = np.empty(len(a), dtype=float)
    ranks[order] = np.arange(1, len(a) + 1)
    # 处理并列：同值取平均秩
    _, inv, counts = np.unique(a, return_inverse=True, return_counts=True)
    # 累计起始秩
    csum = np.cumsum(counts)
    start = csum - counts
    avg = (start + csum + 1) / 2.0  # 每组的平均秩 (1-indexed)
    return avg[inv]


def _spearman_with_p(x, y):
    """不依赖 scipy 的 Spearman ρ 与近似 p 值（t 分布双尾）。
    ρ = 秩变换后的 Pearson 相关；p 用学生 t 近似。"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    if n < 3:
        return float('nan'), float('nan')
    rx = _rankdata(x)
    ry = _rankdata(y)
    # 秩的 Pearson 相关
    sx = rx.std()
    sy = ry.std()
    if sx == 0 or sy == 0:
        return float('nan'), float('nan')
    rho = float(np.corrcoef(rx, ry)[0, 1])
    # 完全相关时 p->0
    if abs(rho) >= 1.0:
        return float(rho), 0.0
    # t = rho * sqrt((n-2)/(1-rho^2)) ~ t(n-2)
    t = rho * np.sqrt((n - 2) / (1 - rho ** 2))
    from math import lgamma, log
    df = n - 2
    def betacf(a, b, xx):
        MAXIT, EPS, FPMIN = 200, 3e-12, 1e-300
        qab, qap, qam = a + b, a + 1.0, a - 1.0
        c = 1.0; d = 1.0 - qab * xx / qap
        if abs(d) < FPMIN: d = FPMIN
        d = 1.0 / d; h = d
        for m in range(1, MAXIT + 1):
            m2 = 2 * m
            aa = m * (b - m) * xx / ((qam + m2) * (a + m2))
            d = 1.0 + aa * d
            if abs(d) < FPMIN: d = FPMIN
            c = 1.0 + aa / c
            if abs(c) < FPMIN: c = FPMIN
            d = 1.0 / d; h *= d * c
            aa = -(a + m) * (qab + m) * xx / ((a + m2) * (qap + m2))
            d = 1.0 + aa * d
            if abs(d) < FPMIN: d = FPMIN
            c = 1.0 + aa / c
            if abs(c) < FPMIN: c = FPMIN
            d = 1.0 / d; delta = d * c; h *= delta
            if abs(delta - 1.0) < EPS: break
        return h
    def betai(a, b, xx):
        if xx <= 0.0: return 0.0
        if xx >= 1.0: return 1.0
        lbeta = lgamma(a + b) - lgamma(a) - lgamma(b)
        bt = np.exp(lbeta + a * log(xx) + b * log(1.0 - xx))
        if xx < (a + 1.0) / (a + b + 2.0):
            return bt * betacf(a, b, xx) / a
        else:
            return 1.0 - bt * betacf(b, a, 1.0 - xx) / b
    x_beta = df / (df + t * t)
    p = betai(df / 2.0, 0.5, x_beta)
    return float(rho), float(min(max(p, 0.0), 1.0))


def compute_gene_drug_correlation(crispr_df, gdsc_df, gene_col, compound_id):
    """按 ModelID 对齐，返回 (merged_df, rho, p, n, err)。
    merged_df 列: dep (Gene Effect), auc (GDSC2 AUC)。"""
    if compound_id not in gdsc_df.columns:
        return None, None, None, 0, 'drug_not_found'
    dep = crispr_df[gene_col].rename('dep')
    auc = gdsc_df[compound_id].rename('auc')
    merged = pd.concat([dep, auc], axis=1).dropna()
    n = len(merged)
    if n < 10:
        return merged, None, None, n, 'too_few'
    rho, p = _spearman_with_p(merged['dep'].values, merged['auc'].values)
    return merged, float(rho), float(p), n, None


def create_correlation_scatter(merged, gene_name, drug_label, rho, p, point_size=4):
    """散点 + 线性拟合线，沿用 Morandi/theme_classic 风格，p 值斜体。"""
    th = get_theme()
    x = merged['dep'].values
    y = merged['auc'].values

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=x, y=y, mode='markers',
        marker=dict(size=point_size * 2.0, color=PLOT_COLORS['interest'],
                    opacity=0.55, line=dict(width=0.35, color=th['plot_axis'])),
        hovertemplate='Gene Effect: %{x:.3f}<br>AUC: %{y:.3f}<extra></extra>',
        showlegend=False,
    ))
    # 线性拟合线
    if len(x) >= 2:
        coef = np.polyfit(x, y, 1)
        xs = np.array([x.min(), x.max()])
        ys = coef[0] * xs + coef[1]
        fig.add_trace(go.Scatter(
            x=xs, y=ys, mode='lines',
            line=dict(color=PLOT_COLORS['essential'], width=1.2, dash='solid'),
            showlegend=False, hoverinfo='skip',
        ))

    p_txt = f"<i>p</i> = {p:.2g}" if p is not None else ""
    rho_txt = f"Spearman ρ = {rho:.2f}" if rho is not None else ""
    subtitle = f"{rho_txt}　{p_txt}　n = {len(merged)}"

    fig.update_layout(
        title=dict(text=f"{gene_name}  ×  {drug_label}<br><sub>{subtitle}</sub>",
                   font=dict(size=16, family=FONT_FAMILY), x=0.5, xanchor='center'),
        xaxis_title=t('corr_axis_x'),
        yaxis_title=t('corr_axis_y'),
        height=520,
        margin=dict(l=70, r=40, t=80, b=60),
    )
    fig.update_xaxes(showgrid=False, zeroline=False, linewidth=0.35, ticks='outside')
    fig.update_yaxes(showgrid=False, zeroline=False, linewidth=0.35, ticks='outside')
    apply_theme_to_fig(fig)
    return fig


# ---- 肿瘤类型分层 ----
# Morandi 风格分类配色（低饱和度）
MORANDI_PALETTE = [
    '#A7B5A0', '#C5A9A0', '#9FA8B0', '#C7B79B', '#B0A1B5',
    '#8FA8A3', '#C2A38F', '#A2A9B8', '#BFAFA0', '#A8B0A0',
    '#B5A0A8', '#9DAEB0', '#C8B5A0', '#A0A5B0', '#B8A8A0',
    '#A5B0A8', '#B0A0A5', '#A0B0B5', '#C0B0A0', '#ABA5B0',
]


def attach_lineage(merged, model_df, lineage_col):
    """给 merged(index=ModelID) 加一列 lineage；缺失或无 model_df 时返回原表+None。"""
    if model_df is None or lineage_col not in (model_df.columns if model_df is not None else []):
        return merged, None
    out = merged.copy()
    out['lineage'] = model_df[lineage_col].reindex(out.index)
    out['lineage'] = out['lineage'].fillna('Unknown')
    return out, 'lineage'


def stratified_correlation(merged_with_lin, min_n=MIN_N_PER_GROUP):
    """对每个 lineage 分别算 Spearman ρ/p；返回按 |ρ| 排序的 DataFrame，含 BH-FDR。
    需要 merged_with_lin 含列 dep/auc/lineage。"""
    rows = []
    for lin, g in merged_with_lin.groupby('lineage'):
        gg = g[['dep', 'auc']].dropna()
        if len(gg) < min_n:
            continue
        rho, p = _spearman_with_p(gg['dep'].values, gg['auc'].values)
        if np.isnan(rho):
            continue
        rows.append({'lineage': lin, 'n': len(gg), 'rho': rho, 'p': p})
    if not rows:
        return pd.DataFrame(columns=['lineage', 'n', 'rho', 'p', 'fdr'])
    df = pd.DataFrame(rows)
    # Benjamini-Hochberg
    df = df.sort_values('p').reset_index(drop=True)
    m = len(df)
    df['fdr'] = (df['p'] * m / (df.index + 1)).clip(upper=1.0)
    df['fdr'] = df['fdr'][::-1].cummin()[::-1]  # 单调化
    df = df.sort_values('rho', ascending=False).reset_index(drop=True)
    return df


def create_lineage_scatter(merged_lin, gene_name, drug_label, rho, p, point_size=4):
    """按 lineage 着色的总体散点（Morandi 调色），保留总体拟合线与总体 ρ。"""
    th = get_theme()
    fig = go.Figure()
    lineages = sorted(merged_lin['lineage'].dropna().unique().tolist())
    for i, lin in enumerate(lineages):
        sub = merged_lin[merged_lin['lineage'] == lin]
        fig.add_trace(go.Scatter(
            x=sub['dep'].values, y=sub['auc'].values, mode='markers', name=str(lin),
            marker=dict(size=point_size * 2.0,
                        color=MORANDI_PALETTE[i % len(MORANDI_PALETTE)],
                        opacity=0.7, line=dict(width=0.35, color=th['plot_axis'])),
            hovertemplate=f'{lin}<br>Gene Effect: %{{x:.3f}}<br>AUC: %{{y:.3f}}<extra></extra>',
        ))
    # 总体拟合线
    x = merged_lin['dep'].values
    y = merged_lin['auc'].values
    if len(x) >= 2:
        coef = np.polyfit(x, y, 1)
        xs = np.array([x.min(), x.max()])
        fig.add_trace(go.Scatter(
            x=xs, y=coef[0] * xs + coef[1], mode='lines', name='Overall fit',
            line=dict(color=th['plot_axis'], width=1.2),
            hoverinfo='skip', showlegend=False,
        ))
    p_txt = f"<i>p</i> = {p:.2g}" if p is not None else ""
    subtitle = f"Overall Spearman ρ = {rho:.2f}　{p_txt}　n = {len(merged_lin)}"
    fig.update_layout(
        title=dict(text=f"{gene_name}  ×  {drug_label}<br><sub>{subtitle}</sub>",
                   font=dict(size=16, family=FONT_FAMILY), x=0.5, xanchor='center'),
        xaxis_title=t('corr_axis_x'), yaxis_title=t('corr_axis_y'),
        height=560, margin=dict(l=70, r=40, t=80, b=60),
        legend=dict(font=dict(size=10, family=FONT_FAMILY)),
    )
    fig.update_xaxes(showgrid=False, zeroline=False, linewidth=0.35, ticks='outside')
    fig.update_yaxes(showgrid=False, zeroline=False, linewidth=0.35, ticks='outside')
    apply_theme_to_fig(fig)
    return fig


def create_forest_plot(strat_df, gene_name, drug_label):
    """分层相关森林图：每个 lineage 一行，点=ρ，按显著性着色。"""
    th = get_theme()
    d = strat_df.iloc[::-1].reset_index(drop=True)  # 从下往上画
    colors = [PLOT_COLORS['essential'] if r < 0.05 else th['plot_axis']
              for r in d['fdr']]
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=d['rho'], y=d['lineage'], mode='markers',
        marker=dict(size=[8 + min(np.sqrt(n), 12) for n in d['n']],
                    color=colors, opacity=0.85,
                    line=dict(width=0.35, color=th['plot_axis'])),
        customdata=np.stack([d['n'], d['p'], d['fdr']], axis=-1),
        hovertemplate=('%{y}<br>ρ = %{x:.2f}<br>n = %{customdata[0]}'
                       '<br>p = %{customdata[1]:.2g}<br>FDR = %{customdata[2]:.2g}<extra></extra>'),
        showlegend=False,
    ))
    fig.add_vline(x=0, line=dict(color=th['plot_axis'], width=0.5, dash='dot'))
    fig.update_layout(
        title=dict(text=f"Per-lineage Spearman ρ: {gene_name} × {drug_label}"
                        f"<br><sub>{t('corr_forest_sub')}</sub>",
                   font=dict(size=15, family=FONT_FAMILY), x=0.5, xanchor='center'),
        xaxis_title="Spearman ρ", yaxis_title="",
        height=max(320, 40 * len(d) + 120), margin=dict(l=160, r=40, t=80, b=50),
    )
    fig.update_xaxes(showgrid=False, zeroline=False, linewidth=0.35, ticks='outside')
    fig.update_yaxes(showgrid=False, linewidth=0.35,
                     tickfont=dict(size=10, family=FONT_FAMILY))
    apply_theme_to_fig(fig)
    return fig


# =============================================================================
# 绘图配置
# =============================================================================
PLOT_COLORS = {
    # Okabe–Ito 色盲友好配色
    'essential': '#D55E00',
    'nonessential': '#009E73',
    'interest': '#0072B2',
    'threshold': '#9BA6B0',
    'boxplot_fill': '#56B4E9',
}
FONT_FAMILY = "Inter, Helvetica Neue, Arial, sans-serif"

PLOT_CONFIG = {
    'displaylogo': False,
    'modeBarButtonsToRemove': [
        'zoom2d', 'pan2d', 'select2d', 'lasso2d',
        'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d',
        'hoverClosestCartesian', 'hoverCompareCartesian', 'toggleSpikelines'
    ],
    'toImageButtonOptions': {
        'format': 'png', 'filename': 'crispr_plot',
        'height': 600, 'width': 1000, 'scale': 3
    }
}


def centered_plot(fig, config=None):
    """使用容器宽度渲染，避免窄屏三列布局造成图表压缩。"""
    st.plotly_chart(fig, config=config or PLOT_CONFIG)


def apply_theme_to_fig(fig):
    th = get_theme()
    fig.update_layout(
        plot_bgcolor=th['plot_bg'],
        paper_bgcolor=th['plot_bg'],
        font=dict(family=FONT_FAMILY, size=13, color=th['plot_text']),
        hoverlabel=dict(font=dict(family=FONT_FAMILY, size=13)),
    )
    fig.update_xaxes(
        linecolor=th['plot_axis'], tickcolor=th['plot_axis'],
        tickfont=dict(size=12, color=th['plot_text']),
        title_font=dict(size=14, color=th['plot_text']),
    )
    fig.update_yaxes(
        linecolor=th['plot_axis'], tickcolor=th['plot_axis'],
        tickfont=dict(size=12, color=th['plot_text']),
        title_font=dict(size=14, color=th['plot_text']),
    )
    return fig


# =============================================================================
# 绘图函数
# =============================================================================
def _rank_label_annotations(points, x_range, y_range, background):
    """Space labels in a conservative viewport, keeping arrows on the true data.

    A small spatial grid bounds collision checks even for a large submitted set.
    Positions use data coordinates so labels and leader lines also remain vector
    objects in exports and scale with the responsive chart.
    """
    plot_width, plot_height = 720, 360
    x_span, y_span = x_range[1] - x_range[0], y_range[1] - y_range[0]
    occupied = {}
    annotations = []

    def cells(rect):
        left, top, right, bottom = rect
        return [(col, row) for col in range(int(left // 64), int(right // 64) + 1)
                for row in range(int(top // 64), int(bottom // 64) + 1)]

    offsets = [(dx, dy) for distance in (24, 48, 72, 96, 120)
               for dx in (0, 48, -48, 96, -96) for dy in (-distance, distance)]
    for point in sorted(points, key=lambda item: (item['x'], item['y'], item['gene'])):
        px = (point['x'] - x_range[0]) / x_span * plot_width
        py = (y_range[1] - point['y']) / y_span * plot_height
        half_width = max(24, len(point['gene']) * 3.7 + 7)
        half_height = 12
        best = None
        for dx, dy in offsets:
            cx = min(plot_width - half_width - 4, max(half_width + 4, px + dx))
            cy = min(plot_height - half_height - 4, max(half_height + 4, py + dy))
            rect = (cx - half_width, cy - half_height, cx + half_width, cy + half_height)
            nearby = {other for cell in cells(rect) for other in occupied.get(cell, [])}
            overlap = sum(max(0, min(rect[2], other[2]) - max(rect[0], other[0]))
                          * max(0, min(rect[3], other[3]) - max(rect[1], other[1]))
                          for other in nearby)
            if best is None or overlap < best[0]:
                best = (overlap, cx, cy, rect)
            if overlap == 0:
                break
        _, cx, cy, rect = best
        for cell in cells(rect):
            occupied.setdefault(cell, []).append(rect)
        annotations.append(dict(
            x=point['x'], y=point['y'], xref='x', yref='y',
            ax=x_range[0] + cx / plot_width * x_span,
            ay=y_range[1] - cy / plot_height * y_span, axref='x', ayref='y',
            text=escape(point['gene']), showarrow=True, arrowhead=0,
            arrowwidth=0.8, arrowcolor=point['color'], standoff=5,
            xanchor='center', yanchor='middle', borderpad=2, bgcolor=background,
            font=dict(size=13, color=point['color'], family=FONT_FAMILY),
        ))
    return annotations


def build_rank_figure(rankings, layers, references, n_cell_lines, point_size=4,
                      show_reference_labels=True):
    th = get_theme()
    fig = go.Figure()
    highlighted = {g.upper() for layer in layers for g in layer['genes']}
    refs = list(dict.fromkeys(g.upper() for g in references if g))
    bg = rankings[~rankings.gene_upper.isin(highlighted | set(refs))]
    hover = ('<b>%{text}</b><br>' + ui('Rank', '排名') + ': %{x:,}<br>'
             + ui('Mean score', '平均分') + ': %{y:.4f}<br>'
             + ui('Valid cell lines', '有效细胞系数') + ': %{customdata[0]}<br>'
             + ui('Rank percentile (lower = stronger)', '排名百分位（越低依赖越强）')
             + ': %{customdata[1]:.2f}%<extra></extra>')
    fig.add_trace(go.Scattergl(x=bg['rank'], y=bg.mean_score, mode='markers',
                              marker=dict(size=2.5, color=th['plot_scatter_bg']), name=ui('All genes', '全部基因'),
                              text=bg.gene, customdata=bg[['n_valid', 'percentile']].values, hovertemplate=hover))
    label_points = []

    def collect_labels(part, color):
        label_points.extend({'x': record.rank, 'y': record.mean_score,
                             'gene': record.gene, 'color': color}
                            for record in part.itertuples(index=False))

    for i, gene in enumerate(refs):
        if gene in highlighted:
            continue
        part = rankings[rankings.gene_upper == gene]
        if part.empty:
            continue
        color = PLOT_COLORS['essential'] if i == 0 else PLOT_COLORS['nonessential']
        fig.add_trace(go.Scatter(x=part['rank'], y=part.mean_score, mode='markers',
                                 marker=dict(size=point_size * 2.6, color=color,
                                             symbol='diamond' if i == 0 else 'square',
                                             line=dict(width=1, color=th['plot_bg'])),
                                 text=part.gene, name=ui(f'Reference · {gene}', f'参考基因 · {gene}'),
                                 customdata=part[['n_valid', 'percentile']].values, hovertemplate=hover))
        if show_reference_labels:
            collect_labels(part, color)
    for layer in layers:
        part = rankings[rankings.gene_upper.isin({g.upper() for g in layer['genes']})].sort_values('rank')
        if part.empty:
            continue
        fig.add_trace(go.Scatter(x=part['rank'], y=part.mean_score,
                                 mode='markers',
                                 marker=dict(size=point_size * layer.get('size', 2.5), color=layer['color'],
                                             symbol=layer.get('symbol', 'circle'),
                                             line=dict(width=1, color=th['plot_bg'])),
                                 text=part.gene,
                                 name=layer['name'],
                                 customdata=part[['n_valid', 'percentile']].values, hovertemplate=hover))
        if layer['labels']:
            collect_labels(part, layer['color'])
    fig.add_hline(y=ESSENTIALITY_THRESHOLD, line=dict(dash='dash', color=PLOT_COLORS['threshold'], width=1))
    fig.add_hline(y=0, line=dict(color=th['plot_reference'], width=0.7))
    lo, hi = rankings.mean_score.min(), rankings.mean_score.max()
    span = max(hi - lo, 0.2)
    x_range = [0, max(len(rankings) * 1.02, 1)]
    y_range = [lo - span * 0.12, hi + span * 0.18]
    fig.update_layout(
        xaxis=dict(title=ui('Gene rank', '基因排名'), tickformat=',d', range=x_range,
                   showgrid=False, automargin=True),
        yaxis=dict(title=ui('Mean CRISPR score', '平均 CRISPR 分数'), range=y_range,
                   showgrid=True, gridcolor=th['plot_grid'], automargin=True),
        annotations=_rank_label_annotations(label_points, x_range, y_range, th['plot_bg']),
        title=dict(text=ui(f'{n_cell_lines:,} cell lines · hover for gene-level sample counts',
                           f'{n_cell_lines:,} 个细胞系 · 悬停查看各基因有效样本数'),
                   font=dict(size=13), x=0, xanchor='left'),
        legend=dict(orientation='h', yanchor='top', y=-0.18, xanchor='center', x=0.5,
                    font=dict(size=12), itemsizing='constant'),
        height=540, margin=dict(l=75, r=30, t=50, b=100))
    return apply_theme_to_fig(fig)


def create_rank_plot(gene_rank_df, genes_of_interest, essential_gene='MYC', nonessential_gene='PTEN',
                     n_cell_lines=0, show_labels=True, point_size=4):
    return build_rank_figure(gene_rank_df,
                            [{'genes': genes_of_interest, 'labels': show_labels, 'color': PLOT_COLORS['interest'],
                              'name': ui('Genes of interest', '目标基因')}],
                            [essential_gene, nonessential_gene], n_cell_lines, point_size,
                            show_reference_labels=show_labels)


def create_multilayer_rank_plot(gene_rank_df, background_genes, highlight_genes, bg_color='#8DA8C0', hl_color='#0072B2',
                               essential_gene='MYC', nonessential_gene='PTEN', n_cell_lines=0, show_labels=True):
    bg_only = [g for g in background_genes if g not in highlight_genes]
    layers = [
        {'genes': bg_only, 'labels': False, 'color': bg_color, 'size': 1.7,
         'name': ui('Background gene set', '背景基因集')},
        {'genes': highlight_genes, 'labels': show_labels, 'color': hl_color, 'symbol': 'circle',
         'name': ui('Highlight genes', '高亮基因')},
    ]
    return build_rank_figure(gene_rank_df, layers, [essential_gene, nonessential_gene], n_cell_lines,
                            st.session_state.get('point_size', 4), show_reference_labels=show_labels)


def create_lineage_boxplot(lineage_data, genes, sort_mode='median', show_points=False):
    th = get_theme()
    # Horizontal groups keep long cancer-type names readable on narrow screens.
    genes = [g for g in genes if g in set(lineage_data.gene)]
    n_genes = len(genes)
    fig = make_subplots(rows=n_genes, cols=1, vertical_spacing=min(0.08, 0.35 / n_genes),
                        subplot_titles=genes)
    group_max = 1
    score_min, score_max = lineage_data.crispr_score.min(), lineage_data.crispr_score.max()
    pad = max((score_max - score_min) * 0.06, 0.15)
    for row, gene in enumerate(genes, 1):
        sub = lineage_data[lineage_data.gene == gene]
        summary = sub.groupby('lineage').crispr_score.agg(['median', 'count']).reset_index()
        summary = summary.sort_values(['median', 'lineage'] if sort_mode == 'median' else ['lineage'])
        group_max = max(group_max, len(summary))
        labels = []
        for record in summary.itertuples(index=False):
            label = f'{record.lineage} (n={record.count})'
            labels.append(label)
            points = sub[sub.lineage == record.lineage]
            fig.add_trace(go.Box(x=points.crispr_score, y=[label] * len(points), orientation='h',
                                 name=record.lineage, showlegend=False,
                                 marker=dict(color=PLOT_COLORS['interest'], size=4, opacity=0.5),
                                 line=dict(color=PLOT_COLORS['interest'], width=1.2),
                                 fillcolor='rgba(0,114,178,0.18)',
                                 boxpoints='all' if show_points else 'outliers', jitter=0.35, pointpos=0,
                                 customdata=points[['cell_line', 'cell_line_name']].values,
                                 hovertemplate=ui('Cell line', '细胞系') + ': %{customdata[1]} (%{customdata[0]})<br>'
                                               + ui('Score', '分数') + ': %{x:.4f}<extra>%{y}</extra>'), row=row, col=1)
        fig.update_yaxes(categoryorder='array', categoryarray=labels, autorange='reversed',
                         automargin=True, showgrid=False, row=row, col=1)
        fig.update_xaxes(title_text=ui('CRISPR score', 'CRISPR 分数'), range=[score_min - pad, score_max + pad],
                         showgrid=True, gridcolor=th['plot_grid'], row=row, col=1)
        fig.add_vline(x=0, line=dict(color=th['plot_reference'], width=0.7), row=row, col=1)
        fig.add_vline(x=ESSENTIALITY_THRESHOLD, line=dict(dash='dash', color=PLOT_COLORS['threshold'], width=1), row=row, col=1)
    fig.update_annotations(font=dict(size=15, family=FONT_FAMILY, color=th['plot_text']))
    fig.update_layout(height=max(300, group_max * 30 + 90) * n_genes + 50,
                      margin=dict(l=145, r=25, t=45, b=55), showlegend=False)
    return apply_theme_to_fig(fig)


# =============================================================================
# 图片导出（始终白底，方便论文用）
# =============================================================================
# =============================================================================
# Citation 渲染
# =============================================================================
def build_citations():
    has_doi = bool(ZENODO_DOI.strip())

    if has_doi:
        apa = (f"{TOOL_AUTHORS} ({TOOL_YEAR}). CRISPR Score Analyzer ({TOOL_VERSION}) "
               f"[Software]. Zenodo. https://doi.org/{ZENODO_DOI}")
        bibtex = f"""@software{{crispr_score_analyzer_{TOOL_YEAR},
  author  = {{{TOOL_AUTHORS}}},
  title   = {{CRISPR Score Analyzer: An interactive platform for DepMap gene essentiality}},
  year    = {{{TOOL_YEAR}}},
  version = {{{TOOL_VERSION}}},
  doi     = {{{ZENODO_DOI}}},
  url     = {{https://doi.org/{ZENODO_DOI}}}
}}"""
    else:
        apa = (f"{TOOL_AUTHORS} ({TOOL_YEAR}). CRISPR Score Analyzer ({TOOL_VERSION}) "
               f"[Software]. {GITHUB_URL}")
        bibtex = f"""@software{{crispr_score_analyzer_{TOOL_YEAR},
  author  = {{{TOOL_AUTHORS}}},
  title   = {{CRISPR Score Analyzer: An interactive platform for DepMap gene essentiality}},
  year    = {{{TOOL_YEAR}}},
  version = {{{TOOL_VERSION}}},
  url     = {{{GITHUB_URL}}}
}}"""

    return apa, bibtex, has_doi


def render_citation_section():
    apa, bibtex, has_doi = build_citations()

    st.markdown(f"### {t('citation')}")

    if has_doi:
        st.markdown(
            f'<a href="https://doi.org/{ZENODO_DOI}" target="_blank" class="doi-badge">'
            f'DOI: {ZENODO_DOI}</a>',
            unsafe_allow_html=True
        )
    else:
        st.caption(f"⏳ {t('doi_pending')}")

    st.markdown(f"**{t('cite_this_tool')} (APA):**")
    st.markdown(f'<div class="cite-box">{apa}</div>', unsafe_allow_html=True)

    st.markdown("**BibTeX:**")
    st.code(bibtex, language="bibtex")

    st.markdown("**DepMap data citation:**")
    st.markdown(
        '<div class="cite-box">Tsherniak, A., et al. (2017). '
        'Defining a Cancer Dependency Map. <i>Cell</i> 170, 564–576. '
        'https://depmap.org/portal/</div>',
        unsafe_allow_html=True
    )


def ui(en, zh):
    return zh if st.session_state.get('lang') == 'zh' else en


def render_download_buttons(fig, filename_base, key_prefix, height=600):
    """Keep at most one current file per format, tied to its exact figure."""
    with st.expander(ui('Figure size', '图形尺寸'), expanded=False):
        c1, c2 = st.columns(2)
        width = c1.number_input(ui('Width (px)', '宽度（px）'), 400, 3000, 1000, 50,
                                key=f'{key_prefix}_width')
        out_height = c2.number_input(ui('Height (px)', '高度（px）'), 300, 12000,
                                    min(12000, max(300, int(height))), 50,
                                    key=f'{key_prefix}_height')
    st.caption(ui(f'PDF / SVG: full vector. PNG: {width * 3} × {out_height * 3} px, 300 DPI.',
                  f'PDF / SVG：完整矢量。PNG：{width * 3} × {out_height * 3} 像素，300 DPI。'))
    for col, fmt, mime in zip(st.columns(3), ['pdf', 'png', 'svg'],
                              ['application/pdf', 'image/png', 'image/svg+xml']):
        cache_key = f'_export_{key_prefix}_{fmt}'
        signature = export_fingerprint(fig, fmt, width, out_height)
        cached = st.session_state.get(cache_key)
        if cached and cached['signature'] != signature:
            del st.session_state[cache_key]
            cached = None
        with col:
            if fmt == 'png' and width * out_height * 9 > MAX_PNG_PIXELS:
                st.caption(ui('PNG is too large. Reduce dimensions or use PDF/SVG.',
                              'PNG 尺寸过大，请缩小宽高或使用 PDF/SVG。'))
                continue
            if cached:
                st.download_button(f'↓ {fmt.upper()}', cached['bytes'],
                                   f'{filename_base}.{fmt}', mime,
                                   key=f'{key_prefix}_{fmt}_dl', on_click='ignore', width='stretch')
            elif st.button(ui(f'Generate {fmt.upper()}', f'生成 {fmt.upper()}'),
                           key=f'{key_prefix}_{fmt}_btn', width='stretch'):
                with st.spinner(ui('Generating figure…', '正在生成图片…')):
                    try:
                        image_bytes = render_figure_bytes(fig, fmt, width, out_height)
                        st.session_state[cache_key] = {'signature': signature, 'bytes': image_bytes}
                        st.rerun()
                    except Exception as exc:
                        st.error(ui('Export failed: ', '导出失败：') + str(exc))


def render_diagnostics():
    with st.expander(ui('Data checks & column diagnostics', '数据校验与列详情')):
        st.caption(ui('Only explicitly recognized gene columns enter the ranking. Constant genes and small cohorts are retained; nonfinite scores are excluded and counted. Means use each gene’s valid scores.',
                      '排名仅使用明确识别的基因列。保留低变异基因和小样本；无效分数单独计数。均值使用各基因的有效分数。'))
        st.dataframe(diagnostics, hide_index=True, width='stretch')
        st.download_button(ui('Column diagnostics CSV', '列校验 CSV'), diagnostics.to_csv(index=False),
                            'column_diagnostics.csv', 'text/csv', key='diagnostics_csv', on_click='ignore')


def read_gene_input(prefix, default, method):
    if method == 'text':
        value = st.text_area(t('gene_list'), default, height=90, key=f'{prefix}_genes',
                             help=ui('Separate with newlines, spaces, tabs, commas or semicolons.',
                                     '支持换行、空格、制表符、中英文逗号和分号。'))
        return value, None
    file = st.file_uploader(t('input_file'), type=['csv', 'txt', 'tsv'], key=f'{prefix}_file',
                            help=ui('First column; a gene/symbol header is optional.',
                                    '读取第一列；可有 gene/symbol 表头，也可无表头。'))
    return '', file


def parse_input(text, file):
    return parse_gene_upload(file.getvalue(), file.name) if file else parse_gene_text(text)


def show_matches(parsed, matched, missing, title=None):
    if title:
        st.markdown(f'**{title}**')
    st.caption(ui(
        f"Input {parsed['input_count']} · unique {len(parsed['genes'])} · matched {len(matched)} · not matched {len(missing)} · duplicates {len(parsed['duplicates'])}",
        f"输入 {parsed['input_count']} · 去重后 {len(parsed['genes'])} · 匹配 {len(matched)} · 未匹配 {len(missing)} · 重复 {len(parsed['duplicates'])}"))
    if matched:
        st.markdown(' '.join(f'<span class="gene-tag">{escape(g)}</span>' for g in matched),
                    unsafe_allow_html=True)
    if missing or parsed['duplicates']:
        with st.expander(ui('Input details / unmatched genes', '输入详情 / 未匹配基因')):
            if missing:
                st.write(ui('Not matched (check column diagnostics for excluded genes):',
                            '未匹配（被排除的基因请查看列校验详情）：'))
                st.code('\n'.join(missing), language=None)
            if parsed['duplicates']:
                st.write(ui('Duplicate entries removed:', '已去除的重复输入：'), ', '.join(parsed['duplicates']))


def build_result(genes, matches, config=None, plot_data=None, draft=None):
    return {
        'genes': genes,
        'draft': copy.deepcopy(draft or {}),
        'input_signature': input_signature(draft or {}),
        'rankings': gene_rankings,
        'matches': matches,
        'config': config or {},
        'plot_data': get_lineage_frame(df_scope, genes, gene_rankings, lineage_col) if plot_data is None else plot_data,
        'metadata': {**dataset_metadata, 'created_at_utc': datetime.now(timezone.utc).isoformat(),
                     'selected_genes': genes, 'parameters': config or {},
                     'input_matching': matches},
    }


def result_downloads(result, fig, prefix):
    selected = result['rankings'].set_index('gene', drop=False).reindex(result['genes']).reset_index(drop=True)
    title_col, figure_col, data_col = st.columns([2, 1, 1.15], vertical_alignment='center')
    with title_col:
        st.markdown(f"### {ui('Results', '分析结果')}")
    with figure_col:
        with st.popover(ui('Export figure', '导出图形'), width='stretch'):
            render_download_buttons(fig, prefix, prefix, height=int(fig.layout.height or 600))
    with data_col:
        st.download_button(ui('Download results', '下载结果'),
                           selected.drop(columns=['gene_upper'], errors='ignore').to_csv(index=False),
                           f'{prefix}_selected_results.csv', 'text/csv', key=f'{prefix}_csv',
                           help=ui('CSV with original field names and full numeric precision.',
                                   'CSV 保留规范字段名与原始数值精度。'),
                           on_click='ignore', width='stretch')
    display_columns = ['gene', 'mean_score', 'rank', 'percentile', 'n_valid']
    table_config = {
        'gene': st.column_config.TextColumn(ui('Gene', '基因')),
        'mean_score': st.column_config.NumberColumn(ui('Mean score', '平均分'), format='%.3f'),
        'rank': st.column_config.NumberColumn(ui('Rank', '排名'), format='%d'),
        'percentile': st.column_config.NumberColumn(ui('Rank percentile', '排名百分位'), format='%.2f%%'),
        'n_valid': st.column_config.NumberColumn(ui('Valid samples', '有效样本数'), format='%d'),
        'n_missing': st.column_config.NumberColumn(ui('Missing samples', '缺失样本数'), format='%d'),
        'missing_fraction': st.column_config.NumberColumn(ui('Missing fraction', '缺失比例'), format='percent'),
    }
    st.dataframe(selected[display_columns], hide_index=True, width='stretch', column_config=table_config)
    st.caption(ui('Lower rank percentile indicates stronger mean dependency. Valid samples count finite scores. CSV downloads retain the original field names and full precision.',
                  '排名百分位越低，平均依赖越强。有效样本数仅统计有限数值。CSV 下载保留规范字段名与原始精度。'))
    if prefix == 'box':
        st.caption(ui('These ranks and means use the sidebar cohort. Group statistics and plot data use the cancer types selected for this comparison.',
                      '上表排名与均值基于侧栏癌种范围；分组统计与绘图数据基于本次比较所选癌种。'))
    metadata = {**result['metadata'], 'display': {
        'language': st.session_state.lang, 'theme': st.session_state.theme,
        'reference_genes': [essential_gene, nonessential_gene],
        'show_labels': show_labels, 'point_size': point_size,
    }, 'export': {'width': st.session_state.get(f'{prefix}_width', 1000),
                  'height': st.session_state.get(f'{prefix}_height', int(fig.layout.height or 600)),
                  'png_scale': 3, 'png_dpi': 300, 'full_vector_pdf_svg': True}}
    with st.expander(ui('Details & reproducibility', '详细数据与复现'), expanded=False):
        st.markdown(f"**{ui('Sample completeness', '样本完整性')}**")
        st.dataframe(selected[['gene', 'n_valid', 'n_missing', 'missing_fraction']],
                     hide_index=True, width='stretch', column_config=table_config)
        if prefix == 'box':
            st.markdown(f"**{ui('Cancer-type summaries', '癌种分组统计')}**")
            group_config = {
                'gene': table_config['gene'],
                'lineage': st.column_config.TextColumn(ui('Cancer type', '癌种')),
                'n_valid': table_config['n_valid'],
                **{key: st.column_config.NumberColumn(label, format='%.3f') for key, label in [
                    ('median', ui('Median', '中位数')), ('mean', ui('Mean', '均值')),
                    ('std', ui('SD', '标准差')), ('min', ui('Minimum', '最小值')),
                    ('max', ui('Maximum', '最大值')),
                ]},
            }
            st.dataframe(lineage_summary(result['plot_data']), hide_index=True, width='stretch',
                         column_config=group_config)
        st.download_button(ui('Plot data CSV', '绘图数据 CSV'), result['plot_data'].to_csv(index=False),
                           f'{prefix}_plot_data.csv', 'text/csv', key=f'{prefix}_data', on_click='ignore')
        bundle_key = f'_export_{prefix}_zip'
        bundle_signature = hashlib.sha256(json.dumps(metadata, sort_keys=True, default=str).encode()).hexdigest()
        bundle = st.session_state.get(bundle_key)
        if bundle and bundle['signature'] != bundle_signature:
            del st.session_state[bundle_key]
            bundle = None
        if bundle is None and st.button(ui('Prepare analysis ZIP', '准备完整分析 ZIP'), key=f'{prefix}_zip_btn'):
            with st.spinner(ui('Preparing analysis bundle…', '正在准备分析文件…')):
                bundle = {'signature': bundle_signature,
                          'bytes': make_analysis_bundle(result['rankings'], result['genes'], result['plot_data'],
                                                        metadata, diagnostics)}
                st.session_state[bundle_key] = bundle
        if bundle:
            st.download_button(ui('Complete analysis ZIP', '完整分析 ZIP'), bundle['bytes'],
                               f'{prefix}_analysis.zip', 'application/zip', key=f'{prefix}_zip', on_click='ignore')
        st.caption(ui('ZIP includes selected results, complete cohort ranking, plot data, input matching, column diagnostics and versioned analysis parameters. Download figures separately.',
                      'ZIP 包含所选结果、当前范围完整排名、绘图数据、输入匹配情况、列校验与版本参数。图形请单独下载。'))
        st.json(metadata, expanded=False)


def input_draft(text='', file=None, **options):
    return {'text': text, 'file_name': file.name if file else None,
            'file_sha256': hashlib.sha256(file.getvalue()).hexdigest() if file else None,
            **options}


def input_signature(draft):
    return hashlib.sha256(json.dumps(draft, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def apply_input_preset(prefix, clear=False):
    presets = {'rank': {'genes': 'E2F1\nE2F2\nE2F3'}, 'box': {'genes': 'E2F1\nE2F2'},
               'multi': {'bg': 'CDK1\nCDK2\nCCNB1\nCCND1\nCCNE1', 'hl': 'PLK1\nAURKA'}}
    for key, value in presets[prefix].items():
        st.session_state[f'{prefix}_{key}'] = '' if clear else value
    if prefix != 'multi':
        st.session_state[f'{prefix}_method'] = 'text'
    if prefix == 'box':
        st.session_state.pop('box_candidates', None)


def input_actions(prefix):
    run_col, example_col, clear_col = st.columns([2, 1, 1])
    with run_col:
        submitted = st.button(ui('Run analysis', 'Run · 生成图表'), type='primary',
                              key=f'{prefix}_run', width='stretch')
    with example_col:
        st.button(ui('Load example', '载入示例'), key=f'{prefix}_example', width='stretch',
                  on_click=apply_input_preset, args=(prefix,))
    with clear_col:
        st.button(ui('Clear', '清空输入'), key=f'{prefix}_clear', width='stretch',
                  on_click=apply_input_preset, args=(prefix, True))
    return submitted


def render_run_status(prefix, draft, result):
    if result is None:
        st.info(ui('Ready to analyze. Enter genes or load an example, then Run.',
                   '准备就绪：输入基因或载入示例，点击 Run 生成结果。'))
        return
    if result.get('input_signature') != input_signature(draft):
        st.warning(ui('Inputs changed. Run to update the result.', '输入已修改，点击 Run 更新结果。'))
    meta = result['metadata']
    scope = ', '.join(meta['cohort_lineages'])
    if scope == 'ALL':
        scope = ui('All cancer types', '全部癌种')
    timestamp = datetime.fromisoformat(meta['created_at_utc']).strftime('%Y-%m-%d %H:%M:%S UTC')
    summary = ui(f"{len(result['genes'])} genes · {meta['cohort_rows']:,} cell lines · {scope}",
                 f"{len(result['genes'])} 个基因 · {meta['cohort_rows']:,} 个细胞系 · {scope}")
    st.markdown(f'<div class="status-strip"><strong>{ui("Displayed result", "当前显示结果")}</strong>'
                f' · {escape(summary)}<br>{escape(timestamp)}</div>', unsafe_allow_html=True)


def update_preferences(preference):
    options = TRANSLATIONS if preference == 'lang' else THEMES
    value = st.session_state.get(f'{preference}_select')
    if value in options:
        st.session_state[preference] = value


def render_preferences():
    with st.sidebar.expander(ui('Preferences', '界面设置')):
        st.selectbox('Language / 语言', ['en', 'zh'], key='lang_select',
                     index=['en', 'zh'].index(st.session_state.lang),
                     format_func={'en': 'English', 'zh': '中文'}.get,
                     on_change=update_preferences, args=('lang',))
        st.selectbox(t('theme'), ['light', 'dark'], key='theme_select',
                     index=['light', 'dark'].index(st.session_state.theme),
                     format_func={x: t(x) for x in ['light', 'dark']}.get,
                     on_change=update_preferences, args=('theme',))


for draft_key in ['rank_genes', 'box_genes', 'multi_bg', 'multi_hl', 'rank_method', 'box_method',
                  'multi_bg_color', 'multi_hl_color', 'box_lineages', 'box_order', 'box_points']:
    if draft_key in st.session_state:
        st.session_state[draft_key] = st.session_state[draft_key]

for preference, options in [('lang', TRANSLATIONS), ('theme', THEMES)]:
    if st.session_state.get(f'{preference}_select') not in options:
        st.session_state[f'{preference}_select'] = st.session_state[preference]
inject_css()
with st.sidebar:
    st.markdown(f"## {t('data_source')}")
    source_mode = st.segmented_control(ui('Dataset', '数据集'), ['builtin', 'custom'], default='builtin',
                                       format_func={'builtin': ui('Built-in', '内置数据'),
                                                    'custom': ui('Custom CSV', '自定义 CSV')}.get,
                                       key='dataset_source', selection_mode='single', width='stretch',
                                       label_visibility='collapsed') or 'builtin'
    uploaded_file = None
    if source_mode == 'custom':
        uploaded_file = st.file_uploader(t('upload_csv'), type=['csv'], key='dataset_file')
    else:
        st.caption(f'{DATA_VERSION} · {SCORE_TYPE}')

crispr_data = None
source_sha = None
source_id = ('upload:' + hashlib.sha256(uploaded_file.getvalue()).hexdigest()
             if uploaded_file is not None else ('hf:' + HF_REVISION if source_mode == 'builtin' else 'custom:empty'))
if st.session_state.get('_source_id') != source_id:
    for key in list(st.session_state):
        if key in {'rank_result', 'box_result', 'multi_result', 'box_candidates', 'box_display_genes', 'gene_columns', 'custom_mapping', 'box_lineages'} or key.startswith('_export_'):
            del st.session_state[key]
    st.session_state['_source_id'] = source_id
if source_mode == 'custom' and uploaded_file is None:
    render_preferences()
    st.markdown(f'# {t("app_title")}')
    st.info(ui('Upload a CRISPR score CSV in the sidebar to begin, or select Built-in to explore the reference dataset.',
               '请在侧栏上传 CRISPR 分数 CSV，或选择“内置数据”开始探索。'))
    st.stop()
try:
    if uploaded_file is not None:
        uploaded_bytes = uploaded_file.getvalue()
        source_sha = hashlib.sha256(uploaded_bytes).hexdigest()
        crispr_data = load_uploaded_data(source_sha, uploaded_bytes)
    else:
        with st.spinner(t('loading_hf')):
            crispr_data, success, err = download_from_huggingface(HF_REPO_ID, HF_FILENAME)
        if not success:
            st.error(err)
    if crispr_data is None or crispr_data.empty:
        render_preferences()
        st.error(ui('No data rows. Upload a nonempty score matrix.', '数据没有有效行，请上传非空分数矩阵。'))
        st.stop()
except (ValueError, UnicodeError, pd.errors.ParserError, pd.errors.EmptyDataError) as exc:
    render_preferences()
    st.error(ui('Cannot read this CSV: ', '无法读取 CSV：') + str(exc))
    st.stop()

df = crispr_data
source_sha = source_sha or df.attrs.get('sha256')
lineage_col = find_lineage_column(df)
gene_columns = None
builtin_metadata = ['depmap_id', 'cell_line_display_name', 'lineage_1', 'lineage_2', 'lineage_3', 'lineage_6', 'lineage_4']
if uploaded_file is None and list(df.columns[:7]) == builtin_metadata:
    gene_columns = df.columns[7:].tolist()
with st.sidebar:
    if uploaded_file is not None:
        with st.expander(ui('Gene column mapping', '基因列映射')):
            st.caption(ui('Automatically recognize SYMBOL (ENTREZ_ID). For plain symbols, explicitly select score columns. Metadata is never inferred from score ranges.',
                          '自动识别 SYMBOL (ENTREZ_ID)。仅有基因名时，请手动选择分数列；不会依据数值范围把元数据当作基因。'))
            if st.checkbox(ui('Select columns manually', '手动选择基因列'), key='custom_mapping'):
                gene_columns = st.multiselect(ui('Score columns', '分数列'), list(dict.fromkeys(df.columns)), key='gene_columns')
    st.markdown(ui('### Analysis scope', '### 分析范围'))
    if lineage_col:
        labels = df[lineage_col].fillna('Unknown').astype(str).str.strip().replace('', 'Unknown')
        all_lineages = sorted(labels.unique().tolist())
        if any(x not in all_lineages for x in st.session_state.get('cohort', [])):
            st.session_state.cohort = []
        cohort = st.multiselect(ui('Cancer types', '癌种'), all_lineages, key='cohort',
                                placeholder=ui('All cancer types', '全部癌种'),
                                help=ui('Leave empty to include all cancer types. Changing the scope clears previous results.',
                                        '留空表示全部癌种；修改分析范围后，需要重新运行。'))
        df_scope = df.loc[labels.isin(cohort)] if cohort else df
    else:
        cohort, all_lineages, df_scope = [], [], df
        st.caption(t('lineage_missing'))
    with st.expander(ui('Figure settings', '图形设置'), expanded=False):
        essential_gene = st.text_input(t('essential'), value='MYC', key='reference_a').strip().upper()
        nonessential_gene = st.text_input(t('nonessential'), value='PTEN', key='reference_b').strip().upper()
        show_labels = st.checkbox(t('show_labels'), value=True, key='show_labels')
        point_size = st.slider(t('point_size'), 2, 8, 4, key='point_size')

context = json.dumps({'sha256': source_sha, 'cohort': sorted(cohort),
                      'gene_columns': gene_columns, 'source': uploaded_file.name if uploaded_file else HF_FILENAME}, sort_keys=True)
context_hash = hashlib.sha256(context.encode()).hexdigest()
if st.session_state.get('_analysis_context') != context_hash:
    for key in list(st.session_state):
        if key in {'rank_result', 'box_result', 'multi_result', 'box_candidates', 'box_display_genes', 'box_lineages'} or key.startswith('_export_'):
            del st.session_state[key]
    st.session_state['_analysis_context'] = context_hash

gene_rankings, diagnostics = compute_gene_rankings(context_hash, df_scope, gene_columns)
n_cell_lines = len(df_scope)
hero_version = t('custom_dataset') if uploaded_file is not None else DATA_VERSION
st.markdown(f'<section class="hero-shell"><div><h1 class="main-header">{t("app_title")}</h1>'
            f'<p class="sub-header">{t("app_subtitle")}</p></div>'
            f'<div class="hero-version">{escape(hero_version)}</div></section>', unsafe_allow_html=True)
scope_name = ', '.join(cohort) if cohort else ui('All cancer types', '全部癌种')
if gene_rankings.empty:
    render_preferences()
    render_diagnostics()
    st.warning(ui('No analyzable gene columns. Check the column diagnostics or manually map your score columns in the sidebar.',
                  '没有可分析的基因列。请检查列详情，或在侧栏手动映射分数列。'))
    st.stop()

st.markdown(f'<div class="context-strip"><span><strong>{n_cell_lines:,}</strong> {t("cell_lines")}</span>'
            f'<span><strong>{len(gene_rankings):,}</strong> {t("gene_count")}</span>'
            f'<span class="context-scope">{escape(scope_name)}</span></div>', unsafe_allow_html=True)
with st.sidebar.expander(ui('Data overview', '数据概况')):
    st.metric(t('essential_genes').format(threshold=ESSENTIALITY_THRESHOLD),
              f"{(gene_rankings.mean_score < ESSENTIALITY_THRESHOLD).sum():,}")
    st.caption(ui('Mean-score range: ', '平均分范围：') +
               f'{gene_rankings.mean_score.min():.2f} – {gene_rankings.mean_score.max():.2f}')
    st.caption(ui('Lower score means stronger dependency. Rankings use means within the analysis scope; they do not establish cancer-type selectivity. The −0.5 line is a screening reference.',
                  '分数越低，依赖越强。排名使用分析范围内的平均分，本身不能证明癌种特异性。−0.5 是筛选参考线。'))
render_preferences()
dataset_metadata = {
    'data_source': uploaded_file.name if uploaded_file else HF_REPO_ID,
    'data_file': uploaded_file.name if uploaded_file else HF_FILENAME,
    'data_release': 'custom' if uploaded_file else DATA_VERSION,
    'data_revision': None if uploaded_file else HF_REVISION,
    'data_sha256': source_sha, 'cohort_lineages': cohort or ['ALL'],
    'lineage_column': lineage_col, 'gene_column_mapping': gene_columns,
    'cohort_rows': n_cell_lines, 'ranked_genes': len(gene_rankings),
    'mean_score_cutoff': ESSENTIALITY_THRESHOLD,
    'software': {name: version(name) for name in ['streamlit', 'pandas', 'numpy', 'plotly', 'kaleido', 'pillow']},
    'analysis_schema': 2,
}

view = st.segmented_control(ui('Analysis', '分析'), ['rank', 'box', 'multi'], default='rank',
                            format_func={'rank': t('tab1'), 'box': t('tab2'), 'multi': t('tab3')}.get,
                            key='analysis_tab', selection_mode='single', label_visibility='collapsed', width='stretch') or 'rank'
if view == 'rank':
    st.caption(ui('Locate genes in the dependency ranking. Lower mean scores indicate stronger dependency.',
                  '定位关注基因的依赖排名；平均分越低，依赖越强。'))
    with st.container(border=True):
        method = st.radio(t('input_method'), ['text', 'file'], horizontal=True, key='rank_method',
                          format_func={'text': t('input_direct'), 'file': t('input_file')}.get,
                          label_visibility='collapsed')
        gene_text, gene_file = read_gene_input('rank', 'E2F1\nE2F2\nE2F3', method)
        draft = input_draft(gene_text, gene_file, method=method)
        submitted = input_actions('rank')
    if submitted:
        try:
            parsed = parse_input(gene_text, gene_file)
            matched, missing = match_genes(gene_rankings, parsed['genes'])
            if matched:
                st.session_state.rank_result = build_result(matched, {'targets': {'parsed': parsed, 'matched': matched, 'missing': missing}}, draft=draft)
            else:
                show_matches(parsed, matched, missing)
                st.warning(ui('No matched genes. The previous result is unchanged.', '没有匹配的基因，保留上次运行结果。'))
        except (ValueError, UnicodeError) as exc:
            st.error(str(exc))
    result = st.session_state.get('rank_result')
    render_run_status('rank', draft, result)
    if result:
        match = result['matches']['targets']
        show_matches(match['parsed'], match['matched'], match['missing'])
        fig = create_rank_plot(result['rankings'], result['genes'], essential_gene, nonessential_gene,
                               n_cell_lines, show_labels, point_size)
        centered_plot(fig)
        result_downloads(result, fig, 'rank')

elif view == 'box':
    st.caption(ui('Compare score distributions across cancer types for each gene.', '比较每个基因在不同癌种中的分数分布。'))
    if not lineage_col:
        st.info(t('lineage_missing'))
    else:
        with st.container(border=True):
            method = st.radio(t('input_method'), ['text', 'file'], horizontal=True, key='box_method',
                              format_func={'text': t('input_direct'), 'file': t('input_file')}.get,
                              label_visibility='collapsed')
            gene_text, gene_file = read_gene_input('box', 'E2F1\nE2F2', method)
            with st.expander(ui('Comparison options', '比较选项')):
                box_lineages = st.multiselect(ui('Cancer types shown in the chart', '图中展示的癌种'),
                                              cohort or all_lineages, key='box_lineages',
                                              placeholder=ui('All types in the analysis scope', '分析范围内的全部癌种'),
                                              help=ui('This only limits the chart. Gene ranks still use the analysis scope in the sidebar.',
                                                      '此处仅限制图中分组；基因排名仍按侧栏的分析范围计算。'))
                sort_mode = st.selectbox(ui('Cancer-type order', '癌种顺序'), ['median', 'alphabetical'], key='box_order',
                                         format_func={'median': ui('Median score', '中位数'), 'alphabetical': ui('Alphabetical', '字母顺序')}.get)
                all_points = st.checkbox(ui('Show all cell lines', '显示全部细胞系散点'), key='box_points')
            draft = input_draft(gene_text, gene_file, method=method, lineages=box_lineages, sort=sort_mode, all_points=all_points)
            submitted = input_actions('box')
        candidate = st.session_state.get('box_candidates')
        if candidate and candidate['input_signature'] != input_signature(draft):
            st.session_state.pop('box_candidates', None)
        if submitted:
            try:
                parsed = parse_input(gene_text, gene_file)
                matched, missing = match_genes(gene_rankings, parsed['genes'])
                match = {'targets': {'parsed': parsed, 'matched': matched, 'missing': missing}}
                config = {'lineages': box_lineages, 'sort': sort_mode, 'all_points': all_points}
                st.session_state.pop('box_candidates', None)
                if len(matched) > 8:
                    st.session_state.box_candidates = {'genes': matched, 'matches': match, 'config': config,
                                                       'draft': copy.deepcopy(draft), 'input_signature': input_signature(draft)}
                    st.session_state.pop('box_display_genes', None)
                elif matched:
                    frame = get_lineage_frame(df_scope, matched, gene_rankings, lineage_col)
                    if box_lineages:
                        frame = frame[frame.lineage.isin(box_lineages)]
                    st.session_state.box_result = build_result(matched, match, config, frame, draft)
                else:
                    show_matches(parsed, matched, missing)
                    st.warning(ui('No matched genes. The previous result is unchanged.', '没有匹配的基因，保留上次运行结果。'))
            except (ValueError, UnicodeError) as exc:
                st.error(str(exc))
        candidate = st.session_state.get('box_candidates')
        if candidate:
            st.info(ui('More than 8 genes matched. Choose up to 8 to plot; the full input stays in the analysis metadata.',
                       '匹配超过 8 个基因，请选择最多 8 个绘图；完整输入仍保存在分析参数中。'))
            with st.form('box_selection'):
                chosen = st.multiselect(ui('Genes to display', '选择展示基因'), candidate['genes'],
                                        default=candidate['genes'][:8], max_selections=8, key='box_display_genes')
                apply_selection = st.form_submit_button(ui('Run selected genes', '运行所选基因'), key='box_selection_run', type='primary')
            if apply_selection and chosen:
                frame = get_lineage_frame(df_scope, chosen, gene_rankings, lineage_col)
                if candidate['config']['lineages']:
                    frame = frame[frame.lineage.isin(candidate['config']['lineages'])]
                st.session_state.box_result = build_result(chosen, candidate['matches'], candidate['config'], frame, candidate['draft'])
                del st.session_state.box_candidates
                st.rerun()
        result = st.session_state.get('box_result')
        render_run_status('box', draft, result)
        if result:
            match = result['matches']['targets']
            show_matches(match['parsed'], match['matched'], match['missing'])
            st.caption(ui('Displayed genes: ', '实际展示基因：') + ', '.join(result['genes']))
            valid_data = result['plot_data'].dropna(subset=['crispr_score'])
            if valid_data.empty:
                st.warning(ui('No finite scores in the selected groups.', '所选分组没有有效分数。'))
            else:
                absent = [g for g in result['genes'] if g not in set(valid_data.gene)]
                if absent:
                    st.info(ui('No valid scores in these groups: ', '这些基因在所选分组没有有效分数：') + ', '.join(absent))
                fig = create_lineage_boxplot(valid_data, result['genes'], result['config']['sort'], result['config']['all_points'])
                st.plotly_chart(fig, config=PLOT_CONFIG, width='stretch')
                result_downloads(result, fig, 'box')

elif view == 'multi':
    st.caption(ui('Mark a background gene set and highlight selected genes within the same ranking.',
                  '在同一排名中标注背景基因集，并突出显示关注基因。'))
    with st.container(border=True):
        c1, c2 = st.columns(2)
        with c1:
            bg_text = st.text_area(t('bg_gene_set'), 'CDK1\nCDK2\nCCNB1\nCCND1\nCCNE1', key='multi_bg', height=110)
        with c2:
            hl_text = st.text_area(t('hl_gene_set'), 'PLK1\nAURKA', key='multi_hl', height=110)
        with st.expander(ui('Annotation colors', '标注颜色')):
            color1, color2 = st.columns(2)
            with color1:
                bg_color = st.color_picker(t('bg_color'), '#8DA8C0', key='multi_bg_color')
            with color2:
                hl_color = st.color_picker(t('hl_color'), '#0072B2', key='multi_hl_color')
        draft = {'background': bg_text, 'highlight': hl_text, 'background_color': bg_color, 'highlight_color': hl_color}
        submitted = input_actions('multi')
    if submitted:
        bg_parsed, hl_parsed = parse_gene_text(bg_text), parse_gene_text(hl_text)
        bg_matched, bg_missing = match_genes(gene_rankings, bg_parsed['genes'])
        hl_matched, hl_missing = match_genes(gene_rankings, hl_parsed['genes'])
        matches = {'background': {'parsed': bg_parsed, 'matched': bg_matched, 'missing': bg_missing},
                   'highlight': {'parsed': hl_parsed, 'matched': hl_matched, 'missing': hl_missing}}
        if bg_matched or hl_matched:
            genes = list(dict.fromkeys(bg_matched + hl_matched))
            st.session_state.multi_result = build_result(genes, matches, {'background': bg_matched, 'highlight': hl_matched,
                                                                          'background_color': bg_color, 'highlight_color': hl_color}, draft=draft)
        else:
            for key, title in [('background', t('bg_gene_set')), ('highlight', t('hl_gene_set'))]:
                m = matches[key]
                show_matches(m['parsed'], m['matched'], m['missing'], title)
            st.warning(ui('No matched genes. The previous result is unchanged.', '没有匹配的基因，保留上次运行结果。'))
    result = st.session_state.get('multi_result')
    render_run_status('multi', draft, result)
    if result:
        for key, title in [('background', t('bg_gene_set')), ('highlight', t('hl_gene_set'))]:
            match = result['matches'][key]
            show_matches(match['parsed'], match['matched'], match['missing'], title)
        cfg = result['config']
        fig = create_multilayer_rank_plot(result['rankings'], cfg['background'], cfg['highlight'],
                                          cfg['background_color'], cfg['highlight_color'],
                                          essential_gene, nonessential_gene, n_cell_lines, show_labels)
        centered_plot(fig)
        result_downloads(result, fig, 'multi')


render_diagnostics()

# ---- Gene × Drug correlation (temporarily hidden from the public UI) ----
if ENABLE_GENE_DRUG_UI:
    st.markdown(f"### {t('corr_title')}")
    st.markdown(t('corr_desc'))

    with st.spinner(t('corr_loading')):
        gdsc_df, crispr26_df, compounds_df, model_df, corr_err = load_corr_datasets(
            HF_REPO_ID, HF_GDSC_FILENAME, HF_CRISPR26Q1_FILENAME,
            HF_COMPOUNDS_FILENAME, HF_MODEL_FILENAME)

    if corr_err is not None:
        st.error(f"❌ {corr_err}")
        st.caption(
            "请确认 HF 数据集仓库内包含："
            f"`{HF_GDSC_FILENAME}` / `{HF_CRISPR26Q1_FILENAME}` / `{HF_COMPOUNDS_FILENAME}`"
        )
    else:
        c_left, c_right = st.columns(2)
        with c_left:
            corr_gene = st.text_input(t('corr_gene_label'),
                                      value=DEFAULT_CORR_GENE,
                                      help=t('corr_gene_help'),
                                      key="corr_gene")
        with c_right:
            drug_query = st.text_input(t('corr_drug_search'),
                                       value=DEFAULT_CORR_DRUG_QUERY,
                                       help=t('corr_drug_search_help'),
                                       key="corr_drug_query")

        gdsc2_only = st.checkbox(t('corr_gdsc2_only'), value=True,
                                 help=t('corr_gdsc2_only_help'),
                                 key="corr_gdsc2_only")

        gdsc_ids = set(gdsc_df.columns)
        matches = search_compounds(compounds_df, drug_query,
                                   gdsc2_only=gdsc2_only, gdsc_drug_ids=gdsc_ids)
        selected_cid = None
        if matches:
            labels = [lbl for _, lbl in matches]
            sel_label = st.selectbox(
                f"{t('corr_drug_select')}  ({len(labels)})",
                labels, key="corr_drug_sel")
            selected_cid = matches[labels.index(sel_label)][0]
        elif drug_query.strip():
            st.warning(t('corr_no_drug'))

        # ---- 肿瘤类型分层控件（仅当 Model.csv 加载成功时显示）----
        lineage_available = model_df is not None
        lin_gran_col = LINEAGE_COL_COARSE
        restrict_lineage = None
        if lineage_available:
            lc1, lc2 = st.columns(2)
            with lc1:
                gran_label = st.radio(
                    t('corr_lineage_gran'),
                    [t('corr_gran_coarse'), t('corr_gran_fine')],
                    horizontal=True, key="corr_gran")
                lin_gran_col = (LINEAGE_COL_COARSE
                                if gran_label == t('corr_gran_coarse')
                                else LINEAGE_COL_FINE)
            with lc2:
                if lin_gran_col in model_df.columns:
                    opts = [t('corr_all_lineages')] + sorted(
                        model_df[lin_gran_col].dropna().unique().tolist())
                    pick = st.selectbox(t('corr_restrict_lineage'), opts,
                                        key="corr_restrict")
                    restrict_lineage = None if pick == t('corr_all_lineages') else pick
        else:
            st.caption(f"ℹ️ {t('corr_no_model')}")

        run = st.button(f"▶️ {t('corr_run')}", key="corr_run_btn",
                        width="content")

        if run and selected_cid:
            gene_col = find_crispr_gene_column(crispr26_df, corr_gene)
            if gene_col is None:
                st.error(f"❌ {t('corr_gene_not_found')} ({corr_gene})")
            else:
                merged, rho, p, n, err = compute_gene_drug_correlation(
                    crispr26_df, gdsc_df, gene_col, selected_cid)
                drug_label = next((lbl for cid, lbl in matches
                                   if cid == selected_cid), selected_cid)
                drug_short = drug_label.split(' (')[0]

                # 附加 lineage
                lin_key = None
                if lineage_available and err != 'drug_not_found':
                    merged, lin_key = attach_lineage(merged, model_df, lin_gran_col)
                    if restrict_lineage is not None and lin_key:
                        merged = merged[merged['lineage'] == restrict_lineage]
                        if len(merged) >= 10:
                            rho, p = _spearman_with_p(
                                merged['dep'].values, merged['auc'].values)
                            n = len(merged)
                        else:
                            err = 'too_few'
                            n = len(merged)

                if err == 'drug_not_found':
                    st.error(f"❌ {t('corr_drug_not_found')}")
                elif err == 'too_few':
                    st.warning(f"⚠️ {t('corr_too_few')} (n = {n})")
                else:
                    m1, m2, m3 = st.columns(3)
                    m1.metric(t('corr_stat_n'), f"{n:,}")
                    m2.metric(t('corr_stat_rho'), f"{rho:.3f}")
                    m3.metric(t('corr_stat_p'), f"{p:.2g}")

                    title_suffix = f" · {restrict_lineage}" if restrict_lineage else ""
                    # A：散点（有 lineage 且未限定单一癌种时按癌种着色）
                    if lin_key and restrict_lineage is None:
                        fig = create_lineage_scatter(
                            merged, corr_gene.upper(), drug_short + title_suffix,
                            rho, p, point_size=point_size)
                    else:
                        fig = create_correlation_scatter(
                            merged, corr_gene.upper(), drug_short + title_suffix,
                            rho, p, point_size=point_size)
                    centered_plot(fig)

                    if p is not None and p < 0.05 and rho > 0:
                        st.success(t('corr_result_dir_pos'))
                    elif p is not None and p < 0.05 and rho < 0:
                        st.info(t('corr_result_dir_neg'))
                    else:
                        st.info(t('corr_result_dir_ns'))

                    # B：分层相关（仅在有 lineage 且看全部癌种时）
                    if lin_key and restrict_lineage is None:
                        strat = stratified_correlation(merged, min_n=MIN_N_PER_GROUP)
                        if len(strat) >= 2:
                            st.markdown(f"#### {t('corr_strat_title')}")
                            st.caption(t('corr_strat_desc'))
                            forest = create_forest_plot(
                                strat, corr_gene.upper(), drug_short)
                            centered_plot(forest)
                            show_tbl = strat.copy()
                            show_tbl['rho'] = show_tbl['rho'].round(3)
                            show_tbl['p'] = show_tbl['p'].map(lambda v: f"{v:.2g}")
                            show_tbl['fdr'] = show_tbl['fdr'].map(lambda v: f"{v:.2g}")
                            show_tbl.columns = [t('corr_tbl_lineage'), 'n',
                                                'Spearman ρ', 'p', 'BH-FDR']
                            st.dataframe(show_tbl, width="stretch",
                                         hide_index=True)
                        else:
                            st.caption(f"ℹ️ {t('corr_strat_insufficient')}")

                    st.warning(t('corr_caveat'))

                    with st.expander(f"📥 {t('export_title')}", expanded=False):
                        render_download_buttons(fig, "gene_drug_correlation",
                                                 "corr", height=export_height)
                    with st.expander(f"📊 {t('download_csv')}", expanded=False):
                        cols_out = ['dep', 'auc'] + (['lineage'] if lin_key else [])
                        csv_out = merged[cols_out].reset_index().rename(
                            columns={'index': 'ModelID'}).to_csv(index=False)
                        st.download_button(
                            t('download_csv'), data=csv_out,
                            file_name=f"{corr_gene}_{selected_cid}_correlation.csv",
                            mime='text/csv', key="corr_csv_dl")


# =============================================================================
# 页脚：致谢 + 引用
# =============================================================================
st.markdown("---")
with st.expander(f"📚 {t('resources')}", expanded=False):
    st.markdown(f"""
    <div class="footer-card">
        <h4>🙏 {t('acknowledgements')}</h4>
        <div style="display: flex; flex-wrap: wrap; gap: 1.5rem;">
            <div style="flex: 1; min-width: 250px;">
                <p>
                    <strong>{t('data_from')}</strong><br>
                    <a href="https://depmap.org" target="_blank">DepMap Portal (Broad Institute)</a><br>
                    <span style="font-size: 0.8rem;">CRISPR Chronos dependency scores</span>
                </p>
            </div>
            <div style="flex: 1; min-width: 250px;">
                <p>
                    <strong>{t('dev_with')}</strong><br>
                    <a href="https://www.anthropic.com/claude" target="_blank">Claude (Anthropic)</a><br>
                    <span style="font-size: 0.8rem;">{t('ai_dev')}</span>
                </p>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)

    render_citation_section()

st.markdown(
    f'<div style="text-align:center; color:{get_theme()["text_muted"]}; '
    f'font-size:0.8rem; padding:1rem;">'
    f'CRISPR Score Analyzer {TOOL_VERSION} | Deng Lab | '
    f'<a href="{GITHUB_URL}" target="_blank" style="color:{get_theme()["accent"]};">GitHub</a>'
    f'</div>',
    unsafe_allow_html=True
)
