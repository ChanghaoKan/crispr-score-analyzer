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
        'app_subtitle': 'Gene essentiality analysis platform powered by DepMap',
        'sidebar_settings': 'Settings',
        'language': 'Language',
        'theme': 'Theme',
        'light': 'Light',
        'dark': 'Dark',
        'data_source': 'Data Source',
        'upload_custom': 'Upload custom data (optional)',
        'upload_csv': 'Upload CRISPR Score CSV',
        'reference_genes': 'Reference Genes',
        'essential': 'Essential',
        'nonessential': 'Non-essential',
        'display_settings': 'Display',
        'show_labels': 'Show gene labels',
        'point_size': 'Point size',
        'export_size': 'Export Size',
        'export_height': 'Image height (px)',
        'cell_lines': 'Cell lines',
        'gene_count': 'Genes',
        'essential_genes': 'Mean score < {threshold}',
        'score_range': 'Score range',
        'custom_dataset': 'Custom uploaded dataset',
        'score_guide': "Lower scores indicate stronger gene dependency. The count uses each gene's mean across all cell lines; {threshold} is a screening aid, not a universal biological cutoff.",
        'threshold_label': 'Screening cutoff {threshold}',
        'essential_status': 'Mean-score screen',
        'essential_yes': 'Mean below cutoff',
        'essential_no': 'Mean not below cutoff',
        'first_eight_only': 'Showing the first 8 genes to keep the figure readable.',
        'lineage_missing': 'No lineage column was found in this dataset.',
        'tab1': '📊 Ranking',
        'tab2': '📦 Lineage',
        'tab3': '🎯 Multi-layer',
        'tab4': '🔗 Gene × Drug',
        'gene_ranking_title': 'Gene Essentiality Ranking',
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
        'multilayer_title': 'Multi-layer Gene Annotation',
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
        'loading_hf': '🔄 Loading data from Hugging Face...',
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
        'app_title': 'CRISPR 基因必需性分析器',
        'app_subtitle': '基于 DepMap 数据的基因必需性分析平台',
        'sidebar_settings': '设置',
        'language': '语言',
        'theme': '主题',
        'light': '浅色',
        'dark': '深色',
        'data_source': '数据来源',
        'upload_custom': '上传自定义数据（可选）',
        'upload_csv': '上传 CRISPR Score CSV',
        'reference_genes': '参考基因',
        'essential': '必需基因',
        'nonessential': '非必需基因',
        'display_settings': '显示设置',
        'show_labels': '显示基因名标签',
        'point_size': '点大小',
        'export_size': '导出尺寸',
        'export_height': '图片高度 (px)',
        'cell_lines': '细胞系',
        'gene_count': '基因数',
        'essential_genes': '平均分 < {threshold}',
        'score_range': 'Score 范围',
        'custom_dataset': '自定义上传数据',
        'score_guide': '分数越低表示基因依赖越强；此处按各基因在全部细胞系中的平均分计数，{threshold} 仅是筛选参考，并非通用的生物学阈值。',
        'threshold_label': '筛选参考线 {threshold}',
        'essential_status': '平均分筛选',
        'essential_yes': '平均分低于阈值',
        'essential_no': '平均分未低于阈值',
        'first_eight_only': '为保证图表清晰，仅展示前 8 个基因。',
        'lineage_missing': '该数据中未找到 lineage 列。',
        'tab1': '📊 基因排名',
        'tab2': '📦 癌种箱线图',
        'tab3': '🎯 多层标注',
        'tab4': '🔗 基因×药物',
        'gene_ranking_title': '基因必需性排名',
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
        'multilayer_title': '多层基因标注',
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
        'loading_hf': '🔄 正在从 Hugging Face 加载数据...',
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
if 'lang' not in st.session_state:
    st.session_state.lang = 'en'
if 'theme' not in st.session_state:
    st.session_state.theme = 'light'


# =============================================================================
# 主题
# =============================================================================
THEMES = {
    'light': {
        'bg': '#ffffff',
        'bg_secondary': '#f5f7fa',
        'bg_card': '#ffffff',
        'text': '#1a1a1a',
        'text_muted': '#666',
        'border': '#e0e4e8',
        'accent': '#0062b3',
        'accent_hover': '#004a8a',
        'success': '#35b300',
        'danger': '#b30035',
        'shadow': '0 1px 3px rgba(0,0,0,0.08), 0 1px 2px rgba(0,0,0,0.04)',
        'shadow_hover': '0 4px 12px rgba(0,0,0,0.10)',
        'plot_bg': '#ffffff',
        'plot_text': '#1a1a1a',
        'plot_axis': '#1a1a1a',
        'plot_grid': 'rgba(0,0,0,0.05)',
        'plot_scatter_bg': 'rgba(140,140,140,0.4)',
    },
    'dark': {
        'bg': '#0e1117',
        'bg_secondary': '#1a1f2e',
        'bg_card': '#1a1f2e',
        'text': '#e8eaed',
        'text_muted': '#9aa0a6',
        'border': '#2d333d',
        'accent': '#4a9eff',
        'accent_hover': '#6fb4ff',
        'success': '#4ade80',
        'danger': '#f87171',
        'shadow': '0 1px 3px rgba(0,0,0,0.3), 0 1px 2px rgba(0,0,0,0.2)',
        'shadow_hover': '0 4px 12px rgba(0,0,0,0.4)',
        'plot_bg': '#1a1f2e',
        'plot_text': '#e8eaed',
        'plot_axis': '#e8eaed',
        'plot_grid': 'rgba(255,255,255,0.08)',
        'plot_scatter_bg': 'rgba(160,160,160,0.35)',
    }
}


def get_theme():
    return THEMES[st.session_state.get('theme', 'light')]


def inject_css():
    th = get_theme()
    st.markdown(f"""
    <style>
        .stApp {{ background-color: {th['bg']}; color: {th['text']}; }}
        section[data-testid="stSidebar"] {{
            background-color: {th['bg_secondary']};
            border-right: 1px solid {th['border']};
        }}
        section[data-testid="stSidebar"] * {{ color: {th['text']}; }}
        .stApp p, .stApp label, .stApp span, .stApp div {{ color: {th['text']}; }}

        .main-header {{
            font-family: 'Inter', 'Helvetica Neue', Arial, sans-serif;
            font-size: 2.5rem; font-weight: 700;
            background: linear-gradient(135deg, {th['accent']} 0%, {th['accent_hover']} 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
            margin-bottom: 0.3rem; letter-spacing: -1px; line-height: 1.1;
        }}
        .sub-header {{
            font-family: 'Inter', 'Helvetica Neue', Arial, sans-serif;
            font-size: 1.05rem; color: {th['text_muted']};
            margin-bottom: 1.8rem; font-weight: 400;
        }}
        h1, h2, h3, h4, h5 {{
            color: {th['text']} !important;
            font-family: 'Inter', 'Helvetica Neue', Arial, sans-serif;
        }}

        .gene-tag {{
            display: inline-block;
            background: linear-gradient(135deg, {th['accent']}22 0%, {th['accent']}11 100%);
            color: {th['accent']}; padding: 0.3rem 0.75rem;
            border-radius: 6px; margin: 0.2rem;
            font-size: 0.82rem; font-weight: 600;
            font-family: 'JetBrains Mono', 'Monaco', 'Consolas', monospace;
            border: 1px solid {th['accent']}33;
            transition: all 0.2s ease;
        }}
        .gene-tag:hover {{
            transform: translateY(-1px);
            box-shadow: {th['shadow_hover']};
        }}

        .custom-divider {{
            height: 1px;
            background: linear-gradient(90deg, transparent, {th['border']}, transparent);
            border: none; margin: 2rem 0;
        }}

        .input-section-title {{
            font-size: 0.8rem; font-weight: 700; color: {th['text_muted']};
            margin-bottom: 0.6rem; text-transform: uppercase; letter-spacing: 1px;
        }}

        div[data-testid="stMetric"] {{
            background: {th['bg_card']};
            padding: 1.2rem 1.4rem; border-radius: 12px;
            border: 1px solid {th['border']}; box-shadow: {th['shadow']};
            transition: all 0.3s ease;
        }}
        div[data-testid="stMetric"]:hover {{
            transform: translateY(-2px);
            box-shadow: {th['shadow_hover']};
            border-color: {th['accent']}66;
        }}
        div[data-testid="stMetricLabel"] {{
            color: {th['text_muted']} !important;
            font-size: 0.8rem !important; font-weight: 600 !important;
            text-transform: uppercase; letter-spacing: 0.5px;
        }}
        div[data-testid="stMetricValue"] {{
            color: {th['text']} !important; font-weight: 700 !important;
        }}

        .data-context {{
            display: flex; align-items: center; gap: 0.9rem;
            background: {th['bg_secondary']};
            border: 1px solid {th['border']};
            border-left: 4px solid {th['accent']};
            border-radius: 10px; padding: 0.8rem 1rem;
            margin: 1rem 0 0.25rem 0;
            color: {th['text_muted']}; font-size: 0.88rem;
            line-height: 1.45;
        }}
        .data-context .data-chip {{
            flex: 0 0 auto; color: {th['accent']} !important;
            background: {th['accent']}14; border: 1px solid {th['accent']}33;
            border-radius: 999px; padding: 0.25rem 0.7rem;
            font-size: 0.78rem; font-weight: 700;
        }}
        .data-context > span:not(.data-chip) {{
            color: {th['text_muted']} !important;
        }}

        .stTabs [data-baseweb="tab-list"] {{
            gap: 6px; background: {th['bg_secondary']};
            padding: 6px; border-radius: 12px;
            border: 1px solid {th['border']};
            overflow-x: auto; scrollbar-width: none;
        }}
        .stTabs [data-baseweb="tab-list"]::-webkit-scrollbar {{ display: none; }}
        .stTabs [data-baseweb="tab"] {{
            border-radius: 8px; padding: 10px 20px;
            background: transparent; color: {th['text_muted']};
            font-weight: 600; transition: all 0.2s ease;
            flex: 1 0 auto; justify-content: center;
        }}
        .stTabs [data-baseweb="tab"]:hover {{
            background: {th['bg_card']}; color: {th['text']};
        }}
        .stTabs [aria-selected="true"] {{
            background: {th['accent']} !important; color: white !important;
            box-shadow: 0 2px 8px {th['accent']}44;
        }}

        .stDownloadButton button, .stButton button {{
            background: {th['accent']}; color: white; border: none;
            border-radius: 8px; padding: 0.5rem 1rem;
            font-weight: 600; transition: all 0.2s ease;
            box-shadow: 0 2px 6px {th['accent']}33;
        }}
        .stDownloadButton button:hover, .stButton button:hover {{
            background: {th['accent_hover']};
            transform: translateY(-1px);
            box-shadow: 0 4px 12px {th['accent']}55;
        }}

        .stTextInput input, .stTextArea textarea {{
            background: {th['bg_card']} !important;
            color: {th['text']} !important;
            border-color: {th['border']} !important;
            border-radius: 8px !important;
        }}
        .stTextInput input:focus, .stTextArea textarea:focus {{
            border-color: {th['accent']} !important;
            box-shadow: 0 0 0 2px {th['accent']}22 !important;
        }}

        .streamlit-expanderHeader {{
            background: {th['bg_card']} !important;
            border-radius: 8px !important;
            border: 1px solid {th['border']} !important;
            font-weight: 600 !important;
        }}

        .stAlert {{
            border-radius: 10px !important;
            border: 1px solid {th['border']} !important;
        }}

        .stDataFrame {{
            border-radius: 10px; overflow: hidden;
            border: 1px solid {th['border']};
        }}

        #MainMenu {{visibility: hidden;}}
        footer {{visibility: hidden;}}
        header[data-testid="stHeader"] {{ background: transparent; }}

        .stRadio label {{ color: {th['text']} !important; }}

        .footer-card {{
            background: {th['bg_secondary']};
            border: 1px solid {th['border']};
            padding: 1.8rem; border-radius: 16px;
            margin-bottom: 1rem;
        }}
        .footer-card h4 {{ color: {th['text']} !important; margin: 0 0 1rem 0; }}
        .footer-card p {{ color: {th['text_muted']} !important; font-size: 0.9rem; margin: 0; }}
        .footer-card a {{
            color: {th['accent']} !important;
            text-decoration: none; font-weight: 600;
        }}
        .footer-card a:hover {{ text-decoration: underline; }}

        .cite-box {{
            background: {th['bg_card']};
            border: 1px solid {th['border']};
            border-left: 4px solid {th['accent']};
            padding: 1rem 1.2rem; border-radius: 8px;
            font-family: 'JetBrains Mono', 'Monaco', monospace;
            font-size: 0.82rem; color: {th['text']};
            margin: 0.5rem 0; line-height: 1.6;
        }}
        .doi-badge {{
            display: inline-block;
            background: {th['accent']};
            color: white; padding: 0.2rem 0.7rem;
            border-radius: 4px; font-size: 0.75rem;
            font-weight: 600; margin-left: 0.5rem;
            text-decoration: none;
        }}
        .doi-badge:hover {{
            background: {th['accent_hover']}; color: white !important;
        }}

        button:focus-visible, input:focus-visible, textarea:focus-visible,
        [role="tab"]:focus-visible {{
            outline: 3px solid {th['accent']}66 !important;
            outline-offset: 2px !important;
        }}

        @media (max-width: 768px) {{
            .block-container {{
                padding: 1.25rem 0.9rem 2rem 0.9rem !important;
            }}
            .main-header {{
                font-size: 1.9rem; letter-spacing: -0.5px;
            }}
            .sub-header {{
                font-size: 0.95rem; margin-bottom: 1.2rem;
            }}
            .data-context {{
                align-items: flex-start; flex-direction: column;
                gap: 0.55rem; padding: 0.8rem 0.9rem;
            }}
            div[data-testid="stHorizontalBlock"]:has(div[data-testid="stMetric"]) {{
                flex-wrap: wrap; gap: 0.7rem;
            }}
            div[data-testid="stHorizontalBlock"]:has(div[data-testid="stMetric"])
            > div[data-testid="stColumn"] {{
                flex: 1 1 calc(50% - 0.7rem); min-width: 140px;
            }}
            div[data-testid="stMetric"] {{
                padding: 0.9rem 1rem;
            }}
            div[data-testid="stMetricValue"] {{
                font-size: 1.75rem !important;
            }}
            .stTabs [data-baseweb="tab"] {{
                padding: 8px 12px; font-size: 0.84rem;
            }}
            .footer-card {{ padding: 1.2rem; }}
        }}
    </style>
    """, unsafe_allow_html=True)


# =============================================================================
# 数据加载
# =============================================================================
@st.cache_resource(show_spinner=False)
def download_from_huggingface(repo_id: str, filename: str):
    try:
        from huggingface_hub import hf_hub_download
    except ImportError:
        return None, False, "huggingface_hub not installed"
    try:
        path = hf_hub_download(repo_id=repo_id, filename=filename, repo_type="dataset")
        df = pd.read_csv(path, low_memory=False, memory_map=True)
        return df, True, None
    except Exception as e:
        return None, False, f"HF error: {str(e)}"


@st.cache_resource(show_spinner=False, max_entries=3)
def load_uploaded_data(file_hash: str, _file_content: bytes):
    """按内容摘要缓存上传数据；返回的数据框在应用中只读使用。"""
    if not file_hash:
        raise ValueError("Missing upload content hash")
    return pd.read_csv(io.BytesIO(_file_content), low_memory=False)


def extract_gene_name(col_name: str) -> str:
    match = re.match(r'^([A-Za-z0-9_.-]+)\s*\(', str(col_name))
    if match:
        return match.group(1)
    return str(col_name)


@st.cache_data(show_spinner=False, max_entries=4)
def compute_gene_rankings(df_hash: str, _df: pd.DataFrame):
    """向量化识别并汇总基因列；_df 不参与 Streamlit 的重复哈希。"""
    if not df_hash:
        return None, 0, "Missing dataset cache key"

    numeric_df = _df.select_dtypes(include=[np.number])
    if numeric_df.empty:
        return None, 0, "No numeric CRISPR score columns detected"

    stats = pd.DataFrame({
        'count': numeric_df.count(),
        'mean': numeric_df.mean(),
        'std': numeric_df.std(),
    })
    valid_mask = (
        (stats['count'] > 10)
        & stats['mean'].between(-5, 2, inclusive='neither')
        & (stats['std'] > 0.01)
    )
    gene_cols = stats.index[valid_mask].tolist()
    if not gene_cols:
        return None, 0, "No CRISPR score columns detected"
    mean_scores = stats.loc[gene_cols, 'mean'].sort_values()
    rankings = pd.DataFrame({
        'gene_raw': mean_scores.index,
        'gene': [extract_gene_name(col) for col in mean_scores.index],
        'mean_score': mean_scores.values,
        'rank': range(1, len(mean_scores) + 1),
        'percentile': [(i / len(mean_scores)) * 100 for i in range(1, len(mean_scores) + 1)]
    })
    rankings['gene_upper'] = rankings['gene'].str.upper()
    return rankings, len(_df), None


def filter_genes_by_list(gene_rank_df, gene_list):
    gene_list_upper = [g.upper() for g in gene_list]
    matched_mask = gene_rank_df['gene_upper'].isin(gene_list_upper)
    matched_genes = gene_rank_df[matched_mask]['gene'].tolist()
    matched_upper = set(gene_rank_df[matched_mask]['gene_upper'])
    not_found = [g for g in gene_list if g.upper() not in matched_upper]
    return matched_genes, not_found


def get_lineage_data(df, genes):
    lineage_col = None
    for col in df.columns:
        if 'lineage' in col.lower() and 'sub' not in col.lower():
            lineage_col = col
            break
    if lineage_col is None:
        return None
    raw_col_map = {}
    for col in df.columns:
        raw_col_map[extract_gene_name(col).upper()] = col
    result = []
    for gene in genes:
        actual_col = raw_col_map.get(gene.upper())
        if actual_col and actual_col in df.columns:
            temp = df[[lineage_col, actual_col]].copy()
            temp.columns = ['lineage', 'crispr_score']
            temp['gene'] = gene
            result.append(temp)
    if result:
        return pd.concat(result, ignore_index=True)
    return None


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
    'threshold': '#E69F00',
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
        font=dict(family=FONT_FAMILY, color=th['plot_text']),
    )
    fig.update_xaxes(
        linecolor=th['plot_axis'], tickcolor=th['plot_axis'],
        tickfont=dict(color=th['plot_text']),
        title_font=dict(color=th['plot_text']),
    )
    fig.update_yaxes(
        linecolor=th['plot_axis'], tickcolor=th['plot_axis'],
        tickfont=dict(color=th['plot_text']),
        title_font=dict(color=th['plot_text']),
    )
    return fig


# =============================================================================
# 绘图函数
# =============================================================================
def create_rank_plot(gene_rank_df, genes_of_interest, essential_gene='MYC',
                     nonessential_gene='PTEN', n_cell_lines=0,
                     show_labels=True, point_size=4):
    th = get_theme()
    fig = go.Figure()
    y_min = gene_rank_df['mean_score'].min()
    y_max = gene_rank_df['mean_score'].max()
    y_range = y_max - y_min

    highlight_set = set(g.upper() for g in genes_of_interest)
    highlight_set.add(essential_gene.upper())
    highlight_set.add(nonessential_gene.upper())
    bg_df = gene_rank_df[~gene_rank_df['gene_upper'].isin(highlight_set)]

    # ✅ CHANGED: go.Scatter → go.Scattergl（WebGL 渲染 18000 点，大幅提升流畅度）
    fig.add_trace(go.Scattergl(
        x=bg_df['rank'], y=bg_df['mean_score'], mode='markers',
        marker=dict(size=3, color=th['plot_scatter_bg']),
        name='All genes',
        hovertemplate='<b>%{text}</b><br>Rank: %{x:,}<br>Score: %{y:.4f}<extra></extra>',
        text=bg_df['gene']
    ))

    fig.add_hline(
        y=ESSENTIALITY_THRESHOLD,
        line=dict(dash="dash", color=PLOT_COLORS['threshold'], width=1.5),
        annotation_text=t('threshold_label').format(threshold=ESSENTIALITY_THRESHOLD),
        annotation_position="top left",
    )
    fig.add_hline(y=0, line=dict(color=th['plot_axis'], width=0.8))

    ess_df = gene_rank_df[gene_rank_df['gene_upper'] == essential_gene.upper()]
    if len(ess_df) > 0:
        fig.add_trace(go.Scatter(
            x=ess_df['rank'], y=ess_df['mean_score'], mode='markers+text',
            marker=dict(size=point_size * 2.2, color=PLOT_COLORS['essential'], symbol='diamond'),
            text=[essential_gene], textposition='bottom center',
            textfont=dict(size=11, color=PLOT_COLORS['essential'], family=FONT_FAMILY),
            name=f'Essential ({essential_gene})',
            hovertemplate=f'<b>{essential_gene}</b><br>Rank: %{{x:,}}<br>Score: %{{y:.4f}}<extra></extra>'
        ))

    noness_df = gene_rank_df[gene_rank_df['gene_upper'] == nonessential_gene.upper()]
    if len(noness_df) > 0:
        fig.add_trace(go.Scatter(
            x=noness_df['rank'], y=noness_df['mean_score'], mode='markers+text',
            marker=dict(size=point_size * 2.2, color=PLOT_COLORS['nonessential'], symbol='diamond'),
            text=[nonessential_gene], textposition='top center',
            textfont=dict(size=11, color=PLOT_COLORS['nonessential'], family=FONT_FAMILY),
            name=f'Non-essential ({nonessential_gene})',
            hovertemplate=f'<b>{nonessential_gene}</b><br>Rank: %{{x:,}}<br>Score: %{{y:.4f}}<extra></extra>'
        ))

    interest_df = gene_rank_df[gene_rank_df['gene'].isin(genes_of_interest)].copy()
    interest_df = interest_df.sort_values('rank').reset_index(drop=True)
    if len(interest_df) > 0:
        text_positions = ['top center' if i % 2 == 0 else 'bottom center'
                          for i in range(len(interest_df))]
        fig.add_trace(go.Scatter(
            x=interest_df['rank'], y=interest_df['mean_score'],
            mode='markers+text' if show_labels else 'markers',
            marker=dict(size=point_size * 2.5, color=PLOT_COLORS['interest'],
                        line=dict(width=1.5, color=th['plot_bg'])),
            text=interest_df['gene'] if show_labels else None,
            textposition=text_positions,
            textfont=dict(size=11, color=PLOT_COLORS['interest'], family=FONT_FAMILY),
            name='Genes of interest',
            hovertemplate='<b>%{text}</b><br>Rank: %{x:,}<br>Score: %{y:.4f}<br>Percentile: %{customdata:.1f}%<extra></extra>',
            customdata=interest_df['percentile']
        ))

    y_label = (f"Mean CRISPR Score<br><span style='font-size:11px'>"
               f"({n_cell_lines} cell lines)</span>" if n_cell_lines > 0
               else "Mean CRISPR Score")
    y_tickvals = np.arange(np.floor(y_min / 0.5) * 0.5,
                           np.ceil(y_max / 0.5) * 0.5 + 0.5, 0.5)

    fig.update_layout(
        title=dict(text='<b>Gene Dependency Ranking</b>',
                   font=dict(size=16, family=FONT_FAMILY),
                   x=0.5, xanchor='center', xref='paper'),
        xaxis=dict(title='Gene Rank', showgrid=False, showline=True, linewidth=1.5,
                   tickformat=',d', ticks='outside', ticklen=5,
                   range=[0, len(gene_rank_df) * 1.02]),
        yaxis=dict(title=y_label, showgrid=False, showline=True, linewidth=1.5,
                   tickvals=y_tickvals, ticks='outside', ticklen=5,
                   range=[y_min - 0.1 * y_range, y_max + 0.15 * y_range]),
        legend=dict(orientation='v', yanchor='bottom', y=0.02,
                    xanchor='right', x=0.98,
                    font=dict(size=11),
                    bgcolor=th['plot_bg'], borderwidth=1,
                    bordercolor=th['border']),
        height=650, margin=dict(l=80, r=80, t=60, b=60)
    )
    return apply_theme_to_fig(fig)


def create_lineage_boxplot(lineage_data, genes):
    th = get_theme()
    n_genes = len(genes)
    v_spacing = min(0.15, 0.6 / n_genes)
    fig = make_subplots(rows=n_genes, cols=1, shared_xaxes=True,
                        vertical_spacing=v_spacing,
                        subplot_titles=[f'<i>{g}</i>' for g in genes])
    lineages = sorted(lineage_data['lineage'].unique())

    for i, gene in enumerate(genes, 1):
        gene_data = lineage_data[lineage_data['gene'] == gene]
        fig.add_trace(go.Box(
            x=gene_data['lineage'], y=gene_data['crispr_score'], name=gene,
            marker=dict(color=PLOT_COLORS['boxplot_fill']),
            line=dict(color=th['plot_axis'], width=1),
            fillcolor=PLOT_COLORS['boxplot_fill'], showlegend=False, boxpoints=False
        ), row=i, col=1)
        fig.add_hline(y=0, line=dict(dash="dot", color=th['plot_axis'], width=1),
                      row=i, col=1)
        threshold_kwargs = {}
        if i == 1:
            threshold_kwargs = {
                'annotation_text': t('threshold_label').format(
                    threshold=ESSENTIALITY_THRESHOLD),
                'annotation_position': 'top left',
            }
        fig.add_hline(
            y=ESSENTIALITY_THRESHOLD,
            line=dict(dash="dash", color=PLOT_COLORS['threshold'], width=1.2),
            row=i, col=1, **threshold_kwargs,
        )

    fig.update_layout(
        title=dict(text='<b>CRISPR Score by Cancer Type</b>',
                   font=dict(size=16, family=FONT_FAMILY),
                   x=0.5, xanchor='center', xref='paper'),
        height=280 * n_genes + 100, showlegend=False,
        margin=dict(l=80, r=80, t=80, b=80)
    )
    fig.update_xaxes(tickangle=-45, categoryarray=lineages,
                     showline=True, linewidth=1.5, ticks='outside', ticklen=5)
    fig.update_yaxes(title_text='CRISPR Score', showgrid=False,
                     showline=True, linewidth=1.5,
                     ticks='outside', ticklen=5, dtick=0.5)
    for i in range(1, n_genes):
        fig.update_xaxes(showticklabels=False, row=i, col=1)
    return apply_theme_to_fig(fig)


def create_multilayer_rank_plot(gene_rank_df, background_genes, highlight_genes,
                                 bg_color='#7FB3D5', hl_color='#E74C3C',
                                 essential_gene='MYC', nonessential_gene='PTEN',
                                 n_cell_lines=0, show_labels=True):
    th = get_theme()
    fig = go.Figure()
    y_min = gene_rank_df['mean_score'].min()
    y_max = gene_rank_df['mean_score'].max()
    y_range = y_max - y_min

    all_highlight = set(g.upper() for g in background_genes + highlight_genes)
    all_highlight.add(essential_gene.upper())
    all_highlight.add(nonessential_gene.upper())
    bg_all_df = gene_rank_df[~gene_rank_df['gene_upper'].isin(all_highlight)]

    # ✅ CHANGED: go.Scatter → go.Scattergl（WebGL 渲染背景散点）
    fig.add_trace(go.Scattergl(
        x=bg_all_df['rank'], y=bg_all_df['mean_score'], mode='markers',
        marker=dict(size=2.5, color=th['plot_scatter_bg']),
        name='All genes',
        hovertemplate='<b>%{text}</b><br>Rank: %{x:,}<br>Score: %{y:.4f}<extra></extra>',
        text=bg_all_df['gene']
    ))
    fig.add_hline(
        y=ESSENTIALITY_THRESHOLD,
        line=dict(dash="dash", color=PLOT_COLORS['threshold'], width=1.5),
        annotation_text=t('threshold_label').format(threshold=ESSENTIALITY_THRESHOLD),
        annotation_position="top left",
    )
    fig.add_hline(y=0, line=dict(color=th['plot_axis'], width=0.8))

    ess_df = gene_rank_df[gene_rank_df['gene_upper'] == essential_gene.upper()]
    if len(ess_df) > 0:
        fig.add_trace(go.Scatter(
            x=ess_df['rank'], y=ess_df['mean_score'], mode='markers+text',
            marker=dict(size=9, color=PLOT_COLORS['essential'], symbol='diamond'),
            text=[essential_gene], textposition='bottom center',
            textfont=dict(size=11, color=PLOT_COLORS['essential']),
            name=f'Essential ({essential_gene})'
        ))
    noness_df = gene_rank_df[gene_rank_df['gene_upper'] == nonessential_gene.upper()]
    if len(noness_df) > 0:
        fig.add_trace(go.Scatter(
            x=noness_df['rank'], y=noness_df['mean_score'], mode='markers+text',
            marker=dict(size=9, color=PLOT_COLORS['nonessential'], symbol='diamond'),
            text=[nonessential_gene], textposition='top center',
            textfont=dict(size=11, color=PLOT_COLORS['nonessential']),
            name=f'Non-essential ({nonessential_gene})'
        ))

    bg_only = [g for g in background_genes if g not in highlight_genes]
    bg_df = gene_rank_df[gene_rank_df['gene'].isin(bg_only)]
    if len(bg_df) > 0:
        fig.add_trace(go.Scatter(
            x=bg_df['rank'], y=bg_df['mean_score'], mode='markers',
            marker=dict(size=7, color=bg_color, opacity=0.65),
            name=f'Gene set (n={len(bg_df)})',
            text=bg_df['gene'],
            hovertemplate='<b>%{text}</b><br>Rank: %{x:,}<br>Score: %{y:.4f}<extra></extra>'
        ))

    hl_df = gene_rank_df[gene_rank_df['gene'].isin(highlight_genes)].copy()
    hl_df = hl_df.sort_values('rank').reset_index(drop=True)
    if len(hl_df) > 0:
        text_positions = ['top center' if i % 2 == 0 else 'bottom center'
                          for i in range(len(hl_df))]
        fig.add_trace(go.Scatter(
            x=hl_df['rank'], y=hl_df['mean_score'],
            mode='markers+text' if show_labels else 'markers',
            marker=dict(size=11, color=hl_color, line=dict(width=1.5, color=th['plot_bg'])),
            text=hl_df['gene'] if show_labels else None,
            textposition=text_positions,
            textfont=dict(size=11, color=th['plot_text'], family=FONT_FAMILY),
            name=f'Highlight (n={len(hl_df)})',
            hovertemplate='<b>%{text}</b><br>Rank: %{x:,}<br>Score: %{y:.4f}<br>Percentile: %{customdata:.1f}%<extra></extra>',
            customdata=hl_df['percentile']
        ))

    y_label = (f"Mean CRISPR Score<br><span style='font-size:11px'>"
               f"({n_cell_lines} cell lines)</span>" if n_cell_lines > 0
               else "Mean CRISPR Score")
    y_tickvals = np.arange(np.floor(y_min / 0.5) * 0.5,
                           np.ceil(y_max / 0.5) * 0.5 + 0.5, 0.5)

    fig.update_layout(
        title=dict(text='<b>Multi-layer Gene Annotation</b>',
                   font=dict(size=16, family=FONT_FAMILY),
                   x=0.5, xanchor='center', xref='paper'),
        xaxis=dict(title='Gene Rank', showgrid=False, showline=True, linewidth=1.5,
                   tickformat=',d', ticks='outside', ticklen=5),
        yaxis=dict(title=y_label, showgrid=False, showline=True, linewidth=1.5,
                   tickvals=y_tickvals, ticks='outside', ticklen=5,
                   range=[y_min - 0.1 * y_range, y_max + 0.15 * y_range]),
        legend=dict(yanchor='bottom', y=0.02, xanchor='right', x=0.98,
                    font=dict(size=11), bgcolor=th['plot_bg'],
                    bordercolor=th['border']),
        height=650, margin=dict(l=80, r=80, t=60, b=60)
    )
    return apply_theme_to_fig(fig)


# =============================================================================
# 图片导出（始终白底，方便论文用）
# =============================================================================
def fig_for_export(fig):
    export_fig = copy.deepcopy(fig)
    export_fig.update_layout(
        plot_bgcolor='white', paper_bgcolor='white',
        font=dict(family=FONT_FAMILY, color='#1a1a1a'),
    )
    export_fig.update_xaxes(
        linecolor='black', tickcolor='black',
        tickfont=dict(color='black'), title_font=dict(color='black'),
    )
    export_fig.update_yaxes(
        linecolor='black', tickcolor='black',
        tickfont=dict(color='black'), title_font=dict(color='black'),
    )
    for trace in export_fig.data:
        if hasattr(trace, 'name') and trace.name == 'All genes':
            trace.marker.color = 'rgba(180,180,180,0.4)'
    return export_fig


# ✅ CHANGED: 导出改为按需生成，缓存到 session_state，不在页面加载时渲染
def render_download_buttons(fig, filename_base: str, key_prefix: str, height: int = 600):
    """按需生成导出文件，点击按钮后才调用 kaleido，生成后缓存避免重复渲染"""
    width = 1000

    st.caption(t('download_hint'))

    formats = [
        ('pdf', t('download_pdf'), 'application/pdf', {}),
        ('png', t('download_png'), 'image/png', {'scale': 3}),
        ('svg', t('download_svg'), 'image/svg+xml', {}),
    ]

    cols = st.columns(len(formats))

    for col, (fmt, label, mime, extra_kwargs) in zip(cols, formats):
        cache_key = f"_export_{key_prefix}_{fmt}"

        with col:
            # 如果缓存里已有，直接显示下载按钮
            if cache_key in st.session_state and st.session_state[cache_key] is not None:
                st.download_button(
                    label=f"⬇️ {label}",
                    data=st.session_state[cache_key],
                    file_name=f"{filename_base}.{fmt}",
                    mime=mime,
                    key=f"{key_prefix}_{fmt}_dl",
                    width="stretch"
                )
            else:
                # 首次：点击按钮才生成
                if st.button(label, key=f"{key_prefix}_{fmt}_btn",
                             width="stretch"):
                    with st.spinner(f"Generating {fmt.upper()}..."):
                        try:
                            export_fig = fig_for_export(fig)
                            img_bytes = export_fig.to_image(
                                format=fmt, width=width, height=height,
                                engine='kaleido', **extra_kwargs
                            )
                            st.session_state[cache_key] = img_bytes
                            st.rerun()
                        except Exception as e:
                            st.warning(f"{fmt.upper()}: {str(e)[:60]}")


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


# =============================================================================
# CSS 注入
# =============================================================================
inject_css()


# =============================================================================
# 侧边栏
# =============================================================================
with st.sidebar:
    st.markdown(f"## ⚙️ {t('sidebar_settings')}")

    lang_options = {'English': 'en', '中文': 'zh'}
    current_lang_label = 'English' if st.session_state.lang == 'en' else '中文'
    selected_lang = st.selectbox(
        f"🌐 {t('language')}",
        options=list(lang_options.keys()),
        index=list(lang_options.keys()).index(current_lang_label)
    )
    if lang_options[selected_lang] != st.session_state.lang:
        st.session_state.lang = lang_options[selected_lang]
        # 清除导出缓存（语言切换后 label 变了）
        for k in list(st.session_state.keys()):
            if k.startswith('_export_'):
                del st.session_state[k]
        st.rerun()

    theme_options = {f"☀️ {t('light')}": 'light', f"🌙 {t('dark')}": 'dark'}
    current_theme_label = (f"☀️ {t('light')}" if st.session_state.theme == 'light'
                           else f"🌙 {t('dark')}")
    selected_theme = st.selectbox(
        f"🎨 {t('theme')}",
        options=list(theme_options.keys()),
        index=list(theme_options.keys()).index(current_theme_label)
    )
    if theme_options[selected_theme] != st.session_state.theme:
        st.session_state.theme = theme_options[selected_theme]
        # 清除导出缓存（主题切换后导出图需要重新生成）
        for k in list(st.session_state.keys()):
            if k.startswith('_export_'):
                del st.session_state[k]
        st.rerun()

    st.markdown("---")
    st.markdown(f"### 📁 {t('data_source')}")
    if USE_HUGGINGFACE:
        st.info(f"🤗 HuggingFace\n`{HF_REPO_ID}`")

    with st.expander(f"📤 {t('upload_custom')}"):
        uploaded_file = st.file_uploader(t('upload_csv'), type=['csv'])

    st.markdown("---")
    st.markdown(f"### 🧬 {t('reference_genes')}")
    col1, col2 = st.columns(2)
    with col1:
        essential_gene = st.text_input(t('essential'), value="MYC")
    with col2:
        nonessential_gene = st.text_input(t('nonessential'), value="PTEN")

    st.markdown("---")
    st.markdown(f"### 🎨 {t('display_settings')}")
    show_labels = st.checkbox(t('show_labels'), value=True)
    point_size = st.slider(t('point_size'), 2, 8, 4)

    st.markdown("---")
    st.markdown(f"### 📐 {t('export_size')}")
    export_height = st.slider(t('export_height'), 400, 1000, 600, step=50)


# =============================================================================
# 数据加载
# =============================================================================
data_loaded = False
crispr_data = None
data_cache_key = None

if uploaded_file is not None:
    with st.spinner(t('loading_upload')):
        uploaded_bytes = uploaded_file.getvalue()
        upload_digest = hashlib.sha256(uploaded_bytes).hexdigest()
        crispr_data = load_uploaded_data(upload_digest, uploaded_bytes)
        data_cache_key = f"upload:{upload_digest}"
        st.success(f"{t('loaded')}: {uploaded_file.name}")
        data_loaded = True
elif USE_HUGGINGFACE:
    with st.spinner(t('loading_hf')):
        df_result, success, err = download_from_huggingface(HF_REPO_ID, HF_FILENAME)
        if success:
            crispr_data = df_result
            data_cache_key = f"hf:{HF_REPO_ID}:{HF_FILENAME}"
            data_loaded = True
        else:
            st.error(f"❌ {err}")


# =============================================================================
# 主界面
# =============================================================================
st.markdown(f'<h1 class="main-header">🧬 {t("app_title")}</h1>', unsafe_allow_html=True)
st.markdown(f'<p class="sub-header">{t("app_subtitle")}</p>', unsafe_allow_html=True)

if not data_loaded:
    st.warning(t('no_data_warn'))
    st.stop()

df = crispr_data
df_hash = f"{data_cache_key}:{df.shape}"
gene_rankings, n_cell_lines, error_msg = compute_gene_rankings(df_hash, df)

if gene_rankings is None:
    st.error(f"❌ {error_msg}")
    with st.expander("🔍 Data diagnostics"):
        st.write("**First 10 columns:**", list(df.columns[:10]))
        st.write("**Shape:**", df.shape)
    st.stop()

# 概览指标
col1, col2, col3, col4 = st.columns(4)
with col1:
    st.metric(t('cell_lines'), f"{n_cell_lines:,}")
with col2:
    st.metric(t('gene_count'), f"{len(gene_rankings):,}")
with col3:
    essential_count = (gene_rankings['mean_score'] < ESSENTIALITY_THRESHOLD).sum()
    st.metric(
        t('essential_genes').format(threshold=ESSENTIALITY_THRESHOLD),
        f"{essential_count:,}",
    )
with col4:
    st.metric(t('score_range'),
              f"{gene_rankings['mean_score'].min():.2f}–{gene_rankings['mean_score'].max():.2f}")

dataset_label = (t('custom_dataset') if uploaded_file is not None
                 else f"{DATA_VERSION} · {SCORE_TYPE}")
score_guide = t('score_guide').format(threshold=ESSENTIALITY_THRESHOLD)
st.markdown(
    f'<div class="data-context"><span class="data-chip">{dataset_label}</span>'
    f'<span>{score_guide}</span></div>',
    unsafe_allow_html=True,
)

st.markdown('<div class="custom-divider"></div>', unsafe_allow_html=True)


# =============================================================================
# Tabs
# =============================================================================
tab1, tab2, tab3 = st.tabs([t('tab1'), t('tab2'), t('tab3')])

# ---- Tab 1 ----
with tab1:
    st.markdown(f"### {t('gene_ranking_title')}")
    st.markdown(t('gene_ranking_desc'))
    st.markdown(f'<p class="input-section-title">📝 {t("input_target_genes")}</p>',
                unsafe_allow_html=True)

    input_method = st.radio(t('input_method'), [t('input_direct'), t('input_file')],
                            horizontal=True, label_visibility="collapsed", key="tab1_radio")

    genes_of_interest = []
    if input_method == t('input_direct'):
        gene_input = st.text_area(
            t('gene_list'),
            value="E2F1\nE2F2\nE2F3\nE2F4\nE2F5\nE2F6\nE2F7\nE2F8",
            height=150, help=t('gene_list_help'))
        genes_of_interest = [g.strip() for g in gene_input.replace(',', '\n').replace(' ', '\n').split('\n') if g.strip()]
    else:
        uploaded_genelist = st.file_uploader(t('input_file'), type=['csv', 'txt'], key="genelist1")
        if uploaded_genelist:
            content = uploaded_genelist.getvalue().decode('utf-8')
            if uploaded_genelist.name.endswith('.csv'):
                genes_of_interest = pd.read_csv(io.StringIO(content)).iloc[:, 0].dropna().astype(str).tolist()
            else:
                genes_of_interest = [g.strip() for g in content.split('\n') if g.strip()]

    if genes_of_interest:
        matched_genes, not_found = filter_genes_by_list(gene_rankings, genes_of_interest)
        col_a, col_b = st.columns([3, 1])
        with col_a:
            if matched_genes:
                st.markdown(f"**✓ {t('matched')}:**")
                st.markdown(' '.join([f'<span class="gene-tag">{g}</span>' for g in matched_genes]),
                            unsafe_allow_html=True)
        with col_b:
            if not_found:
                with st.expander(f"⚠️ {t('not_found')} ({len(not_found)})"):
                    st.write(", ".join(not_found))

        if matched_genes:
            fig = create_rank_plot(gene_rankings, matched_genes,
                                    essential_gene, nonessential_gene,
                                    n_cell_lines, show_labels, point_size)
            centered_plot(fig)

            # ✅ CHANGED: expanded=False，不在加载时展开导出区
            with st.expander(f"📥 {t('export_title')}", expanded=False):
                render_download_buttons(fig, "gene_ranking", "rank_plot", height=export_height)

            with st.expander(f"📋 {t('gene_details')}", expanded=False):
                detail = gene_rankings[gene_rankings['gene'].isin(matched_genes)].sort_values('mean_score').copy()
                status_col = t('essential_status')
                detail[status_col] = detail['mean_score'].apply(
                    lambda x: (f"◆ {t('essential_yes')}" if x < ESSENTIALITY_THRESHOLD
                               else f"◇ {t('essential_no')}")
                )
                st.dataframe(detail[['gene', 'rank', 'percentile', 'mean_score', status_col]].round(4),
                             width="stretch", hide_index=True)
                csv_data = detail[['gene', 'rank', 'percentile', 'mean_score']].to_csv(index=False)
                st.download_button(t('download_csv'), data=csv_data,
                                   file_name="gene_ranking_data.csv", mime="text/csv", key="rank_csv")

# ---- Tab 2 ----
with tab2:
    st.markdown(f"### {t('boxplot_title')}")
    st.markdown(f'<p class="input-section-title">📝 {t("input_target_genes")}</p>',
                unsafe_allow_html=True)

    input_method2 = st.radio(t('input_method'), [t('input_direct'), t('input_file')],
                              horizontal=True, label_visibility="collapsed", key="tab2_radio")
    genes_for_box = []
    if input_method2 == t('input_direct'):
        gene_input2 = st.text_area(t('gene_list'), value="E2F1\nE2F2",
                                    height=120, key="box_text")
        genes_for_box = [g.strip() for g in gene_input2.replace(',', '\n').split('\n') if g.strip()]
    else:
        uploaded2 = st.file_uploader(t('input_file'), type=['csv', 'txt'], key="box_file")
        if uploaded2:
            content = uploaded2.getvalue().decode('utf-8')
            if uploaded2.name.endswith('.csv'):
                genes_for_box = pd.read_csv(io.StringIO(content)).iloc[:, 0].dropna().astype(str).tolist()
            else:
                genes_for_box = [g.strip() for g in content.split('\n') if g.strip()]

    if genes_for_box:
        matched, not_found = filter_genes_by_list(gene_rankings, genes_for_box)
        if not_found:
            st.warning(f"{t('not_found')}: {', '.join(not_found)}")
        if matched:
            if len(matched) > 8:
                st.info(t('first_eight_only'))
                matched = matched[:8]
            lineage_data = get_lineage_data(df, matched)
            if lineage_data is not None:
                fig = create_lineage_boxplot(lineage_data, matched)
                st.plotly_chart(fig, config=PLOT_CONFIG)
                # ✅ CHANGED: expanded=False
                with st.expander(f"📥 {t('export_title')}", expanded=False):
                    box_height = max(280 * len(matched) + 100, 400)
                    render_download_buttons(fig, "lineage_boxplot", "boxplot", height=box_height)
            else:
                st.error(t('lineage_missing'))

# ---- Tab 3 ----
with tab3:
    st.markdown(f"### {t('multilayer_title')}")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"**{t('bg_gene_set')}**")
        bg_input = st.text_area("BG", value="CDK1\nCDK2\nCCNB1\nCCND1\nCCNE1",
                                 height=150, key="bg", label_visibility="collapsed")
        bg_color = st.color_picker(t('bg_color'), "#7FB3D5", key="bg_color")
    with col2:
        st.markdown(f"**{t('hl_gene_set')}**")
        hl_input = st.text_area("HL", value="PLK1\nAURKA", height=150, key="hl",
                                 label_visibility="collapsed")
        hl_color = st.color_picker(t('hl_color'), "#E74C3C", key="hl_color")

    bg_genes = [g.strip() for g in bg_input.replace(',', '\n').split('\n') if g.strip()]
    hl_genes = [g.strip() for g in hl_input.replace(',', '\n').split('\n') if g.strip()]

    if bg_genes or hl_genes:
        bg_matched, _ = filter_genes_by_list(gene_rankings, bg_genes)
        hl_matched, _ = filter_genes_by_list(gene_rankings, hl_genes)
        st.markdown(f"{t('bg_gene_set')}: {len(bg_matched)} | {t('hl_gene_set')}: {len(hl_matched)}")

        if bg_matched or hl_matched:
            fig = create_multilayer_rank_plot(gene_rankings, bg_matched, hl_matched,
                                                bg_color, hl_color,
                                                essential_gene, nonessential_gene,
                                                n_cell_lines, show_labels)
            centered_plot(fig)
            # ✅ CHANGED: expanded=False
            with st.expander(f"📥 {t('export_title')}", expanded=False):
                render_download_buttons(fig, "multilayer_annotation", "multilayer",
                                         height=export_height)


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
