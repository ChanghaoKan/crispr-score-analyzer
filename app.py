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

# =============================================================================
# 页面配置
# =============================================================================
st.set_page_config(
    page_title="CRISPR Score Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# =============================================================================
# 数据源配置 (Hugging Face)
# =============================================================================
HF_REPO_ID = "ChanghaoKan/crispr-depmap"
HF_FILENAME = "CRISPR_(DepMap_Public_25Q3+Score,_Chronos)_subsetted.csv"
USE_HUGGINGFACE = True

# =============================================================================
# Citation / DOI 配置
# =============================================================================
ZENODO_DOI = "10.5281/zenodo.19607603"  
TOOL_VERSION = "v1.0"
TOOL_AUTHORS = "Kan, C"
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
        'essential_genes': 'Essential genes',
        'score_range': 'Score range',
        'tab1': '📊 Gene Ranking',
        'tab2': '📦 Lineage Boxplot',
        'tab3': '🎯 Multi-layer',
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
        'download_pdf': '📄 Download PDF (vector)',
        'download_png': '🖼️ Download PNG (300 DPI)',
        'download_svg': '✏️ Download SVG (editable)',
        'download_csv': '📊 Download data (CSV)',
        'download_hint': '💡 Recommend **PDF** for publications (infinite scalable vector). Use **SVG** if you need to edit in Illustrator/Inkscape.',
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
        'essential_genes': '必需基因数',
        'score_range': 'Score 范围',
        'tab1': '📊 基因排名图',
        'tab2': '📦 Lineage 箱线图',
        'tab3': '🎯 多层标注',
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
        'download_pdf': '📄 下载 PDF (矢量)',
        'download_png': '🖼️ 下载 PNG (300 DPI)',
        'download_svg': '✏️ 下载 SVG (可编辑)',
        'download_csv': '📊 下载数据表 (CSV)',
        'download_hint': '💡 论文投稿推荐使用 **PDF** 格式（矢量图，无限缩放不失真）。需要二次编辑可选 **SVG**（可在 Illustrator/Inkscape 中修改）。',
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

        .stTabs [data-baseweb="tab-list"] {{
            gap: 6px; background: {th['bg_secondary']};
            padding: 6px; border-radius: 12px;
            border: 1px solid {th['border']};
        }}
        .stTabs [data-baseweb="tab"] {{
            border-radius: 8px; padding: 10px 20px;
            background: transparent; color: {th['text_muted']};
            font-weight: 600; transition: all 0.2s ease;
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
        df = pd.read_csv(path)
        return df, True, None
    except Exception as e:
        return None, False, f"HF error: {str(e)}"


def load_uploaded_data(file_content):
    return pd.read_csv(io.StringIO(file_content.decode('utf-8')))


def extract_gene_name(col_name: str) -> str:
    match = re.match(r'^([A-Za-z0-9_.-]+)\s*\(', str(col_name))
    if match:
        return match.group(1)
    return str(col_name)


@st.cache_data(show_spinner=False)
def compute_gene_rankings(df_hash: str, df: pd.DataFrame):
    gene_cols = []
    for col in df.columns:
        try:
            if pd.api.types.is_numeric_dtype(df[col]):
                sample_vals = df[col].dropna()
                if len(sample_vals) > 10:
                    mean_val = sample_vals.mean()
                    std_val = sample_vals.std()
                    if -5 < mean_val < 2 and std_val > 0.01:
                        gene_cols.append(col)
        except Exception:
            pass
    if not gene_cols:
        return None, 0, "No CRISPR score columns detected"
    mean_scores = df[gene_cols].mean().sort_values()
    rankings = pd.DataFrame({
        'gene_raw': mean_scores.index,
        'gene': [extract_gene_name(col) for col in mean_scores.index],
        'mean_score': mean_scores.values,
        'rank': range(1, len(mean_scores) + 1),
        'percentile': [(i / len(mean_scores)) * 100 for i in range(1, len(mean_scores) + 1)]
    })
    rankings['gene_upper'] = rankings['gene'].str.upper()
    return rankings, len(df), None


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
# 绘图配置
# =============================================================================
PLOT_COLORS = {
    'essential': '#b30035',
    'nonessential': '#35b300',
    'interest': '#0062b3',
    'boxplot_fill': '#D4A5A5',
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


def centered_plot(fig, config=None, ratio=(1, 8, 1)):
    """将图表居中显示，两侧留白，改善宽高比"""
    _, center, _ = st.columns(ratio)
    with center:
        st.plotly_chart(fig, use_container_width=True, config=config or PLOT_CONFIG)


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

    fig.add_trace(go.Scatter(
        x=bg_df['rank'], y=bg_df['mean_score'], mode='markers',
        marker=dict(size=3, color=th['plot_scatter_bg']),
        name='All genes',
        hovertemplate='<b>%{text}</b><br>Rank: %{x:,}<br>Score: %{y:.4f}<extra></extra>',
        text=bg_df['gene']
    ))

    fig.add_hline(y=-1, line=dict(dash="dash", color=th['text_muted'], width=1))
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

    # 目标基因：按 rank 排序，交替 top/bottom 避免标签重叠
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
                   font=dict(size=16, family=FONT_FAMILY), x=0.5),
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
        height=600, margin=dict(l=70, r=30, t=60, b=60)
    )
    return apply_theme_to_fig(fig)


def create_lineage_boxplot(lineage_data, genes):
    th = get_theme()
    n_genes = len(genes)
    fig = make_subplots(rows=n_genes, cols=1, shared_xaxes=True,
                        vertical_spacing=0.08,
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

    fig.update_layout(
        title=dict(text='<b>CRISPR Score by Cancer Type</b>',
                   font=dict(size=16, family=FONT_FAMILY), x=0.5),
        height=220 * n_genes + 80, showlegend=False,
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

    fig.add_trace(go.Scatter(
        x=bg_all_df['rank'], y=bg_all_df['mean_score'], mode='markers',
        marker=dict(size=2.5, color=th['plot_scatter_bg']),
        name='All genes',
        hovertemplate='<b>%{text}</b><br>Rank: %{x:,}<br>Score: %{y:.4f}<extra></extra>',
        text=bg_all_df['gene']
    ))
    fig.add_hline(y=-1, line=dict(dash="dash", color=th['text_muted'], width=1))
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
                   font=dict(size=16, family=FONT_FAMILY), x=0.5),
        xaxis=dict(title='Gene Rank', showgrid=False, showline=True, linewidth=1.5,
                   tickformat=',d', ticks='outside', ticklen=5),
        yaxis=dict(title=y_label, showgrid=False, showline=True, linewidth=1.5,
                   tickvals=y_tickvals, ticks='outside', ticklen=5,
                   range=[y_min - 0.1 * y_range, y_max + 0.15 * y_range]),
        legend=dict(yanchor='bottom', y=0.02, xanchor='right', x=0.98,
                    font=dict(size=11), bgcolor=th['plot_bg'],
                    bordercolor=th['border']),
        height=600
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


def render_download_buttons(fig, filename_base: str, key_prefix: str, height: int = 600):
    width = 1000
    export_fig = fig_for_export(fig)
    col1, col2, col3 = st.columns(3)

    try:
        pdf_bytes = export_fig.to_image(format='pdf', width=width, height=height, engine='kaleido')
        with col1:
            st.download_button(
                label=t('download_pdf'), data=pdf_bytes,
                file_name=f"{filename_base}.pdf", mime="application/pdf",
                key=f"{key_prefix}_pdf", use_container_width=True
            )
    except Exception as e:
        with col1:
            st.warning(f"PDF: {str(e)[:60]}")

    try:
        png_bytes = export_fig.to_image(format='png', width=width, height=height,
                                        scale=3, engine='kaleido')
        with col2:
            st.download_button(
                label=t('download_png'), data=png_bytes,
                file_name=f"{filename_base}.png", mime="image/png",
                key=f"{key_prefix}_png", use_container_width=True
            )
    except Exception as e:
        with col2:
            st.warning(f"PNG: {str(e)[:60]}")

    try:
        svg_bytes = export_fig.to_image(format='svg', width=width, height=height, engine='kaleido')
        with col3:
            st.download_button(
                label=t('download_svg'), data=svg_bytes,
                file_name=f"{filename_base}.svg", mime="image/svg+xml",
                key=f"{key_prefix}_svg", use_container_width=True
            )
    except Exception as e:
        with col3:
            st.warning(f"SVG: {str(e)[:60]}")

    st.caption(t('download_hint'))


# =============================================================================
# Citation 渲染
# =============================================================================
def build_citations():
    """生成 APA + BibTeX 两种格式"""
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

if uploaded_file is not None:
    with st.spinner(t('loading_upload')):
        crispr_data = load_uploaded_data(uploaded_file.getvalue())
        st.success(f"{t('loaded')}: {uploaded_file.name}")
        data_loaded = True
elif USE_HUGGINGFACE:
    with st.spinner(t('loading_hf')):
        df_result, success, err = download_from_huggingface(HF_REPO_ID, HF_FILENAME)
        if success:
            crispr_data = df_result
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
df_hash = f"{df.shape}_{hash(tuple(df.columns[:5]))}"
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
    essential_count = (gene_rankings['mean_score'] < -0.5).sum()
    st.metric(t('essential_genes'), f"{essential_count:,}")
with col4:
    st.metric(t('score_range'),
              f"{gene_rankings['mean_score'].min():.2f} ~ {gene_rankings['mean_score'].max():.2f}")

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

            with st.expander(f"📥 {t('export_title')}", expanded=True):
                render_download_buttons(fig, "gene_ranking", "rank_plot", height=export_height)

            with st.expander(f"📋 {t('gene_details')}", expanded=False):
                detail = gene_rankings[gene_rankings['gene'].isin(matched_genes)].sort_values('mean_score').copy()
                detail['Essential'] = detail['mean_score'].apply(lambda x: '🔴 Yes' if x < -0.5 else '⚪ No')
                st.dataframe(detail[['gene', 'rank', 'percentile', 'mean_score', 'Essential']].round(4),
                             use_container_width=True, hide_index=True)
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
                st.info("Only displaying first 8 genes for readability")
                matched = matched[:8]
            lineage_data = get_lineage_data(df, matched)
            if lineage_data is not None:
                fig = create_lineage_boxplot(lineage_data, matched)
                # boxplot 是多子图纵向堆叠，保留全宽
                st.plotly_chart(fig, use_container_width=True, config=PLOT_CONFIG)
                with st.expander(f"📥 {t('export_title')}", expanded=True):
                    box_height = max(220 * len(matched) + 80, 400)
                    render_download_buttons(fig, "lineage_boxplot", "boxplot", height=box_height)
            else:
                st.error("No 'lineage' column found in data")

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
            with st.expander(f"📥 {t('export_title')}", expanded=True):
                render_download_buttons(fig, "multilayer_annotation", "multilayer",
                                         height=export_height)


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
