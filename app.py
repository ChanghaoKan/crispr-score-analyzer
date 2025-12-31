"""
CRISPR Score Analysis Web Application
基因必需性分析交互式平台 - 高质量可视化版本

Author: Kan's Lab
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import io

# =============================================================================
# 页面配置（必须放在最前面）
# =============================================================================
st.set_page_config(
    page_title="CRISPR Score Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# =============================================================================
# ⚠️ 重要配置 - 请在这里填入你的 Google Drive 文件 ID
# =============================================================================
# 分享链接格式: https://drive.google.com/file/d/YOUR_FILE_ID/view?usp=sharing
# 只需要提取 YOUR_FILE_ID 部分填入下方

GOOGLE_DRIVE_FILE_ID = "1NMi9mbF51yJ-DAAskDJY7j6kQqhJsQhV"  # ← 替换为你的文件ID，例如: "1AbCdEfGhIjKlMnOpQrS"

# =============================================================================
# 自定义CSS样式 - Cell Journal 风格
# =============================================================================
st.markdown("""
<style>
    .main-header {
        font-family: 'Helvetica Neue', Arial, sans-serif;
        font-size: 2.2rem;
        font-weight: 600;
        color: #1a1a1a;
        margin-bottom: 0.3rem;
        letter-spacing: -0.5px;
    }
    .sub-header {
        font-family: 'Helvetica Neue', Arial, sans-serif;
        font-size: 1rem;
        color: #666;
        margin-bottom: 1.5rem;
    }
    .gene-tag {
        display: inline-block;
        background-color: #e3f2fd;
        color: #1565c0;
        padding: 0.25rem 0.7rem;
        border-radius: 4px;
        margin: 0.15rem;
        font-size: 0.82rem;
        font-weight: 500;
        font-family: 'Monaco', 'Consolas', monospace;
        border: 1px solid #bbdefb;
    }
    .gene-tag-highlight {
        background-color: #ffebee;
        color: #c62828;
        border-color: #ffcdd2;
    }
    .custom-divider {
        height: 1px;
        background: linear-gradient(90deg, transparent, #ddd, transparent);
        border: none;
        margin: 1.5rem 0;
    }
    .input-section-title {
        font-size: 0.9rem;
        font-weight: 600;
        color: #444;
        margin-bottom: 0.5rem;
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
    .stTabs [data-baseweb="tab"] { border-radius: 6px; padding: 8px 16px; }
</style>
""", unsafe_allow_html=True)

# =============================================================================
# Session State 初始化
# =============================================================================
if 'crispr_data' not in st.session_state:
    st.session_state.crispr_data = None
if 'gene_rankings' not in st.session_state:
    st.session_state.gene_rankings = None
if 'n_cell_lines' not in st.session_state:
    st.session_state.n_cell_lines = 0

# =============================================================================
# 数据加载函数
# =============================================================================
@st.cache_data(show_spinner=False, ttl=86400)
def download_from_gdrive(file_id: str):
    """从Google Drive下载数据（缓存24小时）"""
    import tempfile
    import os
    
    try:
        import gdown
    except ImportError:
        return None, False, "缺少 gdown 库"
    
    url = f"https://drive.google.com/uc?id={file_id}"
    
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as tmp:
            tmp_path = tmp.name
        
        # gdown 自动处理大文件的确认下载
        gdown.download(url, tmp_path, quiet=True, fuzzy=True)
        
        # 检查下载的文件是否是CSV（而非HTML警告页）
        with open(tmp_path, 'r', encoding='utf-8', errors='ignore') as f:
            first_line = f.readline()
            if '<html' in first_line.lower() or '<!doctype' in first_line.lower():
                os.unlink(tmp_path)
                return None, False, "下载失败：可能是文件权限问题，请确保设置为'任何人都可以查看'"
        
        df = pd.read_csv(tmp_path)
        os.unlink(tmp_path)
        return df, True, None
    except Exception as e:
        return None, False, f"下载错误: {str(e)}"

@st.cache_data(show_spinner=False)
def load_uploaded_data(file_content):
    """加载用户上传的数据"""
    return pd.read_csv(io.StringIO(file_content.decode('utf-8')))

def extract_gene_name(col_name: str) -> str:
    """从列名提取基因名，支持 'MYC (4609)' 格式"""
    import re
    # 匹配 "GENE (ID)" 或 "GENE (ENTREZ_ID)" 格式
    match = re.match(r'^([A-Za-z0-9_.-]+)\s*\(', str(col_name))
    if match:
        return match.group(1)
    return str(col_name)

@st.cache_data(show_spinner=False)
def compute_gene_rankings(df: pd.DataFrame):
    """计算基因排名"""
    # 识别基因列（数值型，且看起来像CRISPR score）
    gene_cols = []
    meta_cols = []
    
    for col in df.columns:
        # 尝试转换为数值类型
        try:
            # 检查是否是数值列
            if pd.api.types.is_numeric_dtype(df[col]):
                sample_vals = df[col].dropna()
                if len(sample_vals) > 10:  # 至少有一些数据
                    mean_val = sample_vals.mean()
                    std_val = sample_vals.std()
                    # CRISPR scores 通常在 -3 到 1 之间，标准差 > 0
                    if -5 < mean_val < 2 and std_val > 0.01:
                        gene_cols.append(col)
                        continue
        except:
            pass
        meta_cols.append(col)
    
    if not gene_cols:
        return None, 0, "未找到符合CRISPR score特征的数值列"
    
    # 计算每个基因的平均score
    mean_scores = df[gene_cols].mean().sort_values()
    
    # 创建排名DataFrame，提取纯基因名
    rankings = pd.DataFrame({
        'gene_raw': mean_scores.index,  # 原始列名
        'gene': [extract_gene_name(col) for col in mean_scores.index],  # 提取的基因名
        'mean_score': mean_scores.values,
        'rank': range(1, len(mean_scores) + 1),
        'percentile': [(i / len(mean_scores)) * 100 for i in range(1, len(mean_scores) + 1)]
    })
    
    # 创建大写映射用于匹配
    rankings['gene_upper'] = rankings['gene'].str.upper()
    
    return rankings, len(df), None

def filter_genes_by_list(gene_rank_df: pd.DataFrame, gene_list: list):
    """根据基因列表筛选（大小写不敏感）"""
    gene_list_upper = [g.upper() for g in gene_list]
    matched_mask = gene_rank_df['gene_upper'].isin(gene_list_upper)
    matched_genes = gene_rank_df[matched_mask]['gene'].tolist()
    
    matched_upper = set(gene_rank_df[matched_mask]['gene_upper'])
    not_found = [g for g in gene_list if g.upper() not in matched_upper]
    
    return matched_genes, not_found

def get_lineage_data(df: pd.DataFrame, genes: list):
    """获取lineage数据用于boxplot"""
    lineage_col = None
    for col in df.columns:
        if 'lineage' in col.lower() and 'sub' not in col.lower():
            lineage_col = col
            break
    
    if lineage_col is None:
        return None
    
    # 找到实际的基因列名
    df_cols_upper = {c.upper(): c for c in df.columns}
    actual_genes = [df_cols_upper.get(g.upper(), g) for g in genes if g.upper() in df_cols_upper]
    
    if not actual_genes:
        return None
    
    result = []
    for gene in actual_genes:
        if gene in df.columns:
            temp = df[[lineage_col, gene]].copy()
            temp.columns = ['lineage', 'crispr_score']
            temp['gene'] = gene
            result.append(temp)
    
    if result:
        return pd.concat(result, ignore_index=True)
    return None

# =============================================================================
# 可视化函数 - Cell Journal 风格
# =============================================================================
COLORS = {
    'essential': '#b30035',
    'nonessential': '#35b300',
    'interest': '#0062b3',
    'background': 'rgba(200, 200, 200, 0.4)',
    'boxplot_fill': '#D4A5A5',
}
FONT_FAMILY = "Helvetica Neue, Arial, sans-serif"

def create_rank_plot(gene_rank_df, genes_of_interest, essential_gene='MYC', 
                     nonessential_gene='PTEN', n_cell_lines=0, show_labels=True, point_size=4):
    """创建基因排名散点图"""
    fig = go.Figure()
    
    y_min = gene_rank_df['mean_score'].min()
    y_max = gene_rank_df['mean_score'].max()
    y_range = y_max - y_min
    
    # 背景散点
    fig.add_trace(go.Scattergl(
        x=gene_rank_df['rank'], y=gene_rank_df['mean_score'],
        mode='markers',
        marker=dict(size=3, color=COLORS['background']),
        name='All genes',
        hovertemplate='<b>%{text}</b><br>Rank: %{x:,}<br>Score: %{y:.4f}<extra></extra>',
        text=gene_rank_df['gene']
    ))
    
    # 参考线
    fig.add_hline(y=-1, line=dict(dash="dash", color="rgba(128,128,128,0.5)", width=1))
    fig.add_hline(y=0, line=dict(color="black", width=0.8))
    
    # 高亮基因
    interest_df = gene_rank_df[gene_rank_df['gene'].isin(genes_of_interest)]
    if len(interest_df) > 0:
        fig.add_trace(go.Scatter(
            x=interest_df['rank'], y=interest_df['mean_score'],
            mode='markers+text' if show_labels else 'markers',
            marker=dict(size=point_size*2.5, color=COLORS['interest'], opacity=0.9,
                       line=dict(width=0.5, color='white')),
            text=interest_df['gene'] if show_labels else None,
            textposition='top center',
            textfont=dict(size=10, color=COLORS['interest'], family=FONT_FAMILY),
            name='Genes of interest',
            hovertemplate='<b>%{text}</b><br>Rank: %{x:,}<br>Score: %{y:.4f}<br>Percentile: %{customdata:.1f}%<extra></extra>',
            customdata=interest_df['percentile']
        ))
    
    # Essential参考基因
    ess_df = gene_rank_df[gene_rank_df['gene_upper'] == essential_gene.upper()]
    if len(ess_df) > 0:
        fig.add_trace(go.Scatter(
            x=ess_df['rank'], y=ess_df['mean_score'],
            mode='markers+text', marker=dict(size=point_size*2.2, color=COLORS['essential'], symbol='diamond'),
            text=[essential_gene], textposition='bottom center',
            textfont=dict(size=10, color=COLORS['essential'], family=FONT_FAMILY),
            name=f'Essential ({essential_gene})',
            hovertemplate=f'<b>{essential_gene}</b><br>Rank: %{{x:,}}<br>Score: %{{y:.4f}}<extra></extra>'
        ))
    
    # Non-essential参考基因
    noness_df = gene_rank_df[gene_rank_df['gene_upper'] == nonessential_gene.upper()]
    if len(noness_df) > 0:
        fig.add_trace(go.Scatter(
            x=noness_df['rank'], y=noness_df['mean_score'],
            mode='markers+text', marker=dict(size=point_size*2.2, color=COLORS['nonessential'], symbol='diamond'),
            text=[nonessential_gene], textposition='top center',
            textfont=dict(size=10, color=COLORS['nonessential'], family=FONT_FAMILY),
            name=f'Non-essential ({nonessential_gene})',
            hovertemplate=f'<b>{nonessential_gene}</b><br>Rank: %{{x:,}}<br>Score: %{{y:.4f}}<extra></extra>'
        ))
    
    # 布局
    y_label = f"Mean CRISPR Score<br><span style='font-size:11px'>({n_cell_lines} cell lines)</span>" if n_cell_lines > 0 else "Mean CRISPR Score"
    y_tickvals = np.arange(np.floor(y_min/0.5)*0.5, np.ceil(y_max/0.5)*0.5 + 0.5, 0.5)
    
    fig.update_layout(
        title=dict(text='<b>Gene Dependency Ranking</b>', font=dict(size=16, family=FONT_FAMILY), x=0.5),
        xaxis=dict(title='Gene Rank', showgrid=False, showline=True, linewidth=1, linecolor='black',
                   tickformat=',d', ticks='outside', ticklen=5, range=[0, len(gene_rank_df)*1.02]),
        yaxis=dict(title=y_label, showgrid=False, showline=True, linewidth=1, linecolor='black',
                   tickvals=y_tickvals, ticks='outside', ticklen=5, range=[y_min-0.1*y_range, y_max+0.15*y_range]),
        legend=dict(orientation='v', yanchor='bottom', y=0.02, xanchor='right', x=0.98,
                   font=dict(size=9), bgcolor='rgba(255,255,255,0.9)', borderwidth=1),
        plot_bgcolor='white', paper_bgcolor='white', height=550, margin=dict(l=70, r=30, t=60, b=60)
    )
    return fig

def create_lineage_boxplot(lineage_data, genes):
    """创建Lineage Boxplot"""
    n_genes = len(genes)
    fig = make_subplots(rows=n_genes, cols=1, shared_xaxes=True, vertical_spacing=0.08,
                        subplot_titles=[f'<i>{g}</i>' for g in genes])
    
    lineages = sorted(lineage_data['lineage'].unique())
    
    for i, gene in enumerate(genes, 1):
        gene_data = lineage_data[lineage_data['gene'] == gene]
        fig.add_trace(go.Box(
            x=gene_data['lineage'], y=gene_data['crispr_score'], name=gene,
            marker=dict(color=COLORS['boxplot_fill']), line=dict(color='black', width=1),
            fillcolor=COLORS['boxplot_fill'], showlegend=False, boxpoints=False
        ), row=i, col=1)
        fig.add_hline(y=0, line=dict(dash="dot", color="black", width=1), row=i, col=1)
    
    fig.update_layout(
        title=dict(text='<b>CRISPR Score by Cancer Type</b>', font=dict(size=16, family=FONT_FAMILY), x=0.5),
        height=220*n_genes+80, plot_bgcolor='white', paper_bgcolor='white', showlegend=False
    )
    fig.update_xaxes(tickangle=45, categoryarray=lineages, showline=True, linecolor='black')
    fig.update_yaxes(title_text='CRISPR Score', showgrid=False, showline=True, linecolor='black', dtick=0.5)
    
    for i in range(1, n_genes):
        fig.update_xaxes(showticklabels=False, row=i, col=1)
    
    return fig

def create_multilayer_rank_plot(gene_rank_df, background_genes, highlight_genes, 
                                 bg_color='#7FB3D5', hl_color='#E74C3C',
                                 essential_gene='MYC', nonessential_gene='PTEN', 
                                 n_cell_lines=0, show_labels=True):
    """创建多层标注排名图"""
    fig = go.Figure()
    
    y_min = gene_rank_df['mean_score'].min()
    y_max = gene_rank_df['mean_score'].max()
    y_range = y_max - y_min
    
    # 背景
    fig.add_trace(go.Scattergl(
        x=gene_rank_df['rank'], y=gene_rank_df['mean_score'],
        mode='markers', marker=dict(size=2.5, color='rgba(200,200,200,0.35)'),
        name='All genes', hovertemplate='<b>%{text}</b><br>Rank: %{x:,}<br>Score: %{y:.4f}<extra></extra>',
        text=gene_rank_df['gene']
    ))
    
    fig.add_hline(y=-1, line=dict(dash="dash", color="rgba(128,128,128,0.5)", width=1))
    fig.add_hline(y=0, line=dict(color="black", width=0.8))
    
    # 背景基因集
    bg_only = [g for g in background_genes if g not in highlight_genes]
    bg_df = gene_rank_df[gene_rank_df['gene'].isin(bg_only)]
    if len(bg_df) > 0:
        fig.add_trace(go.Scatter(
            x=bg_df['rank'], y=bg_df['mean_score'], mode='markers',
            marker=dict(size=7, color=bg_color, opacity=0.65),
            name=f'Gene set (n={len(bg_df)})',
            hovertemplate='<b>%{text}</b><br>Rank: %{x:,}<br>Score: %{y:.4f}<extra></extra>',
            text=bg_df['gene']
        ))
    
    # 高亮基因
    hl_df = gene_rank_df[gene_rank_df['gene'].isin(highlight_genes)]
    if len(hl_df) > 0:
        fig.add_trace(go.Scatter(
            x=hl_df['rank'], y=hl_df['mean_score'],
            mode='markers+text' if show_labels else 'markers',
            marker=dict(size=11, color=hl_color, line=dict(width=1.5, color='black')),
            text=hl_df['gene'] if show_labels else None, textposition='top center',
            textfont=dict(size=10, color='#333', family=FONT_FAMILY),
            name=f'Highlight (n={len(hl_df)})',
            hovertemplate='<b>%{text}</b><br>Rank: %{x:,}<br>Score: %{y:.4f}<br>Percentile: %{customdata:.1f}%<extra></extra>',
            customdata=hl_df['percentile']
        ))
    
    # 参考基因
    ess_df = gene_rank_df[gene_rank_df['gene_upper'] == essential_gene.upper()]
    if len(ess_df) > 0:
        fig.add_trace(go.Scatter(
            x=ess_df['rank'], y=ess_df['mean_score'], mode='markers+text',
            marker=dict(size=9, color=COLORS['essential'], symbol='diamond'),
            text=[essential_gene], textposition='bottom center',
            textfont=dict(size=10, color=COLORS['essential']), name=f'Essential ({essential_gene})'
        ))
    
    noness_df = gene_rank_df[gene_rank_df['gene_upper'] == nonessential_gene.upper()]
    if len(noness_df) > 0:
        fig.add_trace(go.Scatter(
            x=noness_df['rank'], y=noness_df['mean_score'], mode='markers+text',
            marker=dict(size=9, color=COLORS['nonessential'], symbol='diamond'),
            text=[nonessential_gene], textposition='top center',
            textfont=dict(size=10, color=COLORS['nonessential']), name=f'Non-essential ({nonessential_gene})'
        ))
    
    y_label = f"Mean CRISPR Score<br><span style='font-size:11px'>({n_cell_lines} cell lines)</span>" if n_cell_lines > 0 else "Mean CRISPR Score"
    y_tickvals = np.arange(np.floor(y_min/0.5)*0.5, np.ceil(y_max/0.5)*0.5 + 0.5, 0.5)
    
    fig.update_layout(
        title=dict(text='<b>Multi-layer Gene Annotation</b>', font=dict(size=16, family=FONT_FAMILY), x=0.5),
        xaxis=dict(title='Gene Rank', showgrid=False, showline=True, linecolor='black', tickformat=',d'),
        yaxis=dict(title=y_label, showgrid=False, showline=True, linecolor='black', tickvals=y_tickvals,
                   range=[y_min-0.1*y_range, y_max+0.15*y_range]),
        legend=dict(yanchor='bottom', y=0.02, xanchor='right', x=0.98, font=dict(size=9)),
        plot_bgcolor='white', paper_bgcolor='white', height=550
    )
    return fig

# =============================================================================
# 侧边栏
# =============================================================================
with st.sidebar:
    st.markdown("## 📁 数据来源")
    
    has_gdrive = GOOGLE_DRIVE_FILE_ID is not None
    
    if has_gdrive:
        st.info("☁️ 云端数据 (Google Drive)")
        st.caption("首次加载需下载，之后秒开")
    else:
        st.warning("⚠️ 未配置数据源")
        st.caption("请在 app.py 中配置 GOOGLE_DRIVE_FILE_ID")
    
    with st.expander("📤 上传自定义数据（可选）"):
        uploaded_file = st.file_uploader("上传 CRISPR Score CSV", type=['csv'])
    
    st.markdown("---")
    
    st.markdown("## ⚙️ 参考基因")
    col1, col2 = st.columns(2)
    with col1:
        essential_gene = st.text_input("Essential", value="MYC")
    with col2:
        nonessential_gene = st.text_input("Non-essential", value="PTEN")
    
    st.markdown("---")
    
    st.markdown("## 🎨 显示设置")
    show_labels = st.checkbox("显示基因名标签", value=True)
    point_size = st.slider("点大小", 2, 8, 4)

# =============================================================================
# 数据加载
# =============================================================================
data_loaded = False

# 优先使用上传的数据
if uploaded_file is not None:
    with st.spinner("正在加载上传的数据..."):
        st.session_state.crispr_data = load_uploaded_data(uploaded_file.getvalue())
        st.success(f"✅ 已加载: {uploaded_file.name}")
        data_loaded = True

# 否则从 Google Drive 下载
elif GOOGLE_DRIVE_FILE_ID and st.session_state.crispr_data is None:
    with st.spinner("🔄 首次加载，正在从云端下载数据（约1-2分钟）..."):
        gdrive_df, success, gdrive_error = download_from_gdrive(GOOGLE_DRIVE_FILE_ID)
        if success:
            st.session_state.crispr_data = gdrive_df
            st.success(f"✅ 数据加载成功！后续访问将秒开")
            data_loaded = True
        else:
            st.error(f"❌ 下载失败：{gdrive_error}")

elif st.session_state.crispr_data is not None:
    data_loaded = True

# =============================================================================
# 主界面
# =============================================================================
st.markdown('<h1 class="main-header">🧬 CRISPR Score Analyzer</h1>', unsafe_allow_html=True)
st.markdown('<p class="sub-header">基于 DepMap 数据的基因必需性分析平台</p>', unsafe_allow_html=True)

if not data_loaded:
    st.warning("⚠️ 请配置数据源或上传数据文件")
    st.markdown("""
    ### 配置方法
    
    1. 将你的 DepMap 数据上传到 Google Drive
    2. 设置为"知道链接的任何人可查看"
    3. 复制分享链接中的文件ID
    4. 在 `app.py` 第 27 行填入文件ID
    
    ```python
    GOOGLE_DRIVE_FILE_ID = "你的文件ID"
    ```
    """)
    st.stop()

# 计算基因排名
df = st.session_state.crispr_data
gene_rankings, n_cell_lines, error_msg = compute_gene_rankings(df)
st.session_state.gene_rankings = gene_rankings
st.session_state.n_cell_lines = n_cell_lines

if gene_rankings is None:
    st.error(f"❌ 无法识别基因列：{error_msg}")
    # 显示调试信息
    with st.expander("🔍 数据诊断信息"):
        st.write("**前10列名称：**", list(df.columns[:10]))
        st.write("**数据形状：**", df.shape)
        st.write("**各列数据类型：**")
        st.dataframe(df.dtypes.head(15).to_frame("dtype"))
    st.stop()

# 显示数据概览
col1, col2, col3, col4 = st.columns(4)
with col1:
    st.metric("细胞系", f"{n_cell_lines:,}")
with col2:
    st.metric("基因数", f"{len(gene_rankings):,}")
with col3:
    essential_count = (gene_rankings['mean_score'] < -0.5).sum()
    st.metric("Essential 基因", f"{essential_count:,}")
with col4:
    st.metric("Score 范围", f"{gene_rankings['mean_score'].min():.2f} ~ {gene_rankings['mean_score'].max():.2f}")

st.markdown('<div class="custom-divider"></div>', unsafe_allow_html=True)

# =============================================================================
# 分析标签页
# =============================================================================
tab1, tab2, tab3 = st.tabs(["📊 基因排名图", "📦 Lineage Boxplot", "🎯 多层标注"])

# ----------- Tab 1: 基因排名图 -----------
with tab1:
    st.markdown("### 基因必需性排名")
    st.markdown("在全基因组CRISPR筛选数据中定位您关注的基因。")
    
    st.markdown('<p class="input-section-title">📝 输入目标基因</p>', unsafe_allow_html=True)
    
    input_method = st.radio("选择输入方式", ["直接输入基因名", "上传基因列表文件"], 
                           horizontal=True, label_visibility="collapsed")
    
    genes_of_interest = []
    
    if input_method == "直接输入基因名":
        gene_input = st.text_area("输入基因列表", value="KIF18A\nE2F1\nE2F8\nCCNB1\nCDK1\nAURKA",
                                  height=150, help="每行一个基因名，或用逗号分隔")
        genes_of_interest = [g.strip() for g in gene_input.replace(',', '\n').replace(' ', '\n').split('\n') if g.strip()]
    else:
        uploaded_genelist = st.file_uploader("上传基因列表 (CSV/TXT)", type=['csv', 'txt'], key="genelist1")
        if uploaded_genelist:
            content = uploaded_genelist.getvalue().decode('utf-8')
            if uploaded_genelist.name.endswith('.csv'):
                genes_of_interest = pd.read_csv(io.StringIO(content)).iloc[:, 0].dropna().astype(str).tolist()
            else:
                genes_of_interest = [g.strip() for g in content.split('\n') if g.strip()]
            st.success(f"✓ 已读取 {len(genes_of_interest)} 个基因")
    
    if genes_of_interest:
        matched_genes, not_found = filter_genes_by_list(gene_rankings, genes_of_interest)
        
        col1, col2 = st.columns([3, 1])
        with col1:
            if matched_genes:
                st.markdown("**✓ 匹配成功:**")
                st.markdown(' '.join([f'<span class="gene-tag">{g}</span>' for g in matched_genes]), unsafe_allow_html=True)
        with col2:
            if not_found:
                with st.expander(f"⚠️ 未找到 ({len(not_found)})"):
                    st.write(", ".join(not_found))
        
        if matched_genes:
            fig = create_rank_plot(gene_rankings, matched_genes, essential_gene, nonessential_gene, 
                                  n_cell_lines, show_labels, point_size)
            st.plotly_chart(fig, use_container_width=True)
            
            with st.expander("📋 基因详细信息", expanded=True):
                detail = gene_rankings[gene_rankings['gene'].isin(matched_genes)].sort_values('mean_score').copy()
                detail['Essential'] = detail['mean_score'].apply(lambda x: '🔴 Yes' if x < -0.5 else '⚪ No')
                st.dataframe(detail[['gene', 'rank', 'percentile', 'mean_score', 'Essential']].round(4), 
                           use_container_width=True, hide_index=True)

# ----------- Tab 2: Lineage Boxplot -----------
with tab2:
    st.markdown("### 按癌症类型的CRISPR Score分布")
    
    st.markdown('<p class="input-section-title">📝 输入目标基因</p>', unsafe_allow_html=True)
    
    input_method2 = st.radio("选择输入方式", ["直接输入基因名", "上传基因列表文件"], 
                            horizontal=True, key="box_input", label_visibility="collapsed")
    
    genes_for_box = []
    if input_method2 == "直接输入基因名":
        gene_input2 = st.text_area("输入基因", value="KIF18A\nE2F1\nMYC", height=120, key="box_text")
        genes_for_box = [g.strip() for g in gene_input2.replace(',', '\n').split('\n') if g.strip()]
    else:
        uploaded2 = st.file_uploader("上传基因列表", type=['csv', 'txt'], key="box_file")
        if uploaded2:
            content = uploaded2.getvalue().decode('utf-8')
            if uploaded2.name.endswith('.csv'):
                genes_for_box = pd.read_csv(io.StringIO(content)).iloc[:, 0].dropna().astype(str).tolist()
            else:
                genes_for_box = [g.strip() for g in content.split('\n') if g.strip()]
    
    if genes_for_box:
        matched, not_found = filter_genes_by_list(gene_rankings, genes_for_box)
        if not_found:
            st.warning(f"未找到: {', '.join(not_found)}")
        if matched:
            if len(matched) > 8:
                st.info("为保证效果，仅展示前8个基因")
                matched = matched[:8]
            
            lineage_data = get_lineage_data(df, matched)
            if lineage_data is not None:
                fig = create_lineage_boxplot(lineage_data, matched)
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.error("未找到 lineage 列")

# ----------- Tab 3: 多层标注 -----------
with tab3:
    st.markdown("### 多层基因标注")
    
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**背景基因集**")
        bg_input = st.text_area("背景基因", value="CDK1\nCDK2\nCCNB1\nCCND1\nCCNE1", height=150, key="bg")
        bg_color = st.color_picker("背景颜色", "#7FB3D5", key="bg_color")
    with col2:
        st.markdown("**高亮基因**")
        hl_input = st.text_area("高亮基因", value="KIF18A\nE2F8", height=150, key="hl")
        hl_color = st.color_picker("高亮颜色", "#E74C3C", key="hl_color")
    
    bg_genes = [g.strip() for g in bg_input.replace(',', '\n').split('\n') if g.strip()]
    hl_genes = [g.strip() for g in hl_input.replace(',', '\n').split('\n') if g.strip()]
    
    if bg_genes or hl_genes:
        bg_matched, _ = filter_genes_by_list(gene_rankings, bg_genes)
        hl_matched, _ = filter_genes_by_list(gene_rankings, hl_genes)
        
        st.markdown(f"背景: {len(bg_matched)} 个 | 高亮: {len(hl_matched)} 个")
        
        if bg_matched or hl_matched:
            fig = create_multilayer_rank_plot(gene_rankings, bg_matched, hl_matched,
                                              bg_color, hl_color, essential_gene, nonessential_gene,
                                              n_cell_lines, show_labels)
            st.plotly_chart(fig, use_container_width=True)

# =============================================================================
# 致谢与页脚
# =============================================================================
st.markdown("---")

st.markdown("""
<div style="background: linear-gradient(135deg, #f5f7fa 0%, #e4e8ec 100%); 
            padding: 1.5rem; border-radius: 12px; margin-bottom: 1rem;">
    <h4 style="margin: 0 0 1rem 0; color: #333;">🙏 Acknowledgements</h4>
    <div style="display: flex; flex-wrap: wrap; gap: 1.5rem;">
        <div style="flex: 1; min-width: 250px;">
            <p style="margin: 0; color: #555; font-size: 0.9rem;">
                <strong>Data Source</strong><br>
                <a href="https://depmap.org" target="_blank" style="color: #1565c0;">DepMap Portal (Broad Institute)</a><br>
                <span style="font-size: 0.8rem; color: #777;">CRISPR Chronos dependency scores</span>
            </p>
        </div>
        <div style="flex: 1; min-width: 250px;">
            <p style="margin: 0; color: #555; font-size: 0.9rem;">
                <strong>Development Assistance</strong><br>
                <a href="https://www.anthropic.com/claude" target="_blank" style="color: #1565c0;">Claude (Anthropic)</a><br>
                <span style="font-size: 0.8rem; color: #777;">AI-assisted development</span>
            </p>
        </div>
    </div>
</div>
""", unsafe_allow_html=True)

with st.expander("📚 Citation"):
    st.markdown("""
    **DepMap:** Tsherniak, A., et al. Defining a Cancer Dependency Map. *Cell* 170, 564-576 (2017).
    
    **Portal:** https://depmap.org/portal/
    """)

st.markdown('<div style="text-align:center; color:#999; font-size:0.8rem; padding:1rem;">CRISPR Score Analyzer v2.0 | Developed by Kan\'s Lab</div>', unsafe_allow_html=True)
