# 🧬 CRISPR Score Analyzer

基于 DepMap 数据的基因必需性分析平台，提供交互式可视化分析。

## ✨ 功能

- **基因排名图**: 全基因组必需性排名，定位目标基因
- **Lineage Boxplot**: 按癌症类型展示CRISPR Score分布
- **多层标注**: 背景基因集 + 高亮基因分层展示

## 🚀 部署步骤

### Step 1: 上传数据到 Google Drive

1. 登录 [Google Drive](https://drive.google.com)
2. 上传你的 DepMap CSV 文件
3. 右键文件 → **共享** → **常规访问** 改为 **"知道链接的任何人"**
4. 点击 **复制链接**

### Step 2: 获取文件 ID

分享链接格式：
```
https://drive.google.com/file/d/1AbCdEfGhIjKlMnOpQrStUvWxYz/view?usp=sharing
                                 └──────────────────────────┘
                                        这部分是 FILE_ID
```

### Step 3: 配置代码

编辑 `app.py` 第 27 行，填入你的文件 ID：

```python
GOOGLE_DRIVE_FILE_ID = "1AbCdEfGhIjKlMnOpQrStUvWxYz"  # 替换为你的ID
```

### Step 4: 推送到 GitHub

```bash
# 如果是新仓库
git init
git add .
git commit -m "Initial commit"
git remote add origin https://github.com/你的用户名/crispr-analyzer.git
git branch -M main
git push -u origin main

# 如果是更新现有仓库
git add .
git commit -m "Update"
git push
```

### Step 5: 部署到 Streamlit Cloud

1. 访问 [share.streamlit.io](https://share.streamlit.io)
2. 用 GitHub 账号登录
3. 点击 **New app**
4. 选择你的仓库、分支 `main`、文件 `app.py`
5. 点击 **Deploy**

等待 2-5 分钟，获得你的公开链接！

---

## 📋 数据格式

DepMap CSV 格式：

| 列名 | 说明 |
|-----|------|
| DepMap_ID | 细胞系ID |
| cell_line_name | 细胞系名称 |
| lineage | 癌症类型 |
| MYC, PTEN, ... | 各基因的CRISPR Score |

**推荐数据**: [DepMap Portal](https://depmap.org/portal/) → CRISPR (DepMap Public, Chronos)

---

## 🙏 致谢

- **DepMap Portal** (Broad Institute) - CRISPR 筛选数据
- **Claude** (Anthropic) - AI 辅助开发

---

Developed by Kan's Lab
