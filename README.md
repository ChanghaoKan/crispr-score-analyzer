# 🧬 CRISPR Score Analyzer

基于 DepMap 数据的基因必需性分析平台。

## 功能

| 功能 | 说明 |
|-----|------|
| **基因排名图** | 在全基因组中定位目标基因的必需性排名 |
| **Lineage Boxplot** | 查看基因在不同癌症类型中的依赖性差异 |
| **多层标注** | 同时展示基因集背景与高亮基因 |

## 使用方法

1. 在输入框输入基因名（每行一个）或上传基因列表文件
2. 点击对应功能标签页查看可视化结果
3. 图表支持交互：悬停查看详情，可导出为 SVG/PNG

## CRISPR Score 解读

| Score | 含义 |
|-------|------|
| < -1.0 | 强必需基因（敲除致死）|
| -1.0 ~ -0.5 | 中度必需 |
| > -0.5 | 非必需 |

## 致谢

- **数据来源**: [DepMap Portal](https://depmap.org) (Broad Institute)
- **开发辅助**: [Claude](https://www.anthropic.com/claude) (Anthropic)
- **团队支持**: [DengLab](https://www.deng-lab.org/) (Team)
- 
---

Developed by Kan's Lab
