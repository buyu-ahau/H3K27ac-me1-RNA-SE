# RNA-seq与组蛋白修饰热图分析

## 🔬 项目简介
这个项目提供了一个完整的R脚本，用于分析RNA-seq数据和组蛋白修饰数据（H3K27ac和H3K4me1），并生成完美对齐的热图可视化。

## 🎯 主要特点
- 🔧 **自动数据读取**：支持多种文件格式
- 🧠 **智能列识别**：自动识别基因名和表达量列  
- 📊 **完美数据对齐**：确保三个数据集基因顺序一致
- 🧹 **基因名清理**：去掉前缀和重复后缀
- 📈 **聚类分析**：样本聚类（上方），基因聚类隐藏

## 🎨 热图特点
- ✅ 上方样本聚类树
- ✅ 正确分组标注（Mhp=红色实验组，Con=蓝色对照组）
- ✅ 三图完美对齐
- ✅ 清晰基因名
- ✅ 高质量输出

## 📦 安装要求
### R包依赖
```r
install.packages(c("pheatmap", "RColorBrewer", "gridExtra", "grid", "dplyr"))
```

## 📁 文件结构
```
rnaseq-heatmap-analysis/
├── rnaseq_heatmap_analysis.R    # 主分析脚本
├── README.md                    # 项目说明文档
├── data/                        # 数据文件夹
│   ├── 0-28表达量.txt
│   ├── Final_H3K27ac_Signal_on_DiffSE50.tsv
│   └── Final_H3K4me1_Signal_on_DiffSE50.tsv
└── results/                     # 输出结果文件夹
    ├── H3K27ac_Final.pdf
    ├── H3K4me1_Final.pdf
    ├── RNAseq_Final.pdf
    ├── Triple_Heatmap_Final_Horizontal.pdf
    └── ...
```

## 🚀 使用方法
1. **准备数据文件**：将以下三个数据文件放在工作目录中
   - `0-28表达量.txt` (RNA-seq数据)
   - `Final_H3K27ac_Signal_on_DiffSE50.tsv`
   - `Final_H3K4me1_Signal_on_DiffSE50.tsv`

2. **运行分析**：
   ```r
   source("rnaseq_heatmap_analysis.R")
   ```

3. **查看结果**：检查生成的PDF热图文件

## 📊 输入文件说明
- **0-28表达量.txt**: RNA-seq表达量数据，包含基因名和D0、D28时间点的FPKM值
- **Final_H3K27ac_Signal_on_DiffSE50.tsv**: H3K27ac组蛋白修饰信号数据
- **Final_H3K4me1_Signal_on_DiffSE50.tsv**: H3K4me1组蛋白修饰信号数据

## 📈 输出文件说明
- **单独热图**：`H3K27ac_Final.pdf`, `H3K4me1_Final.pdf`, `RNAseq_Final.pdf`
- **组合热图**：`Triple_Heatmap_Final_Horizontal.pdf` (推荐)
- **数据文件**：`Gene_Matching_Results.csv`, `Extracted_RNAseq_Data.csv`

## 🔬 分析流程
1. **数据读取**：自动识别文件格式和列名
2. **基因匹配**：基于H3K27ac基因列表提取对应的RNA-seq数据
3. **数据对齐**：确保三个数据集基因顺序完全一致
4. **聚类分析**：进行样本聚类分析
5. **热图生成**：创建专业的热图可视化
6. **结果保存**：输出多种格式的结果文件

## 🎨 可视化特点
- **颜色方案**：蓝-白-红渐变，清晰区分高低表达
- **分组标注**：Mhp(实验组-红色)，Con(对照组-蓝色)
- **聚类树**：上方显示样本间聚类关系
- **基因名**：清理后的简洁基因名（去掉前缀和后缀）

## 👤 作者
**buyu-ahau**

## 📅 更新日期
2025-06-11

## 📄 许可证
MIT License

## 🤝 贡献
欢迎提交Issues和Pull Requests来改进这个项目！

## 📧 联系
如有问题或建议，请通过GitHub Issues联系。
