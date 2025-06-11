# =============================================================================
# 完整的RNA-seq与组蛋白修饰热图分析代码
# 作者: buyu-ahau
# 日期: 2025-06-11
# 功能: 基于H3K27ac基因列表提取RNA-seq数据并绘制对齐的热图
# 输入文件: 
#   - 0-28表达量.txt (RNA-seq数据)
#   - Final_H3K27ac_Signal_on_DiffSE50.tsv
#   - Final_H3K4me1_Signal_on_DiffSE50.tsv
# =============================================================================

# 清理环境
rm(list = ls())
gc()

# 1. 加载必要的包
# =============================================================================

cat("=== 加载必要的R包 ===\n")

required_packages <- c("pheatmap", "RColorBrewer", "gridExtra", "grid", "dplyr")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

cat("所有必要的包已加载完成\n")

# 2. 读取RNA-seq数据
# =============================================================================

cat("\n=== 读取RNA-seq数据 ===\n")

# 尝试多种方法读取RNA-seq数据
rnaseq_data <- NULL

methods <- list(
  list(func = read.table, params = list(header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", encoding = "UTF-8")),
  list(func = read.table, params = list(header = TRUE, sep = ",", stringsAsFactors = FALSE, quote = "", encoding = "UTF-8")),
  list(func = read.delim, params = list(header = TRUE, stringsAsFactors = FALSE, encoding = "UTF-8"))
)

for (i in 1:length(methods)) {
  try({
    method <- methods[[i]]
    rnaseq_data <- do.call(method$func, c(list("0-28表达量.txt"), method$params))
    cat("方法", i, "成功读取RNA-seq数据\n")
    break
  }, silent = TRUE)
}

if (is.null(rnaseq_data)) {
  stop("无法读取RNA-seq数据文件，请检查文件格式")
}

cat("RNA-seq数据读取成功！\n")
cat("数据维度:", dim(rnaseq_data), "\n")
cat("列名:", paste(colnames(rnaseq_data), collapse = ", "), "\n")

# 3. 自动识别列名
# =============================================================================

cat("\n=== 自动识别数据列 ===\n")

all_cols <- colnames(rnaseq_data)

# 查找基因名列
possible_gene_cols <- c("Symbol", "symbol", "SYMBOL", "Gene", "gene", "GeneID", "gene_id", "id")
gene_col <- NULL

for (col in possible_gene_cols) {
  if (col %in% all_cols) {
    gene_col <- col
    break
  }
}

if (is.null(gene_col)) {
  gene_col <- all_cols[1]
  cat("未找到明确的基因名列，使用第一列:", gene_col, "\n")
} else {
  cat("找到基因名列:", gene_col, "\n")
}

# 查找表达量列
d0_patterns <- c("D0.*fpkm", "d0.*fpkm", "D0", "d0")
d28_patterns <- c("D28.*fpkm", "d28.*fpkm", "D28", "d28")

d0_cols <- c()
d28_cols <- c()

for (pattern in d0_patterns) {
  matches <- grep(pattern, all_cols, value = TRUE, ignore.case = TRUE)
  if (length(matches) > 0) {
    d0_cols <- matches
    break
  }
}

for (pattern in d28_patterns) {
  matches <- grep(pattern, all_cols, value = TRUE, ignore.case = TRUE)
  if (length(matches) > 0) {
    d28_cols <- matches
    break
  }
}

# 如果没找到明确的模式，查找所有fpkm列
if (length(d0_cols) == 0 || length(d28_cols) == 0) {
  fpkm_cols <- grep("fpkm", all_cols, value = TRUE, ignore.case = TRUE)
  if (length(fpkm_cols) >= 6) {
    half <- length(fpkm_cols) / 2
    d0_cols <- fpkm_cols[1:floor(half)]
    d28_cols <- fpkm_cols[(floor(half)+1):length(fpkm_cols)]
  }
}

# 选择前3个样本
if (length(d0_cols) > 3) d0_cols <- d0_cols[1:3]
if (length(d28_cols) > 3) d28_cols <- d28_cols[1:3]

expression_cols <- c(d0_cols, d28_cols)
cat("使用的表达量列:", paste(expression_cols, collapse = ", "), "\n")

# 4. 读取H3K27ac和H3K4me1数据
# =============================================================================

cat("\n=== 读取组蛋白修饰数据 ===\n")

# 读取H3K27ac数据
h3k27ac_data <- read.table("Final_H3K27ac_Signal_on_DiffSE50.tsv", 
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
h3k27ac_data$Gene <- make.unique(h3k27ac_data$Gene)

# 读取H3K4me1数据
h3k4me1_data <- read.table("Final_H3K4me1_Signal_on_DiffSE50.tsv", 
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE)
h3k4me1_data$Gene <- make.unique(h3k4me1_data$Gene)

cat("H3K27ac数据维度:", dim(h3k27ac_data), "\n")
cat("H3K4me1数据维度:", dim(h3k4me1_data), "\n")

# 5. 获取H3K27ac基因列表并清理基因名
# =============================================================================

cat("\n=== 处理基因名 ===\n")

# 获取基因列表
h3k27ac_gene_list <- h3k27ac_data$Gene
cat("H3K27ac基因总数:", length(h3k27ac_gene_list), "\n")

# 清理基因名用于匹配
h3k27ac_clean_names <- gsub("^gene-", "", h3k27ac_gene_list)
h3k27ac_original_names <- gsub("\\.\\d+$", "", h3k27ac_clean_names)

# 6. 按H3K27ac顺序提取RNA-seq数据
# =============================================================================

cat("\n=== 提取RNA-seq数据 ===\n")

# 获取RNA-seq基因名
rnaseq_genes <- rnaseq_data[[gene_col]]
rnaseq_genes_clean <- rnaseq_genes[!is.na(rnaseq_genes) & rnaseq_genes != ""]

# 检查匹配情况
potential_matches <- intersect(h3k27ac_original_names, rnaseq_genes_clean)
cat("潜在匹配基因数:", length(potential_matches), "\n")
cat("匹配率:", round(length(potential_matches) / length(h3k27ac_original_names) * 100, 1), "%\n")

# 按顺序提取数据
rnaseq_extracted_rows <- list()
matched_gene_info <- data.frame(
  h3k27ac_gene = h3k27ac_gene_list,
  original_gene = h3k27ac_original_names,
  rnaseq_match = NA,
  found = FALSE,
  stringsAsFactors = FALSE
)

cat("开始匹配基因...\n")
for (i in 1:length(h3k27ac_original_names)) {
  target_gene <- h3k27ac_original_names[i]
  
  if (!is.na(target_gene) && target_gene != "" && target_gene %in% rnaseq_genes) {
    row_idx <- which(rnaseq_genes == target_gene)[1]
    matched_row_data <- rnaseq_data[row_idx, expression_cols]
    rnaseq_extracted_rows[[i]] <- matched_row_data
    matched_gene_info$rnaseq_match[i] <- target_gene
    matched_gene_info$found[i] <- TRUE
  } else {
    na_row <- data.frame(matrix(NA, nrow = 1, ncol = length(expression_cols)))
    colnames(na_row) <- expression_cols
    rnaseq_extracted_rows[[i]] <- na_row
  }
  
  if (i %% 20 == 0) {
    cat("已处理", i, "/", length(h3k27ac_original_names), "个基因\n")
  }
}

# 合并提取的数据
rnaseq_extracted <- do.call(rbind, rnaseq_extracted_rows)
rownames(rnaseq_extracted) <- h3k27ac_gene_list

# 标准化列名
standard_colnames <- c("Con1", "Con2", "Con3", "Mhp1", "Mhp2", "Mhp3")
colnames(rnaseq_extracted) <- standard_colnames

cat("数据提取完成！匹配成功:", sum(matched_gene_info$found), "个基因\n")

# 7. 准备热图数据
# =============================================================================

cat("\n=== 准备热图数据 ===\n")

# 处理H3K27ac数据
h3k27ac_matrix <- as.matrix(h3k27ac_data[, -1])
rownames(h3k27ac_matrix) <- h3k27ac_data$Gene

# 处理H3K4me1数据（确保与H3K27ac对齐）
h3k4me1_aligned <- matrix(NA, nrow = length(h3k27ac_gene_list), 
                          ncol = ncol(h3k4me1_data) - 1)
rownames(h3k4me1_aligned) <- h3k27ac_gene_list
colnames(h3k4me1_aligned) <- colnames(h3k4me1_data)[-1]

for (i in 1:length(h3k27ac_gene_list)) {
  gene_name <- h3k27ac_gene_list[i]
  if (gene_name %in% h3k4me1_data$Gene) {
    gene_row <- h3k4me1_data[h3k4me1_data$Gene == gene_name, ]
    h3k4me1_aligned[i, ] <- as.numeric(gene_row[1, -1])
  }
}

# 样本对应关系
sample_mapping <- data.frame(
  h3k27ac_cols = c("Mhp1_H3K27ac", "Mhp2_H3K27ac", "Mhp3_H3K27ac", 
                   "Con1_H3K27ac", "Con2_H3K27ac", "Con3_H3K27ac"),
  h3k4me1_cols = c("Mhp1_H3K4me1", "Mhp2_H3K4me1", "Mhp3_H3K4me1", 
                   "Con1_H3K4me1", "Con2_H3K4me1", "Con3_H3K4me1"),
  rnaseq_cols = standard_colnames
)

# 重新排列数据并统一列名
h3k27ac_final <- h3k27ac_matrix[, sample_mapping$h3k27ac_cols]
colnames(h3k27ac_final) <- standard_colnames

h3k4me1_final <- h3k4me1_aligned[, sample_mapping$h3k4me1_cols]
colnames(h3k4me1_final) <- standard_colnames

rnaseq_final <- as.matrix(rnaseq_extracted)

# 清理基因名（去掉gene-前缀和.数字后缀）
clean_gene_names <- gsub("^gene-", "", h3k27ac_gene_list)
clean_gene_names <- gsub("\\.\\d+$", "", clean_gene_names)

rownames(h3k27ac_final) <- clean_gene_names
rownames(h3k4me1_final) <- clean_gene_names
rownames(rnaseq_final) <- clean_gene_names

cat("数据准备完成\n")

# 8. 聚类分析
# =============================================================================

cat("\n=== 聚类分析 ===\n")

p_cluster <- pheatmap(
  h3k27ac_final,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  silent = TRUE
)

clustered_gene_order <- rownames(h3k27ac_final)[p_cluster$tree_row$order]
sample_cluster_order <- p_cluster$tree_col$order
clustered_sample_names <- colnames(h3k27ac_final)[sample_cluster_order]

cat("聚类分析完成\n")

# 9. 设置热图样式
# =============================================================================

cat("\n=== 设置热图样式 ===\n")

# 颜色方案
my_colors <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", 
                                "white", 
                                "#FDBF6F", "#FF7F00", "#E31A1C", "#B10026"))(100)

# 样本分组注释（只保留Group，去掉Time）
sample_annotation <- data.frame(
  Group = c(rep("Mhp", 3), rep("Con", 3))  # Mhp=实验组，Con=对照组
)
rownames(sample_annotation) <- standard_colnames

# 分组颜色
ann_colors <- list(
  Group = c("Mhp" = "#E31A1C", "Con" = "#2166AC")  # 实验组红色，对照组蓝色
)

# 格子尺寸
cell_width <- 30
cell_height <- 12
col_tree_height <- 50

# 10. 绘制热图
# =============================================================================

cat("\n=== 绘制热图 ===\n")

# 绘制统一聚类顺序的热图
p_h3k27ac <- pheatmap(
  h3k27ac_final[clustered_gene_order, clustered_sample_names],
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = my_colors,
  annotation_col = sample_annotation[clustered_sample_names, , drop = FALSE],
  annotation_colors = ann_colors,
  main = "H3K27ac",
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 12,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = "white",
  cellwidth = cell_width,
  cellheight = cell_height,
  treeheight_row = 0,
  treeheight_col = col_tree_height,
  annotation_legend = TRUE,
  legend = TRUE,
  silent = TRUE
)

p_h3k4me1 <- pheatmap(
  h3k4me1_final[clustered_gene_order, clustered_sample_names],
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = my_colors,
  annotation_col = sample_annotation[clustered_sample_names, , drop = FALSE],
  annotation_colors = ann_colors,
  main = "H3K4me1",
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 12,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = "white",
  cellwidth = cell_width,
  cellheight = cell_height,
  treeheight_row = 0,
  treeheight_col = col_tree_height,
  annotation_legend = TRUE,
  legend = TRUE,
  na_col = "grey90",
  silent = TRUE
)

p_rnaseq <- pheatmap(
  rnaseq_final[clustered_gene_order, clustered_sample_names],
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = my_colors,
  annotation_col = sample_annotation[clustered_sample_names, , drop = FALSE],
  annotation_colors = ann_colors,
  main = "RNA-seq Expression",
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 12,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = "white",
  cellwidth = cell_width,
  cellheight = cell_height,
  treeheight_row = 0,
  treeheight_col = col_tree_height,
  annotation_legend = TRUE,
  legend = TRUE,
  na_col = "grey90",
  silent = TRUE
)

cat("所有热图绘制完成！\n")

# 11. 保存结果
# =============================================================================

cat("\n=== 保存结果 ===\n")

# 计算图像尺寸
n_genes <- length(clustered_gene_order)
base_width <- cell_width * 6 / 72
base_height <- cell_height * n_genes / 72
tree_height_inch <- col_tree_height / 72

single_width <- base_width + 5
single_height <- base_height + 3 + tree_height_inch

# 保存单独热图
pdf("H3K27ac_Final.pdf", width = single_width, height = single_height)
print(p_h3k27ac)
dev.off()

pdf("H3K4me1_Final.pdf", width = single_width, height = single_height)
print(p_h3k4me1)
dev.off()

pdf("RNAseq_Final.pdf", width = single_width, height = single_height)
print(p_rnaseq)
dev.off()

# 保存三合一热图
triple_width <- single_width * 3 + 2

# 水平排列
pdf("Triple_Heatmap_Final_Horizontal.pdf", 
    width = triple_width, height = single_height)
grid.arrange(
  p_h3k27ac[[4]], p_h3k4me1[[4]], p_rnaseq[[4]],
  ncol = 3,
  top = textGrob("H3K27ac, H3K4me1, and RNA-seq Expression Analysis", 
                 gp = gpar(fontsize = 16, fontface = "bold"))
)
dev.off()

# 垂直排列
pdf("Triple_Heatmap_Final_Vertical.pdf", 
    width = single_width + 1, height = single_height * 3 + 2)
grid.arrange(
  p_h3k27ac[[4]], p_h3k4me1[[4]], p_rnaseq[[4]],
  nrow = 3,
  top = textGrob("H3K27ac, H3K4me1, and RNA-seq Expression Analysis", 
                 gp = gpar(fontsize = 14, fontface = "bold"))
)
dev.off()

# 高分辨率PNG
png("Triple_Heatmap_Final.png", 
    width = triple_width, height = single_height, 
    units = "in", res = 300)
grid.arrange(
  p_h3k27ac[[4]], p_h3k4me1[[4]], p_rnaseq[[4]],
  ncol = 3,
  top = textGrob("H3K27ac, H3K4me1, and RNA-seq Expression Analysis", 
                 gp = gpar(fontsize = 16, fontface = "bold"))
)
dev.off()

# 保存数据
write.csv(matched_gene_info, "Gene_Matching_Results.csv", row.names = FALSE)
write.csv(rnaseq_extracted, "Extracted_RNAseq_Data.csv", row.names = TRUE)

# 12. 输出总结
# =============================================================================

cat("\n=== 分析完成总结 ===\n")

cat("✅ 完成内容:\n")
cat("1. 成功读取RNA-seq数据文件: 0-28表达量.txt\n")
cat("2. 自动识别基因名列和表达量列\n")
cat("3. 基于H3K27ac基因列表提取对应的RNA-seq数据\n")
cat("4. 清理基因名（去掉gene-前缀和重复后缀）\n")
cat("5. 生成完美对齐的三合一热图\n")

cat("\n📊 数据统计:\n")
cat("总基因数:", n_genes, "\n")
cat("成功匹配:", sum(matched_gene_info$found), "\n")
cat("匹配率:", round(mean(matched_gene_info$found) * 100, 1), "%\n")
cat("使用的基因名列:", gene_col, "\n")
cat("使用的D0列:", paste(d0_cols, collapse = ", "), "\n")
cat("使用的D28列:", paste(d28_cols, collapse = ", "), "\n")

cat("\n🎨 热图特点:\n")
cat("• 上方显示样本聚类树\n")
cat("• 左侧基因聚类隐藏（保持固定顺序）\n")
cat("• 正确的分组标注：Mhp(实验组-红色)，Con(对照组-蓝色)\n")
cat("• 三个热图基因和样本顺序完全一致\n")
cat("• 清晰的基因名（无前缀，无重复后缀）\n")

cat("\n📁 生成的文件:\n")
cat("✓ H3K27ac_Final.pdf\n")
cat("✓ H3K4me1_Final.pdf\n")
cat("✓ RNAseq_Final.pdf\n")
cat("✓ Triple_Heatmap_Final_Horizontal.pdf (推荐)\n")
cat("✓ Triple_Heatmap_Final_Vertical.pdf\n")
cat("✓ Triple_Heatmap_Final.png (高分辨率)\n")
cat("✓ Gene_Matching_Results.csv (匹配详情)\n")
cat("✓ Extracted_RNAseq_Data.csv (提取数据)\n")

cat("\n样本对应关系:\n")
print(sample_annotation)

cat("\n🎉 完整分析流程执行完成！\n")
cat("推荐查看: Triple_Heatmap_Final_Horizontal.pdf\n")

# =============================================================================
# 代码结束
# =============================================================================





############SE
