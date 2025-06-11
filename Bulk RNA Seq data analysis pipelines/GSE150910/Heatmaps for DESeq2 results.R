library(pheatmap)
library(dplyr)


### 1) Heatmap for IPF vs Control
# Normalize counts using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE) 
vsd_mat <- assay(vsd)

# Select top 50 significant DEGs (FDR < 0.05 and |log2FC| > 1)
top50_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# Subset normalized expression matrix
heatmap_data <- vsd_mat[top50_genes$Gene_symbol, ]

# Sample annotation
sample_ann <- data.frame(Group = sample_info$condition)
rownames(sample_ann) <- rownames(sample_info)

# Gene annotation
gene_annotation <- data.frame(Regulation = top50_genes$regulation)
rownames(gene_annotation) <- top50_genes$Gene_symbol

# Color palettes
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)
sample_colors <- c("CONTROL" = "grey", "IPF" = "orangered")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  Group = sample_colors,
  Regulation = gene_colors
)

# Plot the heatmap
heatmap <- pheatmap(heatmap_data,
                    scale = "row",
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method = "complete",
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    fontsize = 8,
                    main = "Top 50 DEGs Heatmap: IPF vs Control")

# Save the heatmap
ggsave("IPFvsControl_Top50_DEGs_Heatmap.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)





### 2) Heatmap for CHP vs Control
# Normalize counts using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE) 
vsd_mat <- assay(vsd)

# Select top 50 significant DEGs (FDR < 0.05 and |log2FC| > 1)
top50_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# Subset normalized expression matrix
heatmap_data <- vsd_mat[top50_genes$Gene_symbol, ]

# Sample annotation
sample_ann <- data.frame(Group = sample_info_2$condition)
rownames(sample_ann) <- rownames(sample_info_2)

# Gene annotation
gene_annotation <- data.frame(Regulation = top50_genes$regulation)
rownames(gene_annotation) <- top50_genes$Gene_symbol

# Color palettes
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)
sample_colors <- c("CONTROL" = "grey", "CHP" = "cadetblue")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  Group = sample_colors,
  Regulation = gene_colors
)

# Plot the heatmap
heatmap <- pheatmap(heatmap_data,
                    scale = "row",
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method = "complete",
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    fontsize = 8,
                    main = "Top 50 DEGs Heatmap: HP vs Control")

# Save the heatmap
ggsave("CHPvsControl_Top50_DEGs_Heatmap.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)



### 3) Heatmap for CHP vs IPF
# Normalize counts using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE) 
vsd_mat <- assay(vsd)

# Select top 50 significant DEGs (FDR < 0.05 and |log2FC| > 1)
top50_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# Subset normalized expression matrix
heatmap_data <- vsd_mat[top50_genes$Gene_symbol, ]

# Sample annotation
sample_ann <- data.frame(Group = sample_info_3$condition)
rownames(sample_ann) <- rownames(sample_info_3)

# Gene annotation
gene_annotation <- data.frame(Regulation = top50_genes$regulation)
rownames(gene_annotation) <- top50_genes$Gene_symbol

# Color palettes
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)
sample_colors <- c("IPF" = "orangered", "CHP" = "cadetblue")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  Group = sample_colors,
  Regulation = gene_colors
)

# Plot the heatmap
heatmap <- pheatmap(heatmap_data,
                    scale = "row",
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method = "complete",
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    fontsize = 8,
                    main = "Top 50 DEGs Heatmap: HP vs IPF")

# Save the heatmap
ggsave("CHPvsIPF_Top50_DEGs_Heatmap.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)




### 3) Heatmap for ILD vs Control
# Normalize counts using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE) 
vsd_mat <- assay(vsd)

# Select top 50 significant DEGs (FDR < 0.05 and |log2FC| > 1)
top50_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# Subset normalized expression matrix
heatmap_data <- vsd_mat[top50_genes$Gene_symbol, ]

# Sample annotation
sample_ann <- data.frame(Group = sample_info.ild$condition)
rownames(sample_ann) <- rownames(sample_info.ild)

# Gene annotation
gene_annotation <- data.frame(Regulation = top50_genes$regulation)
rownames(gene_annotation) <- top50_genes$Gene_symbol

# Color palettes
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)
sample_colors <- c("CONTROL" = "grey", "ILD" = "red4")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  Group = sample_colors,
  Regulation = gene_colors
)

# Plot the heatmap
heatmap <- pheatmap(heatmap_data,
                    scale = "row",
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method = "complete",
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    fontsize = 8,
                    main = "Top 50 DEGs Heatmap: ILD vs Control")

# Save the heatmap
ggsave("ILDvsControl_Top50_DEGs_Heatmap.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)




### 4) Heatmap for IPF vs Control (Biopsy)
# Normalize counts using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE) 
vsd_mat <- assay(vsd)

# Select top 50 significant DEGs (FDR < 0.05 and |log2FC| > 1)
top50_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# Subset normalized expression matrix
heatmap_data <- vsd_mat[top50_genes$Gene_symbol, ]

# Sample annotation
sample_ann <- data.frame(Group = sample_info_4$condition)
rownames(sample_ann) <- rownames(sample_info_4)

# Gene annotation
gene_annotation <- data.frame(Regulation = top50_genes$regulation)
rownames(gene_annotation) <- top50_genes$Gene_symbol

# Color palettes
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)
sample_colors <- c("CONTROL" = "grey", "IPF" = "orangered")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  Group = sample_colors,
  Regulation = gene_colors
)

# Plot the heatmap
heatmap <- pheatmap(heatmap_data,
                    scale = "row",
                    show_rownames = TRUE,
                    show_colnames = TRUE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method = "complete",
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    fontsize = 8,
                    main = "Top 50 DEGs Heatmap: IPF vs Control (Biopsy)")

# Save the heatmap
ggsave("IPFvsControl_Biopsy_Top50_DEGs_Heatmap.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)





### 5) Heatmap for CHP vs Control (Biopsy)
# Normalize counts using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE) 
vsd_mat <- assay(vsd)

# Select top 50 significant DEGs (FDR < 0.05 and |log2FC| > 1)
top50_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# Subset normalized expression matrix
heatmap_data <- vsd_mat[top50_genes$Gene_symbol, ]

# Sample annotation
sample_ann <- data.frame(Group = sample_info_5$condition)
rownames(sample_ann) <- rownames(sample_info_5)

# Gene annotation
gene_annotation <- data.frame(Regulation = top50_genes$regulation)
rownames(gene_annotation) <- top50_genes$Gene_symbol

# Color palettes
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)
sample_colors <- c("CONTROL" = "grey", "CHP" = "cadetblue")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  Group = sample_colors,
  Regulation = gene_colors
)

# Plot the heatmap
heatmap <- pheatmap(heatmap_data,
                    scale = "row",
                    show_rownames = TRUE,
                    show_colnames = TRUE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method = "complete",
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    fontsize = 8,
                    main = "Top 50 DEGs Heatmap: CHP vs Control (Biopsy)")

# Save the heatmap
ggsave("CHPvsControl_Biopsy_Top50_DEGs_Heatmap.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)



### 6) Heatmap for CHP vs IPF (Biopsy)
# Normalize counts using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE) 
vsd_mat <- assay(vsd)

# Select top 50 significant DEGs (FDR < 0.05 and |log2FC| > 1)
top50_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# Subset normalized expression matrix
heatmap_data <- vsd_mat[top50_genes$Gene_symbol, ]

# Sample annotation
sample_ann <- data.frame(Group = sample_info_6$condition)
rownames(sample_ann) <- rownames(sample_info_6)

# Gene annotation
gene_annotation <- data.frame(Regulation = top50_genes$regulation)
rownames(gene_annotation) <- top50_genes$Gene_symbol

# Color palettes
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)
sample_colors <- c("CHP" = "cadetblue", "IPF" = "orangered")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  Group = sample_colors,
  Regulation = gene_colors
)

# Plot the heatmap
heatmap <- pheatmap(heatmap_data,
                    scale = "row",
                    show_rownames = TRUE,
                    show_colnames = TRUE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method = "complete",
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    fontsize = 8,
                    main = "Top 50 DEGs Heatmap: CHP vs IPF (Biopsy)")

# Save the heatmap
ggsave("CHPvsIPF_Biopsy_Top50_DEGs_Heatmap.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)





### 7) Heatmap for ILD vs Control (Biopsy)
# Normalize counts using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE) 
vsd_mat <- assay(vsd)

# Select top 50 significant DEGs (FDR < 0.05 and |log2FC| > 1)
top50_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# Subset normalized expression matrix
heatmap_data <- vsd_mat[top50_genes$Gene_symbol, ]

# Sample annotation
sample_ann <- data.frame(Group = sample_info.ild$condition)
rownames(sample_ann) <- rownames(sample_info.ild)

# Gene annotation
gene_annotation <- data.frame(Regulation = top50_genes$regulation)
rownames(gene_annotation) <- top50_genes$Gene_symbol

# Color palettes
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)
sample_colors <- c("CONTROL" = "grey", "ILD" = "red4")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  Group = sample_colors,
  Regulation = gene_colors
)

# Plot the heatmap
heatmap <- pheatmap(heatmap_data,
                    scale = "row",
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method = "complete",
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    fontsize = 8,
                    main = "Top 50 DEGs Heatmap: ILD vs Control (Biopsy)")

# Save the heatmap
ggsave("ILDvsControl_Biopsy_Top50_DEGs_Heatmap.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)




### 8) Heatmap for IPF vs Control (Explant)
# Normalize counts using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE) 
vsd_mat <- assay(vsd)

# Select top 50 significant DEGs (FDR < 0.05 and |log2FC| > 1)
top50_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# Subset normalized expression matrix
heatmap_data <- vsd_mat[top50_genes$Gene_symbol, ]

# Sample annotation
sample_ann <- data.frame(Group = sample_info_7$condition)
rownames(sample_ann) <- rownames(sample_info_7)

# Gene annotation
gene_annotation <- data.frame(Regulation = top50_genes$regulation)
rownames(gene_annotation) <- top50_genes$Gene_symbol

# Color palettes
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)
sample_colors <- c("CONTROL" = "grey", "IPF" = "orangered")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  Group = sample_colors,
  Regulation = gene_colors
)

# Plot the heatmap
heatmap <- pheatmap(heatmap_data,
                    scale = "row",
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method = "complete",
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    fontsize = 8,
                    main = "Top 50 DEGs Heatmap: IPF vs Control (Explant)")

# Save the heatmap
ggsave("IPFvsControl_Explant_Top50_DEGs_Heatmap.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)



### 9) Heatmap for CHP vs Control (Explant)
# Normalize counts using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE) 
vsd_mat <- assay(vsd)

# Select top 50 significant DEGs (FDR < 0.05 and |log2FC| > 1)
top50_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# Subset normalized expression matrix
heatmap_data <- vsd_mat[top50_genes$Gene_symbol, ]

# Sample annotation
sample_ann <- data.frame(Group = sample_info_8$condition)
rownames(sample_ann) <- rownames(sample_info_8)

# Gene annotation
gene_annotation <- data.frame(Regulation = top50_genes$regulation)
rownames(gene_annotation) <- top50_genes$Gene_symbol

# Color palettes
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)
sample_colors <- c("CONTROL" = "grey", "CHP" = "cadetblue")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  Group = sample_colors,
  Regulation = gene_colors
)

# Plot the heatmap
heatmap <- pheatmap(heatmap_data,
                    scale = "row",
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method = "complete",
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    fontsize = 8,
                    main = "Top 50 DEGs Heatmap: CHP vs Control (Explant)")

# Save the heatmap
ggsave("CHPvsControl_Explant_Top50_DEGs_Heatmap.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)



### 10) Heatmap for CHP vs IPF (Explant)
# Normalize counts using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE) 
vsd_mat <- assay(vsd)

# Select top 50 significant DEGs (FDR < 0.05 and |log2FC| > 1)
top50_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# Subset normalized expression matrix
heatmap_data <- vsd_mat[top50_genes$Gene_symbol, ]

# Sample annotation
sample_ann <- data.frame(Group = sample_info_9$condition)
rownames(sample_ann) <- rownames(sample_info_9)

# Gene annotation
gene_annotation <- data.frame(Regulation = top50_genes$regulation)
rownames(gene_annotation) <- top50_genes$Gene_symbol

# Color palettes
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)
sample_colors <- c("IPF" = "orangered", "CHP" = "cadetblue")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  Group = sample_colors,
  Regulation = gene_colors
)

# Plot the heatmap
heatmap <- pheatmap(heatmap_data,
                    scale = "row",
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method = "complete",
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    fontsize = 8,
                    main = "Top 50 DEGs Heatmap: CHP vs IPF (Explant)")

# Save the heatmap
ggsave("CHPvsIPF_Explant_Top50_DEGs_Heatmap.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)




### 11) Heatmap for ILD vs Control (Explant)
# Normalize counts using variance stabilizing transformation
vsd <- vst(dds, blind = FALSE) 
vsd_mat <- assay(vsd)

# Select top 50 significant DEGs (FDR < 0.05 and |log2FC| > 1)
top50_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(50)

# Subset normalized expression matrix
heatmap_data <- vsd_mat[top50_genes$Gene_symbol, ]

# Sample annotation
sample_ann <- data.frame(Group = sample_info.ild$condition)
rownames(sample_ann) <- rownames(sample_info.ild)

# Gene annotation
gene_annotation <- data.frame(Regulation = top50_genes$regulation)
rownames(gene_annotation) <- top50_genes$Gene_symbol

# Color palettes
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)
sample_colors <- c("CONTROL" = "grey", "ILD" = "red4")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  Group = sample_colors,
  Regulation = gene_colors
)

# Plot the heatmap
heatmap <- pheatmap(heatmap_data,
                    scale = "row",
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    clustering_distance_rows = "euclidean",
                    clustering_distance_cols = "euclidean",
                    clustering_method = "complete",
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    fontsize = 8,
                    main = "Top 50 DEGs Heatmap: ILD vs Control (Explant)")

# Save the heatmap
ggsave("ILDvsControl_Explant_Top50_DEGs_Heatmap.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)
