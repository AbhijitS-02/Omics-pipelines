## Heatmap - HP vs Control -------------------

# Select top 50 significant genes with |logFC| > 1
top50_genes <- res.HP.GPL14550 %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  arrange(adj.P.Val) %>%
  slice(1:50) %>%
  rownames_to_column(var = "Gene_Symbol")

# Annotation
filtered.merged.pheno.data <- sample.info.GPL14550 %>%
  filter(condition == "Control" | condition == "HP") %>%
  arrange(condition)

# Extract expression matrix for those genes
heatmap_matrix <- exprs.GPL14550[rownames(exprs.GPL14550) %in% top50_genes$Gene_Symbol, ]
heatmap_matrix <- heatmap_matrix[, colnames(heatmap_matrix) %in% filtered.merged.pheno.data$geo_accession]


# Set sample annotations
sample_ann <- filtered.merged.pheno.data %>%
  select(condition)

# Set gene annotation
gene_annotation <- top50_genes %>%
  select(Gene_Symbol, expression)

gene_annotation <- column_to_rownames(gene_annotation, var = "Gene_Symbol")

# Set colors
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)

sample_colors <- c("Control" = "grey", "HP" = "cadetblue")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  condition = sample_colors,
  expression = gene_colors
)
# Plot and save heatmap
heatmap <- pheatmap(heatmap_matrix,
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    clustering_distance_cols = "euclidean",
                    clustering_distance_rows = "euclidean",
                    clustering_method = "complete",
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    scale = "row",
                    fontsize = 8,
                    main = "Heatmap of DEGs_HP_vs_Control_GPL14550")
ggsave("HPvsControl_heatmap_GPL14550.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)
#END





## Heatmap - IPF vs Control -------------------

# Select top 50 significant genes with |logFC| > 1
top50_genes <- res.IPF.GPL14550 %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  arrange(adj.P.Val) %>%
  slice(1:50) %>%
  rownames_to_column(var = "Gene_Symbol")

# Annotation
filtered.merged.pheno.data <- sample.info.GPL14550 %>%
  filter(condition == "Control" | condition == "IPF") %>%
  arrange(condition)

# Extract expression matrix for those genes
heatmap_matrix <- exprs.GPL14550[rownames(exprs.GPL14550) %in% top50_genes$Gene_Symbol, ]
heatmap_matrix <- heatmap_matrix[, colnames(heatmap_matrix) %in% filtered.merged.pheno.data$geo_accession]


# Set sample annotations
sample_ann <- filtered.merged.pheno.data %>%
  select(condition)

# Set gene annotation
gene_annotation <- top50_genes %>%
  select(Gene_Symbol, expression)

gene_annotation <- column_to_rownames(gene_annotation, var = "Gene_Symbol")

# Set colors
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)

sample_colors <- c("Control" = "grey", "IPF" = "orangered")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  condition = sample_colors,
  expression = gene_colors
)
# Plot and save heatmap
heatmap <- pheatmap(heatmap_matrix,
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    clustering_distance_cols = "euclidean",
                    clustering_distance_rows = "euclidean",
                    clustering_method = "complete",
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    scale = "row",
                    fontsize = 8,
                    main = "Heatmap of top 50 DEGs_IPF_vs_Control_GPL14550")
ggsave("IPFvsControl_heatmap_GPL14550.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)
#END




## Heatmap - ILD vs Control -------------------

# Select top 50 significant genes with |logFC| > 1
top50_genes <- res.ild.GPL14550 %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  arrange(adj.P.Val) %>%
  slice(1:50) %>%
  rownames_to_column(var = "Gene_Symbol")

# Set sample annotations
sample_ann <- sample.info.GPL14550 %>%
  arrange(disease) %>%
  select(disease)

# Extract expression matrix for those genes
heatmap_matrix <- exprs.GPL14550[rownames(exprs.GPL14550) %in% top50_genes$Gene_Symbol, ]
heatmap_matrix <- heatmap_matrix[, colnames(heatmap_matrix) %in% rownames(sample_ann)]

# Set gene annotation
gene_annotation <- top50_genes %>%
  select(Gene_Symbol, expression)

gene_annotation <- column_to_rownames(gene_annotation, var = "Gene_Symbol")

# Set colors
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)

sample_colors <- c("Control" = "grey", "ILD" = "red4")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  disease = sample_colors,
  expression = gene_colors
)
# Plot and save heatmap
heatmap <- pheatmap(heatmap_matrix,
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    clustering_distance_cols = "euclidean",
                    clustering_distance_rows = "euclidean",
                    clustering_method = "complete",
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    scale = "row",
                    fontsize = 8,
                    main = "Heatmap of top 50 DEGs_ILD_vs_Control_GPL14550")
ggsave("ILDvsControl_heatmap_GPL14550.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)
#END



## Heatmap - HP vs IPF -------------------

# Select top 50 significant genes with |logFC| > 1
top50_genes <- res.hp.ipf %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  arrange(adj.P.Val) %>%
  slice(1:50) %>%
  rownames_to_column(var = "Gene_Symbol")

# Set sample annotations
sample_ann <- sample.info.hp.ipf %>%
  arrange(condition) %>%
  select(condition)

# Extract expression matrix for those genes
heatmap_matrix <- exprs.GPL14550[rownames(exprs.GPL14550) %in% top50_genes$Gene_Symbol, ]
heatmap_matrix <- heatmap_matrix[, colnames(heatmap_matrix) %in% rownames(sample_ann)]

# Set gene annotation
gene_annotation <- top50_genes %>%
  select(Gene_Symbol, expression)

gene_annotation <- column_to_rownames(gene_annotation, var = "Gene_Symbol")

# Set colors
heatmap_colors <- colorRampPalette(c("#2D004B", "#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026", 
                                     "#FF7F00"))(100)

sample_colors <- c("IPF" = "orangered", "HP" = "cadetblue")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  condition = sample_colors,
  expression = gene_colors
)
# Plot and save heatmap
heatmap <- pheatmap(heatmap_matrix,
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    clustering_distance_cols = "euclidean",
                    clustering_distance_rows = "euclidean",
                    clustering_method = "complete",
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    scale = "row",
                    fontsize = 8,
                    main = "Heatmap of DEGs_HP_vs_IPF_GPL14550")
ggsave("HPvsIPF_GPL14550_heatmap.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)
#END