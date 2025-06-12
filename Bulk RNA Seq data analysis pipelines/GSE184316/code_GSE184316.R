# Set the working directory
setwd("C:\\IIT KGP\\GSE184316")
install.packages("BiocManager")

BiocManager::install("DESeq2", force = TRUE)
BiocManager::install("tidyverse")
BiocManager::install("GEOquery")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")


# Load required libraries
library(DESeq2)
library(tidyverse)
library(GEOquery)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(HTSFilter)

# Step 1: Import and prepare the data
# Replace 'counts.csv' with your actual file name
counts <- read.csv("GSE184316_raw_counts.csv", row.names = 1)

# Add gene symbols to results
gse <- getGEO("GSE184316", GSEMatrix = TRUE)
pheno.data <- gse$GSE184316_series_matrix.txt.gz@phenoData@data
pheno.data$`disease:ch1` <- ifelse(pheno.data$`disease:ch1` == "DONOR", "Control", pheno.data$`disease:ch1`)
colnames(pheno.data)[colnames(pheno.data) == "disease:ch1"] <- "Group"

# Extract Gene IDs from index
gene_ids <- rownames(counts)  # Ensure gene IDs are character type
gene_ids <- as.character(gene_ids)

# Map Gene IDs to Gene Symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = gene_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

# Replace NA values with original Gene IDs (if mapping fails)
gene_symbols[is.na(gene_symbols)] <- gene_ids[is.na(gene_symbols)]

# Make Gene Symbols unique by appending numbers to duplicates
gene_symbols <- make.unique(gene_symbols)

# Assign Gene Symbols as new row names
rownames(counts) <- gene_symbols

write.csv(counts, "GSE184316_COUNTS_GENESYMBOLS.csv")


filter_rows_with_found_column <- function(data, columns_to_check, values_to_keep) {
    # Ensure the input columns are valid
    if (!all(columns_to_check %in% seq_len(ncol(data)))) {
      stop("Some column indices are out of range.")
    }
  
# Filter rows and add the 'Group' column
  filtered_data <- data[apply(data[, columns_to_check, drop = FALSE], 1, function(row) {
    any(row %in% values_to_keep)  # Check if any value matches
  }), ]
  
# Add the 'Group' column
  filtered_data$Group <- apply(filtered_data[, columns_to_check, drop = FALSE], 1, function(row) {
    paste(row[row %in% values_to_keep], collapse = ", ")
  })
  
  return(filtered_data)
}

#For the DGE analysis between IPF and Control
##Starts
# Define the column indices to check and the values to filter for
counts_1 <- counts
columns_to_check <- c(2, 48)  # Column indices to check
values_to_keep <- c("IPF", "Control")  # Values to look for

filtered_pheno1 <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno1$Platform <- filtered_pheno1$platform_id
filtered_pheno1$Sample <- rownames(filtered_pheno1)
rownames(filtered_pheno1) <- NULL
filtered_pheno1 <- filtered_pheno1[, c(48, 54, 55)]
filtered_pheno1$Group_Encoded <- ifelse(filtered_pheno1$Group == "IPF", 1, 0)
counts_1 <- counts_1[, colnames(counts_1) %in% filtered_pheno1$Sample]

# Create a sample information data frame
sample_info <- data.frame(
  sample = colnames(counts_1),
  condition = factor(c(rep("CONTROL", 24), rep("IPF", 40)))
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info) <- sample_info$sample

# Step 2: Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_1,
                              colData = sample_info,
                              design = ~ condition)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Step 4: Run DESeq2
dds <- DESeq(dds)

# Step 5: Extract results
# Change 'condition_treatment_vs_control' if your levels are named differently
res <- results(dds, name="condition_IPF_vs_CONTROL")


# Step 7: Order results by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Step 8: Summarize results
summary(res)

# Step 9: MA plot
plotMA(res, ylim = c(-5, 5))

# Step 10: Export results with regulation information
res_df <- as.data.frame(res_ordered) %>%
  rownames_to_column("Gene_Symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up-regulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "IPFvsCtrl_deseq2_results.csv", row.names = FALSE)

# Step 11: Extract significantly differentially expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "GSE184316_FDEG_results_IPFvsCtrl.csv", row.names = FALSE)

# Step 12: Calculate and print summary statistics
total_genes <- nrow(res_df)
sig_up <- sum(res_df$regulation == "Up-regulated")
sig_down <- sum(res_df$regulation == "Down-regulated")
sig_total <- sig_up + sig_down

cat("Summary of Differential Expression Analysis:\n")
cat("Total genes analyzed:", total_genes, "\n")
cat("Total significantly differentially expressed genes:", sig_total, "\n")
cat("Up-regulated genes:", sig_up, "(", round(sig_up/total_genes*100, 2), "%)\n")
cat("Down-regulated genes:", sig_down, "(", round(sig_down/total_genes*100, 2), "%)\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Up-regulated", "Down-regulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c(round(sig_up/total_genes*100,4), round(sig_down/total_genes*100,4), round((total_genes - sig_total)/total_genes*100,4))
)

# Round the percentage to two decimal places
summary_df$Percentage <- round(summary_df$Percentage, 2)

# Print the summary dataframe
print(summary_df)

# Optionally, save the summary to a CSV file
write.csv(summary_df, "regulation_summary_IPFvsCtrl.csv", row.names = FALSE)

# Step 13: Create a bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")), 
            position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(title = "Summary of Gene Regulation in IPF vs Control",
       x = "Regulation Category",
       y = "Number of Genes") +
  scale_fill_manual(values = c("Up-regulated" = "red3", "Down-regulated" = "blue3", "Not significant" = "gray"))

ggsave("regulation_summary_plot_IPFvsCtrl.png", width = 10, height = 6)

# Step 14: Create a bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +  # Exclude "Not significant" in the plot
  geom_text(data = subset(summary_df, Category != "Not significant"),
            aes(label = paste0(Count, " (", Percentage, "%)")), 
            position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(title = "Summary of Gene Regulation in IPF vs Control",
       x = "Regulation Category",
       y = "Number of Genes") +
  scale_fill_manual(values = c("Up-regulated" = "red3", "Down-regulated" = "blue3"))
ggsave("regulation_summary_plot_IPFvsCtrl_2.png", width = 10, height = 6)

# Step 15: Create a customized volcano plot
# Set significance thresholds
logFC_threshold <- 1
pvalue_threshold <- 0.05

# Create the expression classification column
res_df$expression <- "Not Significant"
res_df$expression[res_df$log2FoldChange > logFC_threshold & res_df$padj < pvalue_threshold] <- "Upregulated"
res_df$expression[res_df$log2FoldChange < -logFC_threshold & res_df$padj < pvalue_threshold] <- "Downregulated"

# Filter and select top 10 upregulated genes
top_up <- res_df[res_df$expression == "Upregulated" & res_df$padj < pvalue_threshold, ]
top_up <- top_up[order(top_up$padj), ][1:15, ]

# Filter and select top 10 downregulated genes
top_down <- res_df[res_df$expression == "Downregulated" & res_df$padj < pvalue_threshold, ]
top_down <- top_down[order(top_down$padj), ][1:15, ]

# Combine top genes
top_genes <- rbind(top_up, top_down)

# Create the volcano plot
vplot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = expression)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(name = "Regulation",
                     values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not Significant" = "gray")) +
  theme_bw() +
  geom_text_repel(data = top_genes,
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_Symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_IPF vs Control",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_IPFvsCtrl.png", vplot, width = 12, height = 10)




#For the DGE analysis between HP and Control
##Starts
# Define the column indices to check and the values to filter for
counts_2 <- counts
columns_to_check <- c(2, 48)  # Column indices to check
values_to_keep <- c("HP", "Control")  # Values to look for

filtered_pheno2 <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno2$Platform <- filtered_pheno2$platform_id
filtered_pheno2$Sample <- rownames(filtered_pheno2)
rownames(filtered_pheno2) <- NULL
filtered_pheno2 <- filtered_pheno2[, c(48, 54, 55)]
filtered_pheno2$Group_Encoded <- ifelse(filtered_pheno2$Group == "HP", 1, 0)
counts_2 <- counts_2[, colnames(counts_2) %in% filtered_pheno2$Sample]

# Create a sample information data frame
# Replace with your actual sample information
sample_info <- data.frame(
  sample = colnames(counts_2),
  condition = factor(c(rep("CONTROL", 24), rep("HP", 36)))
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info) <- sample_info$sample

# Step 2: Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_2,
                              colData = sample_info,
                              design = ~ condition)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Step 4: Run DESeq2
dds <- DESeq(dds)

# Step 5: Extract results
# Change 'condition_treatment_vs_control' if your levels are named differently
res <- results(dds, name="condition_HP_vs_CONTROL")


# Step 7: Order results by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Step 8: Summarize results
summary(res)

# Step 9: MA plot
plotMA(res, ylim = c(-5, 5))

# Step 10: Export results with regulation information
res_df <- as.data.frame(res_ordered) %>%
  rownames_to_column("Gene_Symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up-regulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "HPvsCtrl_deseq2_results.csv", row.names = FALSE)

# Step 11: Extract significantly differentially expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "GSE184316_FDEG_results_HPvsCtrl.csv", row.names = FALSE)

# Step 12: Calculate and print summary statistics
total_genes <- nrow(res_df)
sig_up <- sum(res_df$regulation == "Up-regulated")
sig_down <- sum(res_df$regulation == "Down-regulated")
sig_total <- sig_up + sig_down

cat("Summary of Differential Expression Analysis:\n")
cat("Total genes analyzed:", total_genes, "\n")
cat("Total significantly differentially expressed genes:", sig_total, "\n")
cat("Up-regulated genes:", sig_up, "(", round(sig_up/total_genes*100, 2), "%)\n")
cat("Down-regulated genes:", sig_down, "(", round(sig_down/total_genes*100, 2), "%)\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Up-regulated", "Down-regulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c(sig_up/total_genes*100, sig_down/total_genes*100, (total_genes - sig_total)/total_genes*100)
)

# Round the percentage to two decimal places
summary_df$Percentage <- round(summary_df$Percentage, 2)

# Print the summary dataframe
print(summary_df)

# Optionally, save the summary to a CSV file
write.csv(summary_df, "regulation_summary_HPvsCtrl.csv", row.names = FALSE)

# Step 13: Create a bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")), 
            position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(title = "Summary of Gene Regulation in HP vs Control",
       x = "Regulation Category",
       y = "Number of Genes") +
  scale_fill_manual(values = c("Up-regulated" = "red3", "Down-regulated" = "blue3", "Not significant" = "gray"))

ggsave("regulation_summary_plot_HPvsCtrl.png", width = 10, height = 6)

# Step 14: Create a bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +  # Exclude "Not significant" in the plot
  geom_text(data = subset(summary_df, Category != "Not significant"),  
            aes(label = paste0(Count, " (", Percentage, "%)")), 
            position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(title = "Summary of Gene Regulation in HP vs Control",
       x = "Regulation Category",
       y = "Number of Genes") +
  scale_fill_manual(values = c("Up-regulated" = "red3", "Down-regulated" = "blue3"))
ggsave("regulation_summary_plot_HPvsCtrl_2.png", width = 10, height = 6)

# Step 15: Create a customized volcano plot
# Set significance thresholds
logFC_threshold <- 1
pvalue_threshold <- 0.05

# Create the expression classification column
res_df$expression <- "Not Significant"
res_df$expression[res_df$log2FoldChange > logFC_threshold & res_df$padj < pvalue_threshold] <- "Upregulated"
res_df$expression[res_df$log2FoldChange < -logFC_threshold & res_df$padj < pvalue_threshold] <- "Downregulated"

# Filter and select top 10 upregulated genes
top_up <- res_df[res_df$expression == "Upregulated" & res_df$padj < pvalue_threshold, ]
top_up <- top_up[order(top_up$padj), ][1:15, ]

# Filter and select top 10 downregulated genes
top_down <- res_df[res_df$expression == "Downregulated" & res_df$padj < pvalue_threshold, ]
top_down <- top_down[order(top_down$padj), ][1:15, ]

# Combine top genes
top_genes <- rbind(top_up, top_down)

# Create the volcano plot
vplot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = expression)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(name = "Regulation",
    values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not Significant" = "gray")) +
  theme_bw() +
  geom_text_repel(data = top_genes,
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_Symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_HP vs Control",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_HPvsCtrl.png", vplot, width = 12, height = 10)
##End




#For the DGE analysis between HP and IPF
##Starts
# Define the column indices to check and the values to filter for
counts_3 <- counts
columns_to_check <- c(2, 48)  # Column indices to check
values_to_keep <- c("HP", "IPF")  # Values to look for

filtered_pheno3 <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno3$Platform <- filtered_pheno3$platform_id
filtered_pheno3$Sample <- rownames(filtered_pheno3)
rownames(filtered_pheno3) <- NULL
filtered_pheno3 <- filtered_pheno3[, c(48, 54, 55)]
filtered_pheno3$Group_Encoded <- ifelse(filtered_pheno3$Group == "IPF", 0, 1)
counts_3 <- counts_3[, colnames(counts_3) %in% filtered_pheno3$Sample]


# Create a sample information data frame
# Replace with your actual sample information
sample_info <- data.frame(
  sample = colnames(counts_3),
  condition = factor(c(rep("HP", 36), rep("IPF", 40)), levels = c("IPF", "HP"))
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info) <- sample_info$sample

# Step 2: Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_3,
                              colData = sample_info,
                              design = ~ condition)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]




### Sample QC

vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)

sample_dist <- dist(t(vst_mat), method = "euclidean")

hc <- hclust(sample_dist, method = "complete")

# label sample names with condition
sample_labels <- paste(colnames(vst_mat), vsd$condition, sep = "_")

# Plot dendrogram with labels
plot(hc,
     labels = sample_labels,
     main = "Sample Clustering by Hierarchical Clustering (VST)",
     xlab = "", sub = "",
     hang = -1,
     cex = 0.8)  # shrink label size if long

# PCA
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pcaData$SampleName <- colnames(vsd)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  stat_ellipse(type = "norm", 
               level = 0.95, 
               aes(fill = condition), 
               alpha = 0.2, geom = "polygon", color = NA) +
  geom_text_repel(data = subset(pcaData, condition == "HP"),
            aes(x = PC1, y = PC2, label = SampleName),
            inherit.aes = FALSE,
            size = 4,
            max.overlaps = Inf,
            box.padding = 0.3,
            point.padding = 0.7,
            segment.color = NA) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: Before outlier removal (HP samples labelled)") +
  theme_minimal()


ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  stat_ellipse(type = "norm", 
               level = 0.95, 
               aes(fill = condition), 
               alpha = 0.2, geom = "polygon", color = NA) +
  geom_text_repel(data = subset(pcaData, condition == "IPF"),
                  aes(x = PC1, y = PC2, label = SampleName),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.5,
                  segment.color = NA) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: Before outlier removal (IPF samples labelled)") +
  theme_minimal()

# Sample distance heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- vsd$condition

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample Distance Heatmap")

## Sample correlation heatmap
cors <- cor(assay(vsd))  # Pearson correlation

# Sample annotation
annotation <- data.frame(Condition = colData(vsd)$condition)
rownames(annotation) <- colnames(vsd)

# Create color mapping for annotation
ann_colors <- list(Condition = c("HP" = "skyblue", "IPF" = "salmon"))

# 5. Plot
pheatmap(cors,
         annotation_col = annotation,
         annotation_row = annotation,
         annotation_colors = ann_colors,
         show_colnames = TRUE,
         show_rownames = TRUE,
         main = "Sample Correlation Heatmap",
         fontsize = 6,
         border_color = NA)

# boxplot
library(edgeR)
summary(colSums(counts_3))       # Library sizes
boxplot(log10(counts_3 + 1))     # Distribution across samples





# Step 4: Run DESeq2
dds <- DESeq(dds)

# Filter weakly expressed genes
filter <- HTSFilter(dds, s.len = 100, plot = FALSE)$filteredData
class(filter) #stays the same class - DESeqDataSet

dim(filter)

# Step 5: Extract results
# Change 'condition_treatment_vs_control' if your levels are named differently
res <- results(filter, independentFiltering = FALSE, name="condition_HP_vs_IPF")
head(res)
hist(res$pvalue)


# Step 7: Order results by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Step 8: Summarize results
summary(res)

# Step 9: MA plot
plotMA(res, ylim = c(-5, 5))

# Step 10: Export results with regulation information
res_df <- as.data.frame(res_ordered) %>%
  rownames_to_column("Gene_Symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up-regulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "HPvsIPF_deseq2_results.csv", row.names = FALSE)

# Step 11: Extract significantly differentially expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "GSE184316_FDEG_results_HPvsIPF.csv", row.names = FALSE)

# Step 12: Calculate and print summary statistics
total_genes <- nrow(res_df)
sig_up <- sum(res_df$regulation == "Up-regulated")
sig_down <- sum(res_df$regulation == "Down-regulated")
sig_total <- sig_up + sig_down

cat("Summary of Differential Expression Analysis:\n")
cat("Total genes analyzed:", total_genes, "\n")
cat("Total significantly differentially expressed genes:", sig_total, "\n")
cat("Up-regulated genes:", sig_up, "(", round(sig_up/total_genes*100, 2), "%)\n")
cat("Down-regulated genes:", sig_down, "(", round(sig_down/total_genes*100, 2), "%)\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Up-regulated", "Down-regulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c(sig_up/total_genes*100, sig_down/total_genes*100, (total_genes - sig_total)/total_genes*100)
)

# Round the percentage to two decimal places
summary_df$Percentage <- round(summary_df$Percentage, 2)

# Print the summary dataframe
print(summary_df)

# Optionally, save the summary to a CSV file
write.csv(summary_df, "regulation_summary_HPvsIPF.csv", row.names = FALSE)

# Step 13: Create a bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")), 
            position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(title = "Summary of Gene Regulation in HP vs IPF",
       x = "Regulation Category",
       y = "Number of Genes") +
  scale_fill_manual(values = c("Up-regulated" = "red3", "Down-regulated" = "blue3", "Not significant" = "gray"))

ggsave("regulation_summary_plot_HPvsIPF.png", width = 10, height = 6)

# Step 14: Create a bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +  # Exclude "Not significant" in the plot
  geom_text(data = subset(summary_df, Category != "Not significant"),  
            aes(label = paste0(Count, " (", Percentage, "%)")), 
            position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(title = "Summary of Gene Regulation in HP vs IPF",
       x = "Regulation Category",
       y = "Number of Genes") +
  scale_fill_manual(values = c("Up-regulated" = "red3", "Down-regulated" = "blue3"))
ggsave("regulation_summary_plot_HPvsIPF_2.png", width = 10, height = 6)

# Step 15: Create a customized volcano plot
# Set significance thresholds
logFC_threshold <- 1
pvalue_threshold <- 0.05

# Create the expression classification column
res_df$expression <- "Not Significant"
res_df$expression[res_df$log2FoldChange > logFC_threshold & res_df$padj < pvalue_threshold] <- "Upregulated"
res_df$expression[res_df$log2FoldChange < -logFC_threshold & res_df$padj < pvalue_threshold] <- "Downregulated"

# Filter and select top 10 upregulated genes
top_up <- res_df[res_df$expression == "Upregulated" & res_df$padj < pvalue_threshold, ]
top_up <- top_up[order(top_up$padj), ][1:15, ]

# Filter and select top 10 downregulated genes
top_down <- res_df[res_df$expression == "Downregulated" & res_df$padj < pvalue_threshold, ]
top_down <- top_down[order(top_down$padj), ][1:15, ]

# Combine top genes
top_genes <- rbind(top_up, top_down)

# Create the volcano plot
vplot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = expression)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(name = "Regulation",
                     values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not Significant" = "gray")) +
  theme_bw() +
  geom_text_repel(data = top_genes,
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_Symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_HP vs IPF",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_HPvsIPF.png", vplot, width = 12, height = 10)
##End





#For the DGE analysis between ILD and Control
##Starts
# Define the column indices to check and the values to filter for
counts_4 <- counts
columns_to_check <- c(2, 48)  # Column indices to check
values_to_keep <- c("HP", "IPF", "Control")  # Values to look for

filtered_pheno4 <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno4$Platform <- filtered_pheno4$platform_id
filtered_pheno4$Sample <- rownames(filtered_pheno4)
rownames(filtered_pheno4) <- NULL
filtered_pheno4 <- filtered_pheno4[, c(48, 54, 55)]

filtered_pheno4 <- filtered_pheno4 %>%
  mutate(Group = case_when(
    Group == "IPF" ~ "ILD",
    Group == "HP" ~ "ILD",
    TRUE ~ "Control"
  ))

filtered_pheno4$Group_Encoded <- ifelse(filtered_pheno4$Group == "ILD", 1, 0)
counts_4 <- counts_4[, colnames(counts_4) %in% filtered_pheno4$Sample]

# Create a sample information data frame
# Replace with your actual sample information
sample_info <- data.frame(
  sample = colnames(counts_4),
  condition = factor(c(rep("Control", 24), rep("ILD", 76)), levels = c("Control", "ILD"))
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info) <- sample_info$sample

# Step 2: Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_4,
                              colData = sample_info,
                              design = ~ condition)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Step 4: Run DESeq2
dds <- DESeq(dds)

# Step 5: Extract results
# Change 'condition_treatment_vs_control' if your levels are named differently
res <- results(dds, name="condition_ILD_vs_Control")


# Step 7: Order results by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Step 8: Summarize results
summary(res)

# Step 9: MA plot
plotMA(res, ylim = c(-5, 5))

# Step 10: Export results with regulation information
res_df <- as.data.frame(res_ordered) %>%
  rownames_to_column("Gene_Symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up-regulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "ILDvsControl_deseq2_results.csv", row.names = FALSE)

# Step 11: Extract significantly differentially expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "GSE184316_FDEG_results_ILDvsCtrl.csv", row.names = FALSE)

# Step 12: Calculate and print summary statistics
total_genes <- nrow(res_df)
sig_up <- sum(res_df$regulation == "Up-regulated")
sig_down <- sum(res_df$regulation == "Down-regulated")
sig_total <- sig_up + sig_down

cat("Summary of Differential Expression Analysis:\n")
cat("Total genes analyzed:", total_genes, "\n")
cat("Total significantly differentially expressed genes:", sig_total, "\n")
cat("Up-regulated genes:", sig_up, "(", round(sig_up/total_genes*100, 2), "%)\n")
cat("Down-regulated genes:", sig_down, "(", round(sig_down/total_genes*100, 2), "%)\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Up-regulated", "Down-regulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c(sig_up/total_genes*100, sig_down/total_genes*100, (total_genes - sig_total)/total_genes*100)
)

# Round the percentage to two decimal places
summary_df$Percentage <- round(summary_df$Percentage, 2)

# Print the summary dataframe
print(summary_df)

# Optionally, save the summary to a CSV file
write.csv(summary_df, "regulation_summary_ILDvsCtrl.csv", row.names = FALSE)

# Step 13: Create a bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")), 
            position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(title = "Summary of Gene Regulation in ILD vs Control",
       x = "Regulation Category",
       y = "Number of Genes") +
  scale_fill_manual(values = c("Up-regulated" = "red3", "Down-regulated" = "blue3", "Not significant" = "gray"))

ggsave("regulation_summary_plot_ILDvsCtrl.png", width = 10, height = 6)

# Step 14: Create a bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +  # Exclude "Not significant" in the plot
  geom_text(data = subset(summary_df, Category != "Not significant"),  
            aes(label = paste0(Count, " (", Percentage, "%)")), 
            position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(title = "Summary of Gene Regulation in ILD vs Control",
       x = "Regulation Category",
       y = "Number of Genes") +
  scale_fill_manual(values = c("Up-regulated" = "red3", "Down-regulated" = "blue3"))
ggsave("regulation_summary_plot_ILDvsCtrl_2.png", width = 10, height = 6)

# Create a customized volcano plot
# Set significance thresholds
logFC_threshold <- 1
pvalue_threshold <- 0.05

# Create the expression classification column
res_df$expression <- "Not Significant"
res_df$expression[res_df$log2FoldChange > logFC_threshold & res_df$padj < pvalue_threshold] <- "Upregulated"
res_df$expression[res_df$log2FoldChange < -logFC_threshold & res_df$padj < pvalue_threshold] <- "Downregulated"

# Filter and select top 10 upregulated genes
top_up <- res_df[res_df$expression == "Upregulated" & res_df$padj < pvalue_threshold, ]
top_up <- top_up[order(top_up$padj), ][1:15, ]

# Filter and select top 10 downregulated genes
top_down <- res_df[res_df$expression == "Downregulated" & res_df$padj < pvalue_threshold, ]
top_down <- top_down[order(top_down$padj), ][1:15, ]

# Combine top genes
top_genes <- rbind(top_up, top_down)

# Create the volcano plot
vplot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = expression)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(name = "Regulation",
                     values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not Significant" = "gray")) +
  theme_bw() +
  geom_text_repel(data = top_genes,
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_Symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_ILD vs Control",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_ILDvsCtrl.png", vplot, width = 12, height = 10)
##End




### HPvsIPF : QC and after QC analytics

# exlcude samples outside 95% CI zone in PCA
samples.to.be.excluded <- c("GSM5583981", "GSM5583973", "GSM5583990") #IPF samples

counts_QC <- counts_3[, !(colnames(counts_3) %in% samples.to.be.excluded)]

# Create a sample information data frame
# Replace with your actual sample information
sample_info <- data.frame(
  sample = colnames(counts_QC),
  condition = factor(c(rep("HP", 36), rep("IPF", 37)), levels = c("IPF", "HP"))
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info) <- sample_info$sample

# Step 2: Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts_QC,
                              colData = sample_info,
                              design = ~ condition)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Step 4: Run DESeq2
dds <- DESeq(dds)

# Step 5: Extract results
# Change 'condition_diseased_vs_control' if your levels are named differently
res <- results(dds, name="condition_HP_vs_IPF")


# Step 7: Order results by adjusted p-value
res_ordered <- res[order(res$padj), ]

# Step 8: Summarize results
summary(res)

# Step 9: MA plot
plotMA(res, ylim = c(-5, 5))

# Step 10: Export results with regulation information
res_df <- as.data.frame(res_ordered) %>%
  rownames_to_column("Gene_Symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up-regulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "HPvsIPF_deseq2_results_afterQC.csv", row.names = FALSE)

# Step 11: Extract significantly differentially expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "GSE184316_FDEG_results_HPvsIPF_afterQC.csv", row.names = FALSE)

# Step 12: Calculate and print summary statistics
total_genes <- nrow(res_df)
sig_up <- sum(res_df$regulation == "Up-regulated")
sig_down <- sum(res_df$regulation == "Down-regulated")
sig_total <- sig_up + sig_down

cat("Summary of Differential Expression Analysis:\n")
cat("Total genes analyzed:", total_genes, "\n")
cat("Total significantly differentially expressed genes:", sig_total, "\n")
cat("Up-regulated genes:", sig_up, "(", round(sig_up/total_genes*100, 2), "%)\n")
cat("Down-regulated genes:", sig_down, "(", round(sig_down/total_genes*100, 2), "%)\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Up-regulated", "Down-regulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c(sig_up/total_genes*100, sig_down/total_genes*100, (total_genes - sig_total)/total_genes*100)
)

# Round the percentage to two decimal places
summary_df$Percentage <- round(summary_df$Percentage, 2)

# Print the summary dataframe
print(summary_df)

# Optionally, save the summary to a CSV file
write.csv(summary_df, "regulation_summary_HPvsIPF_afterQC.csv", row.names = FALSE)

# Step 13: Create a bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")), 
            position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(title = "Summary of Gene Regulation in HP vs IPF_Post QC",
       x = "Regulation Category",
       y = "Number of Genes") +
  scale_fill_manual(values = c("Up-regulated" = "red3", "Down-regulated" = "blue3", "Not significant" = "gray"))

ggsave("regulation_summary_plot_HPvsIPF_afterQC.png", width = 10, height = 6)

# Step 14: Create a bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +  # Exclude "Not significant" in the plot
  geom_text(data = subset(summary_df, Category != "Not significant"),  
            aes(label = paste0(Count, " (", Percentage, "%)")), 
            position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(title = "Summary of Gene Regulation in HP vs IPF_Post QC",
       x = "Regulation Category",
       y = "Number of Genes") +
  scale_fill_manual(values = c("Up-regulated" = "red3", "Down-regulated" = "blue3"))
ggsave("regulation_summary_plot_HPvsIPF_2_afterQC.png", width = 10, height = 6)

# Step 15: Create a customized volcano plot
# Set significance thresholds
logFC_threshold <- 1
pvalue_threshold <- 0.05

# Create the expression classification column
res_df$expression <- "Not Significant"
res_df$expression[res_df$log2FoldChange > logFC_threshold & res_df$padj < pvalue_threshold] <- "Upregulated"
res_df$expression[res_df$log2FoldChange < -logFC_threshold & res_df$padj < pvalue_threshold] <- "Downregulated"

# Filter and select top 10 upregulated genes
top_up <- res_df[res_df$expression == "Upregulated" & res_df$padj < pvalue_threshold, ]
top_up <- top_up[order(top_up$padj), ][1:15, ]

# Filter and select top 10 downregulated genes
top_down <- res_df[res_df$expression == "Downregulated" & res_df$padj < pvalue_threshold, ]
top_down <- top_down[order(top_down$padj), ][1:15, ]

# Combine top genes
top_genes <- rbind(top_up, top_down)

# Create the volcano plot
vplot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = expression)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(name = "Regulation",
                     values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not Significant" = "gray")) +
  theme_bw() +
  geom_text_repel(data = top_genes,
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_Symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_HP vs IPF_Post QC",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_HPvsIPF_afterQC.png", vplot, width = 12, height = 10)
##End
