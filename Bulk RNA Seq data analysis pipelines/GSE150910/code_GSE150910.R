#Set working directory
setwd("C:/IIT KGP/GSE150910")

# Load required libraries
library(DESeq2)
library(tidyverse)
library(GEOquery)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(HTSFilter)

# Import and prepare the data
gene_counts = read.csv("GSE150910_gene_counts.csv", row.names = 1) #change the data file accordingly

# Get phenotypic data
gse <- getGEO("GSE150910", GSEMatrix = TRUE)
pheno.data <- gse$GSE150910_series_matrix.txt.gz@phenoData@data

# Create proper Group column for effective classification
pheno.data <- pheno.data %>%
  mutate(characteristics_ch1.3 = gsub("diagnosis: ", "",characteristics_ch1.3))
colnames(pheno.data)[colnames(pheno.data) == "characteristics_ch1.3"] <- "Group"

# Change classification to CHP, IPF and CONTROL
pheno.data$Group <- ifelse(pheno.data$Group == "chp", "CHP", pheno.data$Group)
pheno.data$Group <- ifelse(pheno.data$Group == "ipf", "IPF", pheno.data$Group)
pheno.data$Group <- ifelse(pheno.data$Group == "control", "CONTROL", pheno.data$Group)

# No gene symbol mapping required as already present

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








# PCA plot for all samples
counts_copy <- gene_counts
filtered_pheno <- pheno.data[, c(1,13,62)]
rownames(filtered_pheno) <- NULL

#Group encoding
filtered_pheno$Group_Encoded <- ifelse(filtered_pheno$Group == "CHP", 2, 
                                       ifelse(filtered_pheno$Group == "IPF", 1, 0))
rownames(filtered_pheno) <- filtered_pheno$title
filtered_pheno <- filtered_pheno[match(colnames(counts_copy), rownames(filtered_pheno)), , drop = FALSE]

#Create sample
sample_info_bulk <- data.frame(
  sample <- colnames(counts_copy),
  condition <- factor(filtered_pheno$Group),
  plate <- factor(pheno.data$`plate:ch1`),
  tissue <- factor(filtered_pheno$`sample type:ch1`)
  )

# Ensure row names of sample_info match column names of counts
rownames(sample_info_bulk) <- sample_info_bulk$sample

colnames(sample_info_bulk)[colnames(sample_info_bulk) == "sample....colnames.counts_copy."] <- "sample"
colnames(sample_info_bulk)[colnames(sample_info_bulk) == "condition....factor.filtered_pheno.Group."] <- "condition"
colnames(sample_info_bulk)[colnames(sample_info_bulk) == "plate....factor.pheno.data..plate.ch1.."] <- "plate"
colnames(sample_info_bulk)[colnames(sample_info_bulk) == "tissue....factor.filtered_pheno..sample.type.ch1.."] <- "tissue"

# Make CONTROL the first level since R by default will take CHP as reference level
sample_info_bulk$condition <- factor(sample_info_bulk$condition, levels = c("CONTROL", "IPF", "CHP"))

#Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts_copy,
  colData = sample_info_bulk,
  design = ~ plate + tissue + condition
)

# Variance Stabilizing Transformation
vsd <- vst(dds, blind=TRUE)

# Basic PCA plot
plotPCA(vsd, intgroup= "condition")
plotPCA(vsd, intgroup= "plate")
plotPCA(vsd, intgroup= "tissue")

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Extract VST matrix
vst_data_mat <- assay(vsd)
vst_data_mat_filtered <- vst_data_mat[apply(vst_data_mat, 1, var) != 0, ]

# Perform PCA
pca_res <- prcomp(t(vst_data_mat_filtered), scale. = TRUE)

# Check how much variance each PC explains
summary(pca_res)

pca_df <- as.data.frame(pca_res$x)
pca_df$condition <- colData(dds)$condition

# Plot PC1 vs PC2 for Condition
ggplot(pca_df, aes(x = PC1, y = PC2, color = as.factor(condition))) +
  geom_point(size = 2) +
  labs(title = "PCA: PC1 vs PC2", color = "condition") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))






## Tissue PCA plot
pcaData <- plotPCA(vsd, intgroup = "tissue", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Extract VST matrix
vst_data_mat <- assay(vsd)
vst_data_mat_filtered <- vst_data_mat[apply(vst_data_mat, 1, var) != 0, ]

# Perform PCA
pca_res <- prcomp(t(vst_data_mat_filtered), scale. = TRUE)

# Check how much variance each PC explains
summary(pca_res)

pca_df <- as.data.frame(pca_res$x)
pca_df$tissue <- colData(dds)$tissue

# Plot PC1 vs PC2 for Tissue
tissue_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = as.factor(tissue))) +
  geom_point(size = 3) +
  labs(title = "PCA: PC1 vs PC2", color = "tissue") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("PCA plot_tissue type.png", tissue_pca, dpi = 300, width = 12, height = 9)










### Check for variation due to "Plate"
pcaData <- plotPCA(vsd, intgroup = "plate", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Extract VST matrix
vst_data_mat <- assay(vsd)
vst_data_mat_filtered <- vst_data_mat[apply(vst_data_mat, 1, var) != 0, ]

# Perform PCA
pca_res <- prcomp(t(vst_data_mat_filtered), scale. = TRUE)

# Check how much variance each PC explains
summary(pca_res)

pca_df <- as.data.frame(pca_res$x)
pca_df$plate <- colData(dds)$plate

# Plot PC1 vs PC2 for Plate
ggplot(pca_df, aes(x = PC1, y = PC2, color = as.factor(plate))) +
  geom_point(size = 2) +
  labs(title = "PCA: PC1 vs PC2", color = "plate") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Plot PC1 vs PC3 for Plate
ggplot(pca_df, aes(x = PC1, y = PC3, color = as.factor(plate))) +
  geom_point(size = 2) +
  labs(title = "PCA: PC1 vs PC3", color = "plate") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Plot PC2 vs PC3 for Plate
ggplot(pca_df, aes(x = PC2, y = PC3, color = as.factor(plate))) +
  geom_point(size = 2) +
  labs(title = "PCA: PC2 vs PC3", color = "plate") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Plot PC3 vs PC4 for Plate
ggplot(pca_df, aes(x = PC3, y = PC4, color = as.factor(plate))) +
  geom_point(size = 2) +
  labs(title = "PCA: PC3 vs PC4", color = "plate") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Plot PC1 vs PC4 for Plate
ggplot(pca_df, aes(x = PC1, y = PC4, color = as.factor(plate))) +
  geom_point(size = 2) +
  labs(title = "PCA: PC1 vs PC4", color = "plate") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))











#For the DGE analysis between IPF and CONTROL
##Starts
# Define the column indices to check and the values to filter for
gene_counts_1 <- gene_counts
columns_to_check <- c(1:64) # Check all 64 columns
values_to_keep <- c("IPF", "CONTROL") # Values to look for

filtered_pheno1 <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno1$Platform <- filtered_pheno1$platform_id
rownames(filtered_pheno1) <- NULL
filtered_pheno1 <- filtered_pheno1[,c(1,13)]

# Group encoding
filtered_pheno1$Group_Encoded = ifelse(filtered_pheno1$Group == "CONTROL", 0, 1)
# Since gene_counts_1 has title as column names, filtered_pheno1$title is used to filter the date for gene_counts_1 
gene_counts_1 <- gene_counts_1[, colnames(gene_counts_1) %in% filtered_pheno1$title] 

# Matching order of column names of counts with row names of filtered pheno data
rownames(filtered_pheno1) <- filtered_pheno1$title
filtered_pheno1 <- filtered_pheno1[match(colnames(gene_counts_1), rownames(filtered_pheno1)), , drop = FALSE]

# Create a sample information data frame
sample_info <- data.frame(
  sample = colnames(gene_counts_1),
  condition = factor(filtered_pheno1$Group)
  )

# Ensure row names of sample_info match column names of counts
rownames(sample_info) <- sample_info$sample

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = gene_counts_1,
                              colData = sample_info,
                              design = ~ condition)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, name= "condition_IPF_vs_CONTROL")

# Order according to adjusted p values
res_ordered <- res[order(res$padj), ]

# Summarize results
summary(res)

# MA plot
plotMA(res, ylim = c(-5,5))

#Export results with regulation information
res_df <- as.data.frame(res) %>%
  rownames_to_column("Gene_symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "IPFvsCtrl_deseq2_results.csv", row.names = FALSE)

#Extract significantly expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "IPFvsCtrl_sig_exp_genes_results.csv", row.names = FALSE)

# Summary statistics
total_genes = nrow(res_df)
sig_up = sum(res_df$regulation == "Upregulated")
sig_down = sum(res_df$regulation == "Downregulated")
sig_total = sig_up + sig_down

# Print summary statistics
cat("Summary of Differential Expression Analysis :\n ")
cat("Total number of genes analyzed = ",total_genes, "\n")
cat("Number of significantly upregulated genes = ",sig_up,"\n")
cat("Number of significantly downregulated genes = ",sig_down,"\n")
cat("Total number of significantly differentially expressed genes = ",sig_total,"\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Upregulated", "Downregulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c((sig_up/total_genes)*100, (sig_down/total_genes)*100, ((total_genes - sig_total)/total_genes)*100)
)

# Round up to 2 decimal places
summary_df$Percentage <- round(summary_df$Percentage,2)

# Print summary dataframe
print(summary_df)

# Export (optional)
write.csv(summary_df, "Summary of Regulation_IPFvsControl.csv", row.names = FALSE)

# Create bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Gene regulation summary in IPF vs Control",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not significant" = "grey" ))

ggsave("Regulation summary plot_IPFvsCtrl.png", width = 10, height = 10)

# Create bar plot of regulation summary but excluding non signicant genes
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +
  geom_text(data = subset(summary_df, Category != "Not significant"),
            aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Summary of gene regulation in IPF vs Control",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3"))

ggsave("Regulation_summary_IPFvsCtrl_2.png", height = 10, width = 10)

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
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_symbol),
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
##End


# For the DGE analysis between HP and Control
# Starts
gene_counts_2 <- gene_counts
columns_to_check <- c(1:64)
values_to_keep <- c("CHP", "CONTROL")

filtered_pheno2 <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno2$Platform <- filtered_pheno1$platform_id
rownames(filtered_pheno2) <- NULL
filtered_pheno2 <- filtered_pheno2[,c(1,13)]

# Group encoding
filtered_pheno2$Group_Encoded = ifelse(filtered_pheno2$Group == "CONTROL", 0, 1)
# Since gene_counts_2 has title as column names, filtered_pheno2$title is used to filter the date for gene_counts_2 
gene_counts_2 <- gene_counts_2[, colnames(gene_counts_2) %in% filtered_pheno2$title] 

# Matching order of column names of counts with row names of filtered pheno data
rownames(filtered_pheno2) <- filtered_pheno2$title
filtered_pheno2 <- filtered_pheno2[match(colnames(gene_counts_2), rownames(filtered_pheno2)), , drop = FALSE]

# Create Sample data frame
sample_info_2 <- data.frame(
  Sample = colnames(gene_counts_2),
  condition = factor(filtered_pheno2$Group)
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info_2) <- sample_info_2$Sample

# Create a DESeq2 object
dds = DESeqDataSetFromMatrix(
  countData = gene_counts_2,
  colData = sample_info_2,
  design = ~ condition
)
dds$condition <- relevel(dds$condition, ref = "CONTROL")

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, name= "condition_CHP_vs_CONTROL")

# Order according to adjusted p values
res_ordered <- res[order(res$padj), ]

# Summarize results
summary(res)

# MA plot
plotMA(res, ylim = c(-5,5))

#Export results with regulation information
res_df <- as.data.frame(res) %>%
  rownames_to_column("Gene_symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "CHPvsCtrl_deseq2_results.csv", row.names = FALSE)

#Extract significantly expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "CHPvsCtrl_sig_exp_genes_results.csv", row.names = FALSE)

# Summary statistics
total_genes = nrow(res_df)
sig_up = sum(res_df$regulation == "Upregulated")
sig_down = sum(res_df$regulation == "Downregulated")
sig_total = sig_up + sig_down

# Print summary statistics
cat("Summary of Differential Expression Analysis :\n ")
cat("Total number of genes analyzed = ",total_genes, "\n")
cat("Number of significantly upregulated genes = ",sig_up,"\n")
cat("Number of significantly downregulated genes = ",sig_down,"\n")
cat("Total number of significantly differentially expressed genes = ",sig_total,"\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Upregulated", "Downregulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c((sig_up/total_genes)*100, (sig_down/total_genes)*100, ((total_genes - sig_total)/total_genes)*100)
)

# Round up to 2 decimal places
summary_df$Percentage <- round(summary_df$Percentage,2)

# Print summary dataframe
print(summary_df)

# Export (optional)
write.csv(summary_df, "Summary of Regulation_CHPvsControl.csv", row.names = FALSE)

# Create bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Gene regulation summary in HP vs Control",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not significant" = "grey" ))

ggsave("Regulation summary plot_CHPvsCtrl.png", width = 10, height = 10)

# Create bar plot of regulation summary but excluding non signicant genes
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +
  geom_text(data = subset(summary_df, Category != "Not significant"),
            aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Summary of gene regulation in HP vs Control",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3"))

ggsave("Regulation_summary_CHPvsCtrl_2.png", height = 10, width = 10)

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
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_symbol),
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


# For the DGE analysis between CHP and IPF
# Starts
gene_counts_3 <- gene_counts
columns_to_check <- c(1:64)
values_to_keep <- c("CHP", "IPF")

filtered_pheno3 <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno3$Platform <- filtered_pheno3$platform_id
rownames(filtered_pheno3) <- NULL
filtered_pheno3 <- filtered_pheno3[,c(1,13)]

# Group encoding
filtered_pheno3$Group_Encoded = ifelse(filtered_pheno3$Group == "IPF", 0, 1)
# Since gene_counts_3 has title as column names, filtered_pheno3$title is used to filter the date for gene_counts_3 
gene_counts_3 <- gene_counts_3[, colnames(gene_counts_3) %in% filtered_pheno3$title] 

# Matching order of column names of counts with row names of filtered pheno data
rownames(filtered_pheno3) <- filtered_pheno3$title
filtered_pheno3 <- filtered_pheno3[match(colnames(gene_counts_3), rownames(filtered_pheno3)), , drop = FALSE]

# Create sample info
sample_info_3 <- data.frame(
  sample = colnames(gene_counts_3),
  condition = factor(filtered_pheno3$Group)
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info_3) <- sample_info_3$sample

# Make IPF as the reference level
sample_info_3$condition <- factor(sample_info_3$condition, levels = c("IPF", "CHP"))
# Create a DESeq2 object
dds = DESeqDataSetFromMatrix(
  countData = gene_counts_3,
  colData = sample_info_3,
  design = ~ condition
)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, name= "condition_CHP_vs_IPF")

# Order according to adjusted p values
res_ordered <- res[order(res$padj), ]

# Summarize results
summary(res)

# MA plot
plotMA(res, ylim = c(-5,5))

#Export results with regulation information
res_df <- as.data.frame(res) %>%
  rownames_to_column("Gene_symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "CHPvsIPF_deseq2_results.csv", row.names = FALSE)

#Extract significantly expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "CHPvsIPF_sig_exp_genes_results.csv", row.names = FALSE)

# Summary statistics
total_genes = nrow(res_df)
sig_up = sum(res_df$regulation == "Upregulated")
sig_down = sum(res_df$regulation == "Downregulated")
sig_total = sig_up + sig_down

# Print summary statistics
cat("Summary of Differential Expression Analysis :\n ")
cat("Total number of genes analyzed = ",total_genes, "\n")
cat("Number of significantly upregulated genes = ",sig_up,"\n")
cat("Number of significantly downregulated genes = ",sig_down,"\n")
cat("Total number of significantly differentially expressed genes = ",sig_total,"\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Upregulated", "Downregulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c((sig_up/total_genes)*100, (sig_down/total_genes)*100, ((total_genes - sig_total)/total_genes)*100)
)

# Round up to 2 decimal places
summary_df$Percentage <- round(summary_df$Percentage,2)

# Print summary dataframe
print(summary_df)

# Export (optional)
write.csv(summary_df, "Summary of Regulation_CHPvsIPF.csv", row.names = FALSE)

# Create bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Gene regulation summary in HP vs IPF",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not significant" = "grey" ))

ggsave("Regulation summary plot_CHPvsIPF.png", width = 10, height = 10)

# Create bar plot of regulation summary but excluding non signicant genes
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +
  geom_text(data = subset(summary_df, Category != "Not significant"),
            aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Summary of gene regulation in HP vs IPF",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3"))

ggsave("Regulation_summary_CHPvsIPF_2.png", height = 10, width = 10)

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
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_symbol),
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








# For the DGE analysis between ILD and Control
# Starts
gene_counts.ild <- gene_counts
columns_to_check <- c(1:64)
values_to_keep <- c("CHP", "IPF", "CONTROL")

filtered_pheno.ild <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno.ild$Platform <- filtered_pheno.ild$platform_id
rownames(filtered_pheno.ild) <- NULL
filtered_pheno.ild <- filtered_pheno.ild[,c(1,13)]

filtered_pheno.ild <- filtered_pheno.ild %>%
  mutate(Group = case_when(
    Group == "IPF" ~ "ILD",
    Group == "CHP" ~ "ILD",
    TRUE ~ "CONTROL"
  ))

# Group encoding
filtered_pheno.ild$Group_Encoded = ifelse(filtered_pheno.ild$Group == "ILD", 1, 0)
# Since gene_counts_3 has title as column names, filtered_pheno3$title is used to filter the date for gene_counts_3 
gene_counts.ild <- gene_counts.ild[, colnames(gene_counts.ild) %in% filtered_pheno.ild$title] 

# Matching order of column names of counts with row names of filtered pheno data
rownames(filtered_pheno.ild) <- filtered_pheno.ild$title
filtered_pheno.ild <- filtered_pheno.ild[match(colnames(gene_counts.ild), rownames(filtered_pheno.ild)), , 
                                         drop = FALSE]

# Create sample info
sample_info.ild <- data.frame(
  sample = colnames(gene_counts.ild),
  condition = factor(filtered_pheno.ild$Group)
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info.ild) <- sample_info.ild$sample

# Make IPF as the reference level
sample_info.ild$condition <- factor(sample_info.ild$condition, levels = c("CONTROL", "ILD"))
# Create a DESeq2 object
dds = DESeqDataSetFromMatrix(
  countData = gene_counts.ild,
  colData = sample_info.ild,
  design = ~ condition
)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, name= "condition_ILD_vs_CONTROL")

# Order according to adjusted p values
res_ordered <- res[order(res$padj), ]

# Summarize results
summary(res)

# MA plot
plotMA(res, ylim = c(-5,5))

#Export results with regulation information
res_df <- as.data.frame(res) %>%
  rownames_to_column("Gene_symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "ILDvsCONTROL_deseq2_results.csv", row.names = FALSE)

#Extract significantly expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "ILDvsCONTROL_sig_exp_genes_results.csv", row.names = FALSE)

# Summary statistics
total_genes = nrow(res_df)
sig_up = sum(res_df$regulation == "Upregulated")
sig_down = sum(res_df$regulation == "Downregulated")
sig_total = sig_up + sig_down

# Print summary statistics
cat("Summary of Differential Expression Analysis :\n ")
cat("Total number of genes analyzed = ",total_genes, "\n")
cat("Number of significantly upregulated genes = ",sig_up,"\n")
cat("Number of significantly downregulated genes = ",sig_down,"\n")
cat("Total number of significantly differentially expressed genes = ",sig_total,"\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Upregulated", "Downregulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c((sig_up/total_genes)*100, (sig_down/total_genes)*100, ((total_genes - sig_total)/total_genes)*100)
)

# Round up to 2 decimal places
summary_df$Percentage <- round(summary_df$Percentage,2)

# Print summary dataframe
print(summary_df)

# Export (optional)
write.csv(summary_df, "Summary of Regulation_ILDvsCONTROL.csv", row.names = FALSE)

# Create bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Gene regulation summary for ILD vs Control",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not significant" = "grey" ))

ggsave("Regulation summary plot_ILDvsCONTROL.png", width = 10, height = 10)

# Create bar plot of regulation summary but excluding non signicant genes
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +
  geom_text(data = subset(summary_df, Category != "Not significant"),
            aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Summary of gene regulation for ILD vs Control",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3"))

ggsave("Regulation_summary_ILDvsCONTROL_2.png", height = 10, width = 10)

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
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_symbol),
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



# BIOPSY SAMPLES
setwd("C:/IIT KGP/GSE150910/Biopsy samples")

# DGE analysis of IPF vs Control biopsy samples
##Start

# Create proper Sample type column for effective classification
pheno.data <- pheno.data %>%
  mutate(characteristics_ch1.7 = gsub("sample type: ", "",characteristics_ch1.7))
colnames(pheno.data)[colnames(pheno.data) == "characteristics_ch1.7"] <- "Sample type"

# Change classification to CHP, IPF and CONTROL
pheno.data$`Sample type` <- ifelse(pheno.data$`Sample type` == "explant", "Explant", pheno.data$`Sample type`)
pheno.data$`Sample type` <- ifelse(pheno.data$`Sample type` == "biopsy", "Biopsy", pheno.data$`Sample type`)

# Define the column indices to check and the values to filter for
gene_counts_4 <- gene_counts
columns_to_check <- c(1:64) # Check all 64 columns
values_to_keep <- c("IPF", "CONTROL") # Values to look for

filtered_pheno4 <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno4$Platform <- filtered_pheno4$platform_id
rownames(filtered_pheno4) <- NULL
filtered_pheno4 <- filtered_pheno4[,c(1,13,17)]

# Filter out the explant samples
filtered_pheno4 <- filtered_pheno4[filtered_pheno4$`Sample type` == "Biopsy",]

# Group encoding
filtered_pheno4$Group_Encoded = ifelse(filtered_pheno4$Group == "CONTROL", 0, 1)
# Since gene_counts_4 has title as column names, filtered_pheno4$title is used to filter the date for gene_counts_4 
gene_counts_4 <- gene_counts_4[, colnames(gene_counts_4) %in% filtered_pheno4$title] 

# Matching order of column names of counts with row names of filtered pheno data
rownames(filtered_pheno4) <- filtered_pheno4$title
filtered_pheno4 <- filtered_pheno4[match(colnames(gene_counts_4), rownames(filtered_pheno4)), , drop = FALSE]

# Create sample info
sample_info_4 <- data.frame(
  sample = colnames(gene_counts_4),
  condition = factor(filtered_pheno4$Group)
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info_4) <- sample_info_4$sample

# Create a DESeq2 object
dds = DESeqDataSetFromMatrix(
  countData = gene_counts_4,
  colData = sample_info_4,
  design = ~ condition
)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, name= "condition_IPF_vs_CONTROL")

# Order according to adjusted p values
res_ordered <- res[order(res$padj), ]

# Summarize results
summary(res)

# MA plot
plotMA(res, ylim = c(-5,5))

#Export results with regulation information
res_df <- as.data.frame(res) %>%
  rownames_to_column("Gene_symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "IPFvsCtrl(Biopsy)_deseq2_results.csv", row.names = FALSE)

#Extract significantly expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "IPFvsCtrl(Biopsy)_sig_exp_genes_results.csv", row.names = FALSE)

# Summary statistics
total_genes = nrow(res_df)
sig_up = sum(res_df$regulation == "Upregulated")
sig_down = sum(res_df$regulation == "Downregulated")
sig_total = sig_up + sig_down

# Print summary statistics
cat("Summary of Differential Expression Analysis :\n ")
cat("Total number of genes analyzed = ",total_genes, "\n")
cat("Number of significantly upregulated genes = ",sig_up,"\n")
cat("Number of significantly downregulated genes = ",sig_down,"\n")
cat("Total number of significantly differentially expressed genes = ",sig_total,"\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Upregulated", "Downregulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c((sig_up/total_genes)*100, (sig_down/total_genes)*100, ((total_genes - sig_total)/total_genes)*100)
)

# Round up to 2 decimal places
summary_df$Percentage <- round(summary_df$Percentage,2)

# Print summary dataframe
print(summary_df)

# Export (optional)
write.csv(summary_df, "Summary of Regulation_IPFvsCtrl_Biopsy.csv", row.names = FALSE)

# Create bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Gene regulation summary in IPF vs Control (Biopsy samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not significant" = "grey" ))

ggsave("Regulation summary plot_IPFvsControl_Biopsy.png", width = 10, height = 10)

# Create bar plot of regulation summary but excluding non signicant genes
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +
  geom_text(data = subset(summary_df, Category != "Not significant"),
            aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Summary of gene regulation in IPF vs Control (Biopsy samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3"))

ggsave("Regulation_summary_IPFvsControl_Biopsy_2.png", height = 10, width = 10)

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
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_IPF vs Control_Biopsy Samples",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_IPFvs_Ctrl_Biopsy_samples.png", vplot, width = 12, height = 10)
##End



# DGE analysis of CHP vs Control biopsy samples
##Start

# Define the column indices to check and the values to filter for
gene_counts_5 <- gene_counts
columns_to_check <- c(1:64) # Check all 64 columns
values_to_keep <- c("CHP", "CONTROL") # Values to look for

filtered_pheno5 <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno5$Platform <- filtered_pheno5$platform_id
rownames(filtered_pheno5) <- NULL
filtered_pheno5 <- filtered_pheno5[,c(1,13,17)]

# Filter out the explant samples
filtered_pheno5 <- filtered_pheno5[filtered_pheno5$`Sample type` == "Biopsy",]

# Group encoding
filtered_pheno5$Group_Encoded = ifelse(filtered_pheno5$Group == "CONTROL", 0, 1)
# Since gene_counts_5 has title as column names, filtered_pheno5$title is used to filter the date for gene_counts_5 
gene_counts_5 <- gene_counts_5[, colnames(gene_counts_5) %in% filtered_pheno5$title] 

# Matching order of column names of counts with row names of filtered pheno data
rownames(filtered_pheno5) <- filtered_pheno5$title
filtered_pheno5 <- filtered_pheno5[match(colnames(gene_counts_5), rownames(filtered_pheno5)), , drop = FALSE]

# Create sample info
sample_info_5 <- data.frame(
  sample = colnames(gene_counts_5),
  condition = factor(filtered_pheno5$Group)
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info_5) <- sample_info_5$sample

# Make CONTROL the first level since R by default will take CHP as reference level
sample_info_5$condition <- factor(sample_info_5$condition, levels = c("CONTROL", "CHP"))

# Create a DESeq2 object
dds = DESeqDataSetFromMatrix(
  countData = gene_counts_5,
  colData = sample_info_5,
  design = ~ condition
)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, name= "condition_CHP_vs_CONTROL")

# Order according to adjusted p values
res_ordered <- res[order(res$padj), ]

# Summarize results
summary(res)

# MA plot
plotMA(res, ylim = c(-5,5))

#Export results with regulation information
res_df <- as.data.frame(res) %>%
  rownames_to_column("Gene_symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "CHPvsCtrl(Biopsy)_deseq2_results.csv", row.names = FALSE)

#Extract significantly expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "CHPvsCtrl(Biopsy)_sig_exp_genes_results.csv", row.names = FALSE)

# Summary statistics
total_genes = nrow(res_df)
sig_up = sum(res_df$regulation == "Upregulated")
sig_down = sum(res_df$regulation == "Downregulated")
sig_total = sig_up + sig_down

# Print summary statistics
cat("Summary of Differential Expression Analysis :\n ")
cat("Total number of genes analyzed = ",total_genes, "\n")
cat("Number of significantly upregulated genes = ",sig_up,"\n")
cat("Number of significantly downregulated genes = ",sig_down,"\n")
cat("Total number of significantly differentially expressed genes = ",sig_total,"\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Upregulated", "Downregulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c((sig_up/total_genes)*100, (sig_down/total_genes)*100, ((total_genes - sig_total)/total_genes)*100)
)

# Round up to 2 decimal places
summary_df$Percentage <- round(summary_df$Percentage,2)

# Print summary dataframe
print(summary_df)

# Export (optional)
write.csv(summary_df, "Summary of Regulation_CHPvsCtrl_Biopsy.csv", row.names = FALSE)

# Create bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Gene regulation summary in HP vs Control (Biopsy samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not significant" = "grey" ))

ggsave("Regulation summary plot_CHPvsControl_Biopsy.png", width = 10, height = 10)

# Create bar plot of regulation summary but excluding non signicant genes
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +
  geom_text(data = subset(summary_df, Category != "Not significant"),
            aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Summary of gene regulation in HP vs Control (Biopsy samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3"))

ggsave("Regulation_summary_CHPvsControl_Biopsy_2.png", height = 10, width = 10)

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
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_HP vs Control_Biopsy_Samples",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_HPvsCtrl_Biopsy_Samples.png", vplot, width = 12, height = 10)
##End




# DGE analysis of IPF vs CHP biopsy samples
##Start

# Define the column indices to check and the values to filter for
gene_counts_6 <- gene_counts
columns_to_check <- c(1:64) # Check all 64 columns
values_to_keep <- c("CHP", "IPF") # Values to look for

filtered_pheno6 <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno6$Platform <- filtered_pheno6$platform_id
rownames(filtered_pheno6) <- NULL
filtered_pheno6 <- filtered_pheno6[,c(1,13,17)]

# Filter out the explant samples
filtered_pheno6 <- filtered_pheno6[filtered_pheno6$`Sample type` == "Biopsy",]

# Group encoding
filtered_pheno6$Group_Encoded = ifelse(filtered_pheno6$Group == "IPF", 0, 1)
# Since gene_counts_5 has title as column names, filtered_pheno5$title is used to filter the date for gene_counts_5 
gene_counts_6 <- gene_counts_6[, colnames(gene_counts_6) %in% filtered_pheno6$title] 

# Matching order of column names of counts with row names of filtered pheno data
rownames(filtered_pheno6) <- filtered_pheno6$title
filtered_pheno6 <- filtered_pheno6[match(colnames(gene_counts_6), rownames(filtered_pheno6)), , drop = FALSE]

# Create sample info
sample_info_6 <- data.frame(
  sample = colnames(gene_counts_6),
  condition = factor(filtered_pheno6$Group)
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info_6) <- sample_info_6$sample

# Make IPF as reference level
sample_info_6$condition <- factor(sample_info_6$condition, levels = c("IPF", "CHP"))

# Create a DESeq2 object
dds = DESeqDataSetFromMatrix(
  countData = gene_counts_6,
  colData = sample_info_6,
  design = ~ condition
)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

#Filter weakly expressed genes
filter <- HTSFilter(dds, s.len = 100, plot = TRUE)$filteredData
class(filter)

dim(dds)
dim(filter)

# Extract results
res <- results(dds, independentFiltering = FALSE,name= "condition_CHP_vs_IPF")
head(res)

# Order according to adjusted p values
res_ordered <- res[order(res$padj), ]

# Summarize results
summary(res)

# MA plot
plotMA(res, ylim = c(-5,5))

#Export results with regulation information
res_df <- as.data.frame(res) %>%
  rownames_to_column("Gene_symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "CHPvsIPF(Biopsy)_deseq2_results.csv", row.names = FALSE)

#Extract significantly expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "CHPvsIPF(Biopsy)_sig_exp_genes_results.csv", row.names = FALSE)

# Summary statistics
total_genes = nrow(res_df)
sig_up = sum(res_df$regulation == "Upregulated")
sig_down = sum(res_df$regulation == "Downregulated")
sig_total = sig_up + sig_down

# Print summary statistics
cat("Summary of Differential Expression Analysis :\n ")
cat("Total number of genes analyzed = ",total_genes, "\n")
cat("Number of significantly upregulated genes = ",sig_up,"\n")
cat("Number of significantly downregulated genes = ",sig_down,"\n")
cat("Total number of significantly differentially expressed genes = ",sig_total,"\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Upregulated", "Downregulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c((sig_up/total_genes)*100, (sig_down/total_genes)*100, ((total_genes - sig_total)/total_genes)*100)
)

# Round up to 2 decimal places
summary_df$Percentage <- round(summary_df$Percentage,2)

# Print summary dataframe
print(summary_df)

# Export (optional)
write.csv(summary_df, "Summary of Regulation_CHPvsIPF_Biopsy.csv", row.names = FALSE)

# Create bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Gene regulation summary in HP vs IPF (Biopsy samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not significant" = "grey" ))

ggsave("Regulation summary plot_CHPvsIPF_Biopsy.png", width = 10, height = 10)

# Create bar plot of regulation summary but excluding non signicant genes
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +
  geom_text(data = subset(summary_df, Category != "Not significant"),
            aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Summary of gene regulation in HP vs Control (Biopsy samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3"))

ggsave("Regulation_summary_CHPvsIPF_Biopsy_2.png", height = 10, width = 10)

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
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_HP vs IPF_Biopsy_Samples",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_HPvsIPF_Biopsy_samples.png", vplot, width = 12, height = 10)
##End






# DGE analysis of ILD vs Control biopsy samples
##Start
# Define the column indices to check and the values to filter for
gene_counts.ild <- gene_counts
columns_to_check <- c(1:64) # Check all 64 columns
values_to_keep <- c("IPF", "CONTROL", "CHP") # Values to look for

filtered_pheno.ild <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno.ild$Platform <- filtered_pheno.ild$platform_id
rownames(filtered_pheno.ild) <- NULL
filtered_pheno.ild <- filtered_pheno.ild[,c(1,13,17)]

# Filter out the explant samples
filtered_pheno.ild <- filtered_pheno.ild[filtered_pheno.ild$`Sample type` == "Biopsy",]

filtered_pheno.ild <- filtered_pheno.ild %>%
  mutate(Group = case_when(
    Group == "IPF" ~ "ILD",
    Group == "CHP" ~ "ILD",
    TRUE ~ "CONTROL"
  ))
# Group encoding
filtered_pheno.ild$Group_Encoded = ifelse(filtered_pheno.ild$Group == "CONTROL", 0, 1)
# Since gene_counts_4 has title as column names, filtered_pheno4$title is used to filter the date for gene_counts_4 
gene_counts.ild <- gene_counts.ild[, colnames(gene_counts.ild) %in% filtered_pheno.ild$title] 

# Matching order of column names of counts with row names of filtered pheno data
rownames(filtered_pheno.ild) <- filtered_pheno.ild$title
filtered_pheno.ild <- filtered_pheno.ild[match(colnames(gene_counts.ild), rownames(filtered_pheno.ild)), , drop = FALSE]

# Create sample info
sample_info.ild <- data.frame(
  sample = colnames(gene_counts.ild),
  condition = factor(filtered_pheno.ild$Group)
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info.ild) <- sample_info.ild$sample

# Create a DESeq2 object
dds = DESeqDataSetFromMatrix(
  countData = gene_counts.ild,
  colData = sample_info.ild,
  design = ~ condition
)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, name= "condition_ILD_vs_CONTROL")

# Order according to adjusted p values
res_ordered <- res[order(res$padj), ]

# Summarize results
summary(res)

# MA plot
plotMA(res, ylim = c(-5,5))

#Export results with regulation information
res_df <- as.data.frame(res) %>%
  rownames_to_column("Gene_symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "ILDvsCtrl(Biopsy)_deseq2_results.csv", row.names = FALSE)

#Extract significantly expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "ILDvsCtrl(Biopsy)_sig_exp_genes_results.csv", row.names = FALSE)

# Summary statistics
total_genes = nrow(res_df)
sig_up = sum(res_df$regulation == "Upregulated")
sig_down = sum(res_df$regulation == "Downregulated")
sig_total = sig_up + sig_down

# Print summary statistics
cat("Summary of Differential Expression Analysis :\n ")
cat("Total number of genes analyzed = ",total_genes, "\n")
cat("Number of significantly upregulated genes = ",sig_up,"\n")
cat("Number of significantly downregulated genes = ",sig_down,"\n")
cat("Total number of significantly differentially expressed genes = ",sig_total,"\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Upregulated", "Downregulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c((sig_up/total_genes)*100, (sig_down/total_genes)*100, ((total_genes - sig_total)/total_genes)*100)
)

# Round up to 2 decimal places
summary_df$Percentage <- round(summary_df$Percentage,2)

# Print summary dataframe
print(summary_df)

# Export (optional)
write.csv(summary_df, "Summary of Regulation_ILDvsCtrl_Biopsy.csv", row.names = FALSE)

# Create bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Gene regulation summary in ILD vs Control (Biopsy samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not significant" = "grey" ))

ggsave("Regulation summary plot_ILDvsControl_Biopsy.png", width = 10, height = 10)

# Create bar plot of regulation summary but excluding non signicant genes
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +
  geom_text(data = subset(summary_df, Category != "Not significant"),
            aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Summary of gene regulation in ILD vs Control (Biopsy samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3"))

ggsave("Regulation_summary_ILDvsControl_Biopsy_2.png", height = 10, width = 10)

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
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_ILD vs Control_Biopsy_Samples",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_ILDvsCtrl_Biopsy_samples.png", vplot, width = 12, height = 10)
##End







#EXPLANT SAMPLES
setwd("C:/IIT KGP/GSE150910/Explant samples")

# 1) IPF vs Control

##START
# Define the column indices to check and the values to filter for
gene_counts_7 <- gene_counts
columns_to_check <- c(1:64) # Check all 64 columns
values_to_keep <- c("IPF", "CONTROL") # Values to look for

filtered_pheno7 <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno7$Platform <- filtered_pheno7$platform_id
rownames(filtered_pheno7) <- NULL
filtered_pheno7 <- filtered_pheno7[,c(1,13,17)]

# Filter out the biopsy samples
filtered_pheno7 <- filtered_pheno7[filtered_pheno7$`Sample type` == "Explant",]

# Group encoding
filtered_pheno7$Group_Encoded = ifelse(filtered_pheno7$Group == "CONTROL", 0, 1)
# Since gene_counts_7 has title as column names, filtered_pheno7$title is used to filter the date for gene_counts_7
gene_counts_7 <- gene_counts_7[, colnames(gene_counts_7) %in% filtered_pheno7$title] 

# Matching order of column names of counts with row names of filtered pheno data
rownames(filtered_pheno7) <- filtered_pheno7$title
filtered_pheno7 <- filtered_pheno7[match(colnames(gene_counts_7), rownames(filtered_pheno7)), , drop = FALSE]

# Create sample info
sample_info_7 <- data.frame(
  sample = colnames(gene_counts_7),
  condition = factor(filtered_pheno7$Group)
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info_7) <- sample_info_7$sample

# Create a DESeq2 object
dds = DESeqDataSetFromMatrix(
  countData = gene_counts_7,
  colData = sample_info_7,
  design = ~ condition
)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, name= "condition_IPF_vs_CONTROL")

# Order according to adjusted p values
res_ordered <- res[order(res$padj), ]

# Summarize results
summary(res)

# MA plot
plotMA(res, ylim = c(-5,5))

#Export results with regulation information
res_df <- as.data.frame(res) %>%
  rownames_to_column("Gene_symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "IPFvsCtrl(Explant)_deseq2_results.csv", row.names = FALSE)

#Extract significantly expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "IPFvsCtrl(Explant)_sig_exp_genes_results.csv", row.names = FALSE)

# Summary statistics
total_genes = nrow(res_df)
sig_up = sum(res_df$regulation == "Upregulated")
sig_down = sum(res_df$regulation == "Downregulated")
sig_total = sig_up + sig_down

# Print summary statistics
cat("Summary of Differential Expression Analysis :\n ")
cat("Total number of genes analyzed = ",total_genes, "\n")
cat("Number of significantly upregulated genes = ",sig_up,"\n")
cat("Number of significantly downregulated genes = ",sig_down,"\n")
cat("Total number of significantly differentially expressed genes = ",sig_total,"\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Upregulated", "Downregulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c((sig_up/total_genes)*100, (sig_down/total_genes)*100, ((total_genes - sig_total)/total_genes)*100)
)

# Round up to 2 decimal places
summary_df$Percentage <- round(summary_df$Percentage,2)

# Print summary dataframe
print(summary_df)

# Export (optional)
write.csv(summary_df, "Summary of Regulation_IPFvsCtrl_Explant.csv", row.names = FALSE)

# Create bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Gene regulation summary in IPF vs Control (Explant samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not significant" = "grey" ))

ggsave("Regulation summary plot_IPFvsControl_Explant.png", width = 10, height = 10)

# Create bar plot of regulation summary but excluding non signicant genes
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +
  geom_text(data = subset(summary_df, Category != "Not significant"),
            aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Summary of gene regulation in IPF vs Control (Explant samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3"))

ggsave("Regulation_summary_IPFvsControl_Explant_2.png", height = 10, width = 10)

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
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_IPF vs Control_Explant_Samples",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_IPFvsCtrl_Explant_Samples.png", vplot, width = 12, height = 10)
##End


# 2) CHP vs Control

##START
# Define the column indices to check and the values to filter for
gene_counts_8 <- gene_counts
columns_to_check <- c(1:64) # Check all 64 columns
values_to_keep <- c("CHP", "CONTROL") # Values to look for

filtered_pheno8 <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno8$Platform <- filtered_pheno8$platform_id
rownames(filtered_pheno8) <- NULL
filtered_pheno8 <- filtered_pheno8[,c(1,13,17)]

# Filter out the biopsy samples
filtered_pheno8 <- filtered_pheno8[filtered_pheno8$`Sample type` == "Explant",]

# Group encoding
filtered_pheno8$Group_Encoded = ifelse(filtered_pheno8$Group == "CONTROL", 0, 1)
# Since gene_counts_7 has title as column names, filtered_pheno7$title is used to filter the date for gene_counts_7
gene_counts_8 <- gene_counts_8[, colnames(gene_counts_8) %in% filtered_pheno8$title] 

# Matching order of column names of counts with row names of filtered pheno data
rownames(filtered_pheno8) <- filtered_pheno8$title
filtered_pheno8 <- filtered_pheno8[match(colnames(gene_counts_8), rownames(filtered_pheno8)), , drop = FALSE]

# Create sample info
sample_info_8 <- data.frame(
  sample = colnames(gene_counts_8),
  condition = factor(filtered_pheno8$Group)
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info_8) <- sample_info_8$sample

# Make CONTROL the first level since R by default will take CHP as reference level
sample_info_8$condition <- factor(sample_info_8$condition, levels = c("CONTROL", "CHP"))

# Create a DESeq2 object
dds = DESeqDataSetFromMatrix(
  countData = gene_counts_8,
  colData = sample_info_8,
  design = ~ condition
)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, name= "condition_CHP_vs_CONTROL")

# Order according to adjusted p values
res_ordered <- res[order(res$padj), ]

# Summarize results
summary(res)

# MA plot
plotMA(res, ylim = c(-5,5))

#Export results with regulation information
res_df <- as.data.frame(res) %>%
  rownames_to_column("Gene_symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "CHPvsCtrl(Explant)_deseq2_results.csv", row.names = FALSE)

#Extract significantly expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "CHPvsCtrl(Explant)_sig_exp_genes_results.csv", row.names = FALSE)

# Summary statistics
total_genes = nrow(res_df)
sig_up = sum(res_df$regulation == "Upregulated")
sig_down = sum(res_df$regulation == "Downregulated")
sig_total = sig_up + sig_down

# Print summary statistics
cat("Summary of Differential Expression Analysis :\n ")
cat("Total number of genes analyzed = ",total_genes, "\n")
cat("Number of significantly upregulated genes = ",sig_up,"\n")
cat("Number of significantly downregulated genes = ",sig_down,"\n")
cat("Total number of significantly differentially expressed genes = ",sig_total,"\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Upregulated", "Downregulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c((sig_up/total_genes)*100, (sig_down/total_genes)*100, ((total_genes - sig_total)/total_genes)*100)
)

# Round up to 2 decimal places
summary_df$Percentage <- round(summary_df$Percentage,2)

# Print summary dataframe
print(summary_df)

# Export (optional)
write.csv(summary_df, "Summary of Regulation_CHPvsCtrl_Explant.csv", row.names = FALSE)

# Create bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Gene regulation summary in HP vs Control (Explant samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not significant" = "grey" ))

ggsave("Regulation summary plot_CHPvsControl_Explant.png", width = 10, height = 10)

# Create bar plot of regulation summary but excluding non signicant genes
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +
  geom_text(data = subset(summary_df, Category != "Not significant"),
            aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Summary of gene regulation in HP vs Control (Explant samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3"))

ggsave("Regulation_summary_CHPvsControl_Explant_2.png", height = 10, width = 10)

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
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_HP vs Control_Explant_samples",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_HPvsCtrl_Explant_samples.png", vplot, width = 12, height = 10)
##End


# 3) CHP vs IPF

##START
# Define the column indices to check and the values to filter for
gene_counts_9 <- gene_counts
columns_to_check <- c(1:64) # Check all 64 columns
values_to_keep <- c("CHP", "IPF") # Values to look for

filtered_pheno9 <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno9$Platform <- filtered_pheno9$platform_id
rownames(filtered_pheno9) <- NULL
filtered_pheno9 <- filtered_pheno9[,c(1,13,17)]

# Filter out the biopsy samples
filtered_pheno9 <- filtered_pheno9[filtered_pheno9$`Sample type` == "Explant",]

# Group encoding
filtered_pheno9$Group_Encoded = ifelse(filtered_pheno9$Group == "IPF", 0, 1)
# Since gene_counts_9 has title as column names, filtered_pheno9$title is used to filter the date for gene_counts_9
gene_counts_9 <- gene_counts_9[, colnames(gene_counts_9) %in% filtered_pheno9$title] 

# Matching order of column names of counts with row names of filtered pheno data
rownames(filtered_pheno9) <- filtered_pheno9$title
filtered_pheno9 <- filtered_pheno9[match(colnames(gene_counts_9), rownames(filtered_pheno9)), , drop = FALSE]

# Create sample info
sample_info_9 <- data.frame(
  sample = colnames(gene_counts_9),
  condition = factor(filtered_pheno9$Group)
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info_9) <- sample_info_9$sample

# Make IPF as reference level
sample_info_9$condition <- factor(sample_info_9$condition, levels = c("IPF", "CHP"))

# Create a DESeq2 object
dds = DESeqDataSetFromMatrix(
  countData = gene_counts_9,
  colData = sample_info_9,
  design = ~ condition
)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

#Filter out lowly expressed genes
filter <- HTSFilter(dds, s.len = 100, plot = FALSE)$filteredData
class(filter)

dim(dds)
dim(filter)

# Extract results
res <- results(dds, independentFiltering = FALSE, name= "condition_CHP_vs_IPF")

# Order according to adjusted p values
res_ordered <- res[order(res$padj), ]

# Summarize results
summary(res)

# MA plot
plotMA(res, ylim = c(-5,5))

#Export results with regulation information
res_df <- as.data.frame(res) %>%
  rownames_to_column("Gene_symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "CHPvsIPF(Explant)_deseq2_results.csv", row.names = FALSE)

#Extract significantly expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "CHPvsIPF(Explant)_sig_exp_genes_results.csv", row.names = FALSE)

# Summary statistics
total_genes = nrow(res_df)
sig_up = sum(res_df$regulation == "Upregulated")
sig_down = sum(res_df$regulation == "Downregulated")
sig_total = sig_up + sig_down

# Print summary statistics
cat("Summary of Differential Expression Analysis :\n ")
cat("Total number of genes analyzed = ",total_genes, "\n")
cat("Number of significantly upregulated genes = ",sig_up,"\n")
cat("Number of significantly downregulated genes = ",sig_down,"\n")
cat("Total number of significantly differentially expressed genes = ",sig_total,"\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Upregulated", "Downregulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c((sig_up/total_genes)*100, (sig_down/total_genes)*100, ((total_genes - sig_total)/total_genes)*100)
)

# Round up to 2 decimal places
summary_df$Percentage <- round(summary_df$Percentage,2)

# Print summary dataframe
print(summary_df)

# Export (optional)
write.csv(summary_df, "Summary of Regulation_CHPvsIPF_Explant.csv", row.names = FALSE)

# Create bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Gene regulation summary in HP vs IPF (Explant samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not significant" = "grey" ))

ggsave("Regulation summary plot_CHPvsIPF_Explant.png", width = 10, height = 10)

# Create bar plot of regulation summary but excluding non signicant genes
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +
  geom_text(data = subset(summary_df, Category != "Not significant"),
            aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Summary of gene regulation in HP vs IPF (Explant samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3"))

ggsave("Regulation_summary_CHPvsIPF_Explant_2.png", height = 10, width = 10)

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
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_HP vs IPF_Explant_Samples",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_HPvsIPF_Explant_Samples.png", vplot, width = 12, height = 10)
##End





# DGE analysis of ILD vs Control explant samples
##Start
# Define the column indices to check and the values to filter for
gene_counts.ild <- gene_counts
columns_to_check <- c(1:64) # Check all 64 columns
values_to_keep <- c("IPF", "CONTROL", "CHP") # Values to look for

filtered_pheno.ild <- filter_rows_with_found_column(pheno.data, columns_to_check, values_to_keep)

filtered_pheno.ild$Platform <- filtered_pheno.ild$platform_id
rownames(filtered_pheno.ild) <- NULL
filtered_pheno.ild <- filtered_pheno.ild[,c(1,13,17)]

# Filter out the explant samples
filtered_pheno.ild <- filtered_pheno.ild[filtered_pheno.ild$`Sample type` == "Explant",]

filtered_pheno.ild <- filtered_pheno.ild %>%
  mutate(Group = case_when(
    Group == "IPF" ~ "ILD",
    Group == "CHP" ~ "ILD",
    TRUE ~ "CONTROL"
  ))
# Group encoding
filtered_pheno.ild$Group_Encoded = ifelse(filtered_pheno.ild$Group == "CONTROL", 0, 1)
# Since gene_counts_4 has title as column names, filtered_pheno4$title is used to filter the date for gene_counts_4 
gene_counts.ild <- gene_counts.ild[, colnames(gene_counts.ild) %in% filtered_pheno.ild$title] 

# Matching order of column names of counts with row names of filtered pheno data
rownames(filtered_pheno.ild) <- filtered_pheno.ild$title
filtered_pheno.ild <- filtered_pheno.ild[match(colnames(gene_counts.ild), rownames(filtered_pheno.ild)), , drop = FALSE]

# Create sample info
sample_info.ild <- data.frame(
  sample = colnames(gene_counts.ild),
  condition = factor(filtered_pheno.ild$Group)
)

# Ensure row names of sample_info match column names of counts
rownames(sample_info.ild) <- sample_info.ild$sample

# Create a DESeq2 object
dds = DESeqDataSetFromMatrix(
  countData = gene_counts.ild,
  colData = sample_info.ild,
  design = ~ condition
)

# Step 3: Filter out low count genes
# Keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Extract results
res <- results(dds, name= "condition_ILD_vs_CONTROL")

# Order according to adjusted p values
res_ordered <- res[order(res$padj), ]

# Summarize results
summary(res)

# MA plot
plotMA(res, ylim = c(-5,5))

#Export results with regulation information
res_df <- as.data.frame(res) %>%
  rownames_to_column("Gene_symbol") %>%
  mutate(
    significant = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Yes", "No"),
    regulation = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

write.csv(res_df, "ILDvsCtrl(Explant)_deseq2_results.csv", row.names = FALSE)

#Extract significantly expressed genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_genes, "ILDvsCtrl(Explant)_sig_exp_genes_results.csv", row.names = FALSE)

# Summary statistics
total_genes = nrow(res_df)
sig_up = sum(res_df$regulation == "Upregulated")
sig_down = sum(res_df$regulation == "Downregulated")
sig_total = sig_up + sig_down

# Print summary statistics
cat("Summary of Differential Expression Analysis :\n ")
cat("Total number of genes analyzed = ",total_genes, "\n")
cat("Number of significantly upregulated genes = ",sig_up,"\n")
cat("Number of significantly downregulated genes = ",sig_down,"\n")
cat("Total number of significantly differentially expressed genes = ",sig_total,"\n")

# Create a summary dataframe
summary_df <- data.frame(
  Category = c("Upregulated", "Downregulated", "Not significant"),
  Count = c(sig_up, sig_down, total_genes - sig_total),
  Percentage = c((sig_up/total_genes)*100, (sig_down/total_genes)*100, ((total_genes - sig_total)/total_genes)*100)
)

# Round up to 2 decimal places
summary_df$Percentage <- round(summary_df$Percentage,2)

# Print summary dataframe
print(summary_df)

# Export (optional)
write.csv(summary_df, "Summary of Regulation_ILDvsCtrl_Explant.csv", row.names = FALSE)

# Create bar plot of regulation summary
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Gene regulation summary in ILD vs Control (Explant samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not significant" = "grey" ))

ggsave("Regulation summary plot_ILDvsControl_Explant.png", width = 10, height = 10)

# Create bar plot of regulation summary but excluding non signicant genes
ggplot(summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", data = subset(summary_df, Category != "Not significant")) +
  geom_text(data = subset(summary_df, Category != "Not significant"),
            aes(label = paste0(Count, " (", Percentage, "%)")),
            position = position_stack(vjust = 0.5),
            size = 5) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = "Summary of gene regulation in ILD vs Control (Explant samples)",
       x = "Regulation category",
       y = "Number of genes (count)") +
  scale_fill_manual(values = c("Upregulated" = "red3", "Downregulated" = "blue3"))

ggsave("Regulation_summary_ILDvsControl_Explant_2.png", height = 10, width = 10)

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
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Gene_symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_ILD vs Control_Explant_Samples",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_ILDvsCtrl_Explant_Samples.png", vplot, width = 12, height = 10)
##End








### edgeR pipeline

# Set working directory
setwd("C:/Users/Admin/Desktop/Abhijit/GSE150910/edgeR pipeline")

#Import library
library(edgeR)
library(pheatmap)
library(gridExtra)
library(tidyverse)
library(ggplot2)
library(ggrepel)

# Prepare the data
counts <- read.csv("C:/Users/Admin/Desktop/Abhijit/GSE150910/GSE150910_gene_counts.csv", row.names = 1)

# Get pheno data
gse <- getGEO("GSE150910", GSEMatrix = TRUE)
pheno.data <- gse$GSE150910_series_matrix.txt.gz@phenoData@data

pheno.data$`diagnosis:ch1` <- ifelse(
  pheno.data$`diagnosis:ch1` == "chp", "CHP",
  ifelse(pheno.data$`diagnosis:ch1`== "ipf", "IPF", "CONTROL")
)

filtered_pheno_data <- data.frame(
  Sample <- pheno.data$title,
  Group <- pheno.data$`diagnosis:ch1`
)

# Modify column names
colnames(filtered_pheno_data)[colnames(filtered_pheno_data) == "Sample....pheno.data.title"] <- "Sample"
colnames(filtered_pheno_data)[colnames(filtered_pheno_data) == "Group....pheno.data..diagnosis.ch1."] <- "Group"

# Match the row names with column names of counts
filtered_pheno_data <- filtered_pheno_data[match(colnames(counts),filtered_pheno_data$Sample),]
rownames(filtered_pheno_data) <- NULL

### DATA is now prepared

# Create group information
group <- factor(filtered_pheno_data$Group)

# Create DGEList object 
dge <- DGEList(counts = counts, group = group, remove.zeros = TRUE) #remove genes with 0 counts across all samples

# Filter lowly expressed genes
keep <- filterByExpr(dge, min.total.count = 10) 
dge <- dge[keep, , keep.lib.sizes=FALSE]
# each kept gene is required to have at least min.total.count reads across all the samples. 
# genes having total count >= 10 are kept

# Normalize (TMM normalization)
dge <- calcNormFactors(dge)

# Visualize MDS plot (optional, checks batch effect/clustering)
plotMDS(dge, col=as.numeric(group), cex=1.2, pch=16)
legend("topleft", legend=levels(group), col=1:length(levels(group)), pch=16)

# Create design matrix for multi group comparisons
design <- model.matrix(~0 + group) # "~0" removes intercept so you get estimates for all groups
colnames(design) <- levels(group)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit the data to glmQLFit model
fit_qlf <- glmQLFit(dge, design)

# Create contrasts (to compare)
contrasts <- makeContrasts(
  IPFvsCTRL = IPF - CONTROL,
  CHPvsCTRL = CHP - CONTROL,
  CHPvsIPF = CHP - IPF,
  levels = design)


## Get QLF result
qlf_ipf_ctrl <- glmQLFTest(fit_qlf, contrast = contrasts[,"IPFvsCTRL"])
qlf_chp_ctrl <- glmQLFTest(fit_qlf, contrast = contrasts[,"CHPvsCTRL"])
qlf_chp_ipf <- glmQLFTest(fit_qlf, contrast = contrasts[,"CHPvsIPF"])







### Contrast 1 - IPF vs Control

# Create and export summary
deg_results_ipf_ctrl <- topTags(qlf_ipf_ctrl, n=Inf)
deg_table <- deg_results_ipf_ctrl$table
deg_table <- deg_table %>%
  mutate(
    significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"),
    regulation = case_when(
      FDR < 0.05 & logFC > 1 ~ "Upregulated",
      FDR < 0.05 & logFC < -1 ~ "Downregulated",
      TRUE ~ "Not significant")
  )
write.csv(deg_table, "IPF_vs_Control_Summary_DEG_QLF.csv")

# Extract significantly expressed genes
deg_table_sig <- deg_table %>%
  filter(FDR < 0.05 & abs(logFC) > 1)
write.csv(deg_table_sig, "IPF_vs_Control_Significant_DEG_QLF.csv")

#MA plot

qlf_ipf_ctrl_deg <- topTags(qlf_ipf_ctrl, n=Inf)$table

# Define DE tags based on FDR threshold and logFC
upregulated <- rownames(qlf_ipf_ctrl_deg)[qlf_ipf_ctrl_deg$FDR < 0.05 & qlf_ipf_ctrl_deg$logFC > 1]
downregulated <- rownames(qlf_ipf_ctrl_deg)[qlf_ipf_ctrl_deg$FDR < 0.05 & qlf_ipf_ctrl_deg$logFC < -1]

# Base MA plot
plotSmear(qlf_ipf_ctrl, de.tags=NULL, main="MA Plot: IPF vs Control")

# Add points for upregulated and downregulated genes
points(qlf_ipf_ctrl$table$logCPM[rownames(qlf_ipf_ctrl$table) %in% upregulated],
       qlf_ipf_ctrl$table$logFC[rownames(qlf_ipf_ctrl$table) %in% upregulated],
       col="red3", pch=16, cex = 0.5)
points(qlf_ipf_ctrl$table$logCPM[rownames(qlf_ipf_ctrl$table) %in% downregulated],
       qlf_ipf_ctrl$table$logFC[rownames(qlf_ipf_ctrl$table) %in% downregulated],
       col="blue3", pch=16, cex = 0.5)

# Add horizontal lines for fold change threshold
abline(h=c(-1,1), col="darkgreen", lty=2)

# Add legend
legend("topright", legend=c("Upregulated", "Downregulated"),
       col=c("red3", "blue3"), pch=16)

# Volcano plot for IPF vs Control
# Set thresholds
logFC_threshold <- 1
FDR_threshold <- 0.05

# Add expression category
res_df <- qlf_ipf_ctrl_deg %>%
  mutate(expression = case_when(
    FDR < FDR_threshold & logFC > logFC_threshold ~ "Upregulated",
    FDR < FDR_threshold & logFC < -logFC_threshold ~ "Downregulated",
    TRUE ~ "Not Significant"
  ),
  gene = rownames(qlf_ipf_ctrl_deg)
  )

# Select top 15 up and top 15 down regulated genes by FDR
top_up <- res_df %>% 
  filter(expression == "Upregulated") %>%
  arrange(FDR) %>%
  slice(1:15)

top_down <- res_df %>%
  filter(expression == "Downregulated") %>%
  arrange(FDR) %>%
  slice(1:15)

top_genes <- bind_rows(top_up, top_down)

# Volcano plot
volcano_plot <- ggplot(res_df, aes(x = logFC, y = -log10(FDR), color = expression)) +
  geom_point(alpha = 0.8, size = 2.5) +
  geom_text_repel(data = top_genes,
                  aes(x = logFC, y = -log10(FDR),label = gene),
                  inherit.aes = FALSE,
                  size = 3.5,
                  segment.color = NA) + 
  scale_color_manual(name = "Regulation", 
                     values = c("Upregulated" = "red3", 
                                "Downregulated" = "blue3", 
                                "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(FDR_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot: IPF vs Control",
       x = "Log2 Fold Change",
       y = "-Log10 FDR") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("IPF_vs_Control_CustomVolcanoPlot.png", volcano_plot, width = 12, height = 10, dpi = 300)

############ Heatmap

# Filter normalized expression matrix (logCPM)
dge_norm <- cpm(dge, log=TRUE)
samples_ipf_ctrl <- which(group %in% c("IPF", "CONTROL"))
dge_norm_subset <- dge_norm[, samples_ipf_ctrl]

deg_table_sig_top.50 <- deg_table %>%
  filter(FDR < 0.05 & abs(logFC) > 1) %>%
  arrange(FDR) %>%
  head(50)  # Top 50 DEGs by FDR

# Match rows in expression matrix
heatmap_data <- dge_norm_subset[rownames(deg_table_sig_top.50), ]

# Sample annotation
sample_ann <- data.frame(condition = group[samples_ipf_ctrl])
rownames(sample_ann) <- colnames(heatmap_data)

# Gene annotation
gene_annotation <- data.frame(Regulation = deg_table_sig$regulation)
rownames(gene_annotation) <- rownames(deg_table_sig)

# Custom color palette
heatmap_colors <- colorRampPalette(c("#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026"))(100)

sample_colors <- c("CONTROL" = "grey", "IPF" = "orangered")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  condition = sample_colors,
  Regulation = gene_colors
)

# Plot heatmap
heatmap <- pheatmap(heatmap_data,
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    show_rownames = TRUE,
                    clustering_distance_cols = "euclidean",
                    clustering_distance_rows = "euclidean",
                    clustering_method = "complete",
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    show_colnames = FALSE,
                    scale = "row",
                    fontsize = 6,
                    main = "Heatmap of DEGs: IPF vs Control")

ggsave("IPFvsControl_Annotated_Heatmap.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)





### Contrast 2 - CHP vs Control

# Create and export summary
deg_results_chp_ctrl <- topTags(qlf_chp_ctrl, n=Inf)
deg_table <- deg_results_chp_ctrl$table
deg_table <- deg_table %>%
  mutate(
    significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"),
    regulation = case_when(
      FDR < 0.05 & logFC > 1 ~ "Upregulated",
      FDR < 0.05 & logFC < -1 ~ "Downregulated",
      TRUE ~ "Not significant")
  )
write.csv(deg_table, "CHP_vs_Control_Summary_DEG_QLF.csv")

# Extract significantly expressed genes
deg_table_sig <- deg_table %>%
  filter(FDR < 0.05 & abs(logFC) > 1)
write.csv(deg_table_sig, "CHP_vs_Control_Significant_DEG_QLF.csv")

#MA plot

qlf_chp_ctrl_deg <- topTags(qlf_chp_ctrl, n=Inf)$table

# Define DE tags based on FDR threshold and logFC
upregulated <- rownames(qlf_chp_ctrl_deg)[qlf_chp_ctrl_deg$FDR < 0.05 & qlf_chp_ctrl_deg$logFC > 1]
downregulated <- rownames(qlf_chp_ctrl_deg)[qlf_chp_ctrl_deg$FDR < 0.05 & qlf_chp_ctrl_deg$logFC < -1]

# Base MA plot
plotSmear(qlf_chp_ctrl, de.tags=NULL, main="MA Plot: HP vs Control")

# Add points for upregulated and downregulated genes
points(qlf_chp_ctrl$table$logCPM[rownames(qlf_chp_ctrl$table) %in% upregulated],
       qlf_chp_ctrl$table$logFC[rownames(qlf_chp_ctrl$table) %in% upregulated],
       col="red3", pch=16, cex=0.5)
points(qlf_chp_ctrl$table$logCPM[rownames(qlf_chp_ctrl$table) %in% downregulated],
       qlf_chp_ctrl$table$logFC[rownames(qlf_chp_ctrl$table) %in% downregulated],
       col="blue3", pch=16, cex=0.5)

# Add horizontal lines for fold change threshold
abline(h=c(-1,1), col="darkgreen", lty=2)

# Add legend
legend("topright", legend=c("Upregulated", "Downregulated"),
       col=c("red3", "blue3"), pch=16)

# Volcano plot for IPF vs Control
# Set thresholds
logFC_threshold <- 1
FDR_threshold <- 0.05

# Add expression category
res_df <- qlf_chp_ctrl_deg %>%
  mutate(expression = case_when(
    FDR < FDR_threshold & logFC > logFC_threshold ~ "Upregulated",
    FDR < FDR_threshold & logFC < -logFC_threshold ~ "Downregulated",
    TRUE ~ "Not Significant"
  ),
  gene = rownames(qlf_chp_ctrl_deg)
  )

# Select top 15 up and top 15 down regulated genes by FDR
top_up <- res_df %>% 
  filter(expression == "Upregulated") %>%
  arrange(FDR) %>%
  slice(1:15)

top_down <- res_df %>%
  filter(expression == "Downregulated") %>%
  arrange(FDR) %>%
  slice(1:15)

top_genes <- bind_rows(top_up, top_down)

# Volcano plot
volcano_plot <- ggplot(res_df, aes(x = logFC, y = -log10(FDR), color = expression)) +
  geom_point(alpha = 0.8, size = 2.5) +
  geom_text_repel(data = top_genes,
                  aes(x = logFC, y = -log10(FDR), label = gene),
                  inherit.aes = FALSE,
                  size = 3.5,
                  segment.color = NA) + 
  scale_color_manual(name = "Regulation", 
                     values = c("Upregulated" = "red3", 
                                "Downregulated" = "blue3", 
                                "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(FDR_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot: HP vs Control",
       x = "Log2 Fold Change",
       y = "-Log10 FDR") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("HP_vs_Control_CustomVolcanoPlot.png", volcano_plot, width = 12, height = 10, dpi = 300)


############ Heatmap

# Filter normalized expression matrix (logCPM)
dge_norm <- cpm(dge, log=TRUE)
samples_chp_ctrl <- which(group %in% c("CHP", "CONTROL"))
dge_norm_subset <- dge_norm[, samples_chp_ctrl]

deg_table_sig_top.50 <- deg_table_sig %>%
  filter(FDR < 0.05 & abs(logFC) > 1) %>%
  arrange(FDR) %>%
  head(50)  # Top 50 DEGs by FDR

# Match rows in expression matrix
heatmap_data <- dge_norm_subset[rownames(deg_table_sig_top.50), ]

# Sample annotation
sample_ann <- data.frame(condition = group[samples_chp_ctrl])
rownames(sample_ann) <- colnames(heatmap_data)

# Gene annotation
gene_annotation <- data.frame(Regulation = deg_table_sig$regulation)
rownames(gene_annotation) <- rownames(deg_table_sig)

# Custom color palette
heatmap_colors <- colorRampPalette(c("#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026"))(100)

sample_colors <- c("CONTROL" = "grey", "CHP" = "cadetblue")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  condition = sample_colors,
  Regulation = gene_colors
)

# Plot heatmap
heatmap <- pheatmap(heatmap_data,
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    show_rownames = TRUE,
                    clustering_distance_cols = "euclidean",
                    clustering_distance_rows = "euclidean",
                    clustering_method = "complete",
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    #show_colnames = FALSE,
                    scale = "row",
                    fontsize = 6,
                    main = "Heatmap of DEGs: IPF vs Control")

ggsave("CHPvsControl_Annotated_Heatmap_withsamples.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)





### Contrast 3 - CHP vs IPF

# Create and export summary
deg_results_chp_ipf <- topTags(qlf_chp_ipf, n=Inf)
deg_table <- deg_results_chp_ipf$table
deg_table <- deg_table %>%
  mutate(
    significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"),
    regulation = case_when(
      FDR < 0.05 & logFC > 1 ~ "Upregulated",
      FDR < 0.05 & logFC < -1 ~ "Downregulated",
      TRUE ~ "Not significant")
  )
write.csv(deg_table, "CHP_vs_IPF_Summary_DEG_QLF.csv")

# Extract significantly expressed genes
deg_table_sig <- deg_table %>%
  filter(FDR < 0.05 & abs(logFC) > 1)
write.csv(deg_table_sig, "CHP_vs_IPF_Significant_DEG_QLF.csv")

#MA plot

qlf_chp_ipf_deg <- topTags(qlf_chp_ipf, n=Inf)$table

# Define DE tags based on FDR threshold and logFC
upregulated <- rownames(qlf_chp_ipf_deg)[qlf_chp_ipf_deg$FDR < 0.05 & qlf_chp_ipf_deg$logFC > 1]
downregulated <- rownames(qlf_chp_ipf_deg)[qlf_chp_ipf_deg$FDR < 0.05 & qlf_chp_ipf_deg$logFC < -1]

# Base MA plot
plotSmear(qlf_chp_ipf, de.tags=NULL, main="MA Plot: HP vs IPF")

# Add points for upregulated and downregulated genes
points(qlf_chp_ipf$table$logCPM[rownames(qlf_chp_ipf$table) %in% upregulated],
       qlf_chp_ipf$table$logFC[rownames(qlf_chp_ipf$table) %in% upregulated],
       col="red3", pch=16, cex=0.5 )
points(qlf_chp_ipf$table$logCPM[rownames(qlf_chp_ipf$table) %in% downregulated],
       qlf_chp_ipf$table$logFC[rownames(qlf_chp_ipf$table) %in% downregulated],
       col="blue3", pch=16, cex=0.5)

# Add horizontal lines for fold change threshold
abline(h=c(-1,1), col="darkgreen", lty=2)

# Add legend
legend("topright", legend=c("Upregulated", "Downregulated"),
       col=c("red3", "blue3"), pch=16)

# Volcano plot for CHP vs IPF
# Set thresholds
logFC_threshold <- 1
FDR_threshold <- 0.05

# Add expression category
res_df <- qlf_chp_ipf_deg %>%
  mutate(expression = case_when(
    FDR < FDR_threshold & logFC > logFC_threshold ~ "Upregulated",
    FDR < FDR_threshold & logFC < -logFC_threshold ~ "Downregulated",
    TRUE ~ "Not Significant"
  ),
  gene = rownames(qlf_chp_ipf_deg)
  )

# Select top 15 up and top 15 down regulated genes by FDR
top_up <- res_df %>% 
  filter(expression == "Upregulated") %>%
  arrange(FDR) %>%
  slice(1:15)

top_down <- res_df %>%
  filter(expression == "Downregulated") %>%
  arrange(FDR) %>%
  slice(1:15)

top_genes <- bind_rows(top_up, top_down)

# Volcano plot
volcano_plot <- ggplot(res_df, aes(x = logFC, y = -log10(FDR), color = expression)) +
  geom_point(alpha = 0.8, size = 2.5) +
  geom_text_repel(data = top_genes,
                  aes(x = logFC, y = -log10(FDR), label = gene),
                  inherit.aes = FALSE,
                  size = 3.5,
                  segment.color = NA) + 
  scale_color_manual(name = "Regulation", 
                     values = c("Upregulated" = "red3", 
                                "Downregulated" = "blue3", 
                                "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(FDR_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot: HP vs IPF",
       x = "Log2 Fold Change",
       y = "-Log10 FDR") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("HP_vs_IPF_CustomVolcanoPlot.png", volcano_plot, width = 12, height = 10, dpi = 300)

############ Heatmap

# Filter normalized expression matrix (logCPM)
dge_norm <- cpm(dge, log=TRUE)
samples_chp_ipf <- which(group %in% c("CHP", "IPF"))
dge_norm_subset <- dge_norm[, samples_chp_ipf]

deg_table_sig_top.50 <- deg_table %>%
  filter(FDR < 0.05 & abs(logFC) > 1) %>%
  arrange(FDR) %>%
  head(50)  # Top 50 DEGs by FDR

# Match rows in expression matrix
heatmap_data <- dge_norm_subset[rownames(deg_table_sig_top.50), ]

# Sample annotation
sample_ann <- data.frame(condition = group[samples_chp_ipf])
rownames(sample_ann) <- colnames(heatmap_data)

# Gene annotation
gene_annotation <- data.frame(Regulation = deg_table_sig$regulation)
rownames(gene_annotation) <- rownames(deg_table_sig)

# Custom color palette
heatmap_colors <- colorRampPalette(c("#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026"))(100)

sample_colors <- c("CHP" = "cadetblue", "IPF" = "orangered")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  condition = sample_colors,
  Regulation = gene_colors
)

# Plot heatmap
heatmap <- pheatmap(heatmap_data,
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    show_rownames = TRUE,
                    clustering_distance_cols = "euclidean",
                    clustering_distance_rows = "euclidean",
                    clustering_method = "complete",
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    show_colnames = FALSE,
                    scale = "row",
                    fontsize = 6,
                    main = "Heatmap of DEGs: CHP vs IPF")

ggsave("CHPvsIPF_Annotated_Heatmap_withsamples.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)
#END


#### ILD vs CONTROL
filtered_pheno_data <- filtered_pheno_data %>%
  mutate(condition = ifelse(Group == "CONTROL", "CONTROL", "ILD"))

#Create group information
group.ild <- factor(filtered_pheno_data$condition, levels = c("CONTROL", "ILD"))

# Create DGEList object 
dge <- DGEList(counts = counts, group = group.ild, remove.zeros = TRUE) #remove genes with 0 counts across all samples

# Filter lowly expressed genes
keep <- filterByExpr(dge, min.total.count = 10) 
dge <- dge[keep, , keep.lib.sizes=FALSE]
# each kept gene is required to have at least min.total.count reads across all the samples. 
# genes having total count >= 10 are kept

# Normalize (TMM normalization)
dge <- calcNormFactors(dge)

# Visualize MDS plot (optional, checks batch effect/clustering)
plotMDS(dge, col=as.numeric(group), cex=1.2, pch=16)
legend("topleft", legend=levels(group), col=1:length(levels(group)), pch=16)

# Create design matrix for multi group comparisons
design <- model.matrix(~0 + group.ild) # "~0" removes intercept so you get estimates for all groups
colnames(design) <- levels(group.ild)

# Estimate dispersion
dge <- estimateDisp(dge, design)

# Fit the data to glmQLFit model
fit_qlf.ild <- glmQLFit(dge, design)

# Create contrasts (to compare)
contrasts <- makeContrasts(
  ILDvsCTRL = ILD - CONTROL,
  levels = design)


## Get QLF result
qlf_ild_ctrl <- glmQLFTest(fit_qlf.ild, contrast = contrasts[,"ILDvsCTRL"])



### Contrast - ILD vs Control

# Create and export summary
deg_results_ild_ctrl <- topTags(qlf_ild_ctrl, n=Inf)
deg_table.ild <- deg_results_ild_ctrl$table
deg_table.ild <- deg_table.ild %>%
  mutate(
    significant = ifelse(FDR < 0.05 & abs(logFC) > 1, "Yes", "No"),
    regulation = case_when(
      FDR < 0.05 & logFC > 1 ~ "Upregulated",
      FDR < 0.05 & logFC < -1 ~ "Downregulated",
      TRUE ~ "Not significant")
  )
write.csv(deg_table.ild, "ILD_vs_Control_Summary_DEG_QLF.csv")

# Extract significantly expressed genes
deg_table_sig.ild <- deg_table.ild %>%
  filter(FDR < 0.05 & abs(logFC) > 1)
write.csv(deg_table_sig.ild, "ILD_vs_Control_Significant_DEG_QLF.csv")

#MA plot

qlf_ild_ctrl_deg <- topTags(qlf_ild_ctrl, n=Inf)$table

# Define DE tags based on FDR threshold and logFC
upregulated <- rownames(qlf_ild_ctrl_deg)[qlf_ild_ctrl_deg$FDR < 0.05 & qlf_ild_ctrl_deg$logFC > 1]
downregulated <- rownames(qlf_ild_ctrl_deg)[qlf_ild_ctrl_deg$FDR < 0.05 & qlf_ild_ctrl_deg$logFC < -1]

# Base MA plot
plotSmear(qlf_ild_ctrl, de.tags=NULL, main="MA Plot: ILD vs Control")

# Add points for upregulated and downregulated genes
points(qlf_ild_ctrl$table$logCPM[rownames(qlf_ild_ctrl$table) %in% upregulated],
       qlf_ild_ctrl$table$logFC[rownames(qlf_ild_ctrl$table) %in% upregulated],
       col="red3", pch=16, cex = 0.5)
points(qlf_ild_ctrl$table$logCPM[rownames(qlf_ild_ctrl$table) %in% downregulated],
       qlf_ild_ctrl$table$logFC[rownames(qlf_ild_ctrl$table) %in% downregulated],
       col="blue3", pch=16, cex = 0.5)

# Add horizontal lines for fold change threshold
abline(h=c(-1,1), col="darkgreen", lty=2)

# Add legend
legend("topright", legend=c("Upregulated", "Downregulated"),
       col=c("red3", "blue3"), pch=16)

# Volcano plot for IPF vs Control
# Set thresholds
logFC_threshold <- 1
FDR_threshold <- 0.05

# Add expression category
res_df <- qlf_ild_ctrl_deg %>%
  mutate(expression = case_when(
    FDR < FDR_threshold & logFC > logFC_threshold ~ "Upregulated",
    FDR < FDR_threshold & logFC < -logFC_threshold ~ "Downregulated",
    TRUE ~ "Not Significant"
  ),
  gene = rownames(qlf_ild_ctrl_deg)
  )

# Select top 15 up and top 15 down regulated genes by FDR
top_up <- res_df %>% 
  filter(expression == "Upregulated") %>%
  arrange(FDR) %>%
  slice(1:15)

top_down <- res_df %>%
  filter(expression == "Downregulated") %>%
  arrange(FDR) %>%
  slice(1:15)

top_genes <- bind_rows(top_up, top_down)

# Volcano plot
volcano_plot <- ggplot(res_df, aes(x = logFC, y = -log10(FDR), color = expression)) +
  geom_point(alpha = 0.8, size = 2.5) +
  geom_text_repel(data = top_genes,
                  aes(x = logFC, y = -log10(FDR), label = gene),
                  inherit.aes = FALSE,
                  size = 3.5,
                  segment.color = NA) + 
  scale_color_manual(name = "Regulation", 
                     values = c("Upregulated" = "red3", 
                                "Downregulated" = "blue3", 
                                "Not Significant" = "gray")) +
  geom_hline(yintercept = -log10(FDR_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot: ILD vs Control",
       x = "Log2 Fold Change",
       y = "-Log10 FDR") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("ILD_vs_Control_CustomVolcanoPlot.png", volcano_plot, width = 12, height = 10, dpi = 300)

############ Heatmap

# Filter normalized expression matrix (logCPM)
dge_norm <- cpm(dge, log=TRUE)
samples_ild_ctrl <- which(group.ild %in% c("ILD", "CONTROL"))
dge_norm_subset <- dge_norm[, samples_ild_ctrl]

deg_table_sig_top.50 <- deg_table.ild %>%
  filter(FDR < 0.05 & abs(logFC) > 1) %>%
  arrange(FDR) %>%
  head(50)  # Top 50 DEGs by FDR

# Match rows in expression matrix
heatmap_data <- dge_norm_subset[rownames(deg_table_sig_top.50), ]

# Sample annotation
sample_ann <- data.frame(condition = group.ild[samples_ild_ctrl])
rownames(sample_ann) <- colnames(heatmap_data)

# Gene annotation
gene_annotation <- data.frame(Regulation = deg_table_sig.ild$regulation)
rownames(gene_annotation) <- rownames(deg_table_sig.ild)

# Custom color palette
heatmap_colors <- colorRampPalette(c("#313695", "#4575B4", "#74ADD1", 
                                     "#FFFFFF", "#FDAE61", "#F46D43", "#A50026"))(100)

sample_colors <- c("CONTROL" = "grey", "ILD" = "red4")
gene_colors <- c("Upregulated" = "red3", "Downregulated" = "blue3")

annotation_colors <- list(
  condition = sample_colors,
  Regulation = gene_colors
)

# Plot heatmap
heatmap <- pheatmap(heatmap_data,
                    annotation_col = sample_ann,
                    annotation_row = gene_annotation,
                    annotation_colors = annotation_colors,
                    color = heatmap_colors,
                    show_rownames = TRUE,
                    clustering_distance_cols = "euclidean",
                    clustering_distance_rows = "euclidean",
                    clustering_method = "complete",
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    show_colnames = FALSE,
                    scale = "row",
                    fontsize = 6,
                    main = "Heatmap of DEGs: ILD vs Control")

ggsave("ILDvsControl_Annotated_Heatmap.png", plot = heatmap$gtable, width = 15, height = 10, dpi = 300)
