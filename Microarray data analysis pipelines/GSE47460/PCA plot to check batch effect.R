library(ggplot2)
library(limma)  # or limma for microarray
library(ggfortify)

# Perform PCA before batch correction
pca_before <- prcomp(t(combined_expr), scale. = TRUE)

# Perform PCA after batch correction
pca_after <- prcomp(t(combat_exprs), scale. = TRUE)

# Convert to data frame
pca_df_before <- as.data.frame(pca_before$x)
pca_df_after <- as.data.frame(pca_after$x)

# Add sample metadata
pca_df_before$Batch <- platform_info
pca_df_before$Condition <- combined.sample.info$condition
pca_df_after$Batch <- platform_info
pca_df_after$Condition <- combined.sample.info$condition

# Plot PCA
before.plot <- ggplot(pca_df_before, aes(x = PC1, y = PC2, color = Batch, shape = Condition)) +
  geom_point(size = 3) +
  ggtitle("PCA Before ComBat") +
  theme_minimal()

after.plot <- ggplot(pca_df_after, aes(x = PC1, y = PC2, color = Batch, shape = Condition)) +
  geom_point(size = 3) +
  ggtitle("PCA After ComBat") +
  theme_minimal()


