##### HP vs CONTROL

# Create a customized volcano plot
# Set significance thresholds
logFC_threshold <- 1
pvalue_threshold <- 0.05

# Create the expression classification column
res.HP.GPL14550$expression <- "Not Significant"
res.HP.GPL14550$expression[res.HP.GPL14550$logFC > logFC_threshold & res.HP.GPL14550$adj.P.Val < pvalue_threshold] <- "Upregulated"
res.HP.GPL14550$expression[res.HP.GPL14550$logFC < -logFC_threshold & res.HP.GPL14550$adj.P.Val < pvalue_threshold] <- "Downregulated"

# Filter and select top 10 upregulated genes
top_up <- res.HP.GPL14550[res.HP.GPL14550$expression == "Upregulated" & res.HP.GPL14550$adj.P.Val < pvalue_threshold, ]
top_up <- top_up[order(top_up$adj.P.Val), ][1:15, ]

top_up <- top_up %>%
  rownames_to_column(var = "Gene_Symbol")

# Filter and select top 10 downregulated genes
top_down <- res.HP.GPL14550[res.HP.GPL14550$expression == "Downregulated" & res.HP.GPL14550$adj.P.Val < pvalue_threshold, ]
top_down <- top_down[order(top_down$adj.P.Val), ][1:15, ]

top_down <- top_down %>%
  rownames_to_column(var = "Gene_Symbol")

# Combine top genes
top_genes <- rbind(top_up, top_down)

# Create the volcano plot
vplot <- ggplot(res.HP.GPL14550, aes(x = logFC, y = -log10(adj.P.Val), color = expression)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(name = "Regulation",
                     values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not Significant" = "gray")) +
  theme_bw() +
  geom_text_repel(data = top_genes,
                  aes(x = logFC, y = -log10(adj.P.Val), label = Gene_Symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_HP vs Control_GPL14550",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_HPvsCtrl_GPL14550.png", vplot, width = 12, height = 10)





##### IPF vs CONTROL

# Create a customized volcano plot
# Set significance thresholds
logFC_threshold <- 1
pvalue_threshold <- 0.05

# Create the expression classification column
res.IPF.GPL14550$expression <- "Not Significant"
res.IPF.GPL14550$expression[res.IPF.GPL14550$logFC > logFC_threshold & res.IPF.GPL14550$adj.P.Val < pvalue_threshold] <- "Upregulated"
res.IPF.GPL14550$expression[res.IPF.GPL14550$logFC < -logFC_threshold & res.IPF.GPL14550$adj.P.Val < pvalue_threshold] <- "Downregulated"

# Filter and select top 10 upregulated genes
top_up <- res.IPF.GPL14550[res.IPF.GPL14550$expression == "Upregulated" & res.IPF.GPL14550$adj.P.Val < pvalue_threshold, ]
top_up <- top_up[order(top_up$adj.P.Val), ][1:15, ]

top_up <- top_up %>%
  rownames_to_column(var = "Gene_Symbol")

# Filter and select top 10 downregulated genes
top_down <- res.IPF.GPL14550[res.IPF.GPL14550$expression == "Downregulated" & res.IPF.GPL14550$adj.P.Val < pvalue_threshold, ]
top_down <- top_down[order(top_down$adj.P.Val), ][1:15, ]

top_down <- top_down %>%
  rownames_to_column(var = "Gene_Symbol")

# Combine top genes
top_genes <- rbind(top_up, top_down)

# Create the volcano plot
vplot <- ggplot(res.IPF.GPL14550, aes(x = logFC, y = -log10(adj.P.Val), color = expression)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(name = "Regulation",
                     values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not Significant" = "gray")) +
  theme_bw() +
  geom_text_repel(data = top_genes,
                  aes(x = logFC, y = -log10(adj.P.Val), label = Gene_Symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_IPF vs Control_GPL14550",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_IPFvsCtrl_GPL14550.png", vplot, width = 12, height = 10)




##### ILD vs CONTROL

# Create a customized volcano plot
# Set significance thresholds
logFC_threshold <- 1
pvalue_threshold <- 0.05

# Create the expression classification column
res.ild.GPL14550$expression <- "Not Significant"
res.ild.GPL14550$expression[res.ild.GPL14550$logFC > logFC_threshold & res.ild.GPL14550$adj.P.Val < pvalue_threshold] <- "Upregulated"
res.ild.GPL14550$expression[res.ild.GPL14550$logFC < -logFC_threshold & res.ild.GPL14550$adj.P.Val < pvalue_threshold] <- "Downregulated"

# Filter and select top 10 upregulated genes
top_up <- res.ild.GPL14550[res.ild.GPL14550$expression == "Upregulated" & res.ild.GPL14550$adj.P.Val < pvalue_threshold, ]
top_up <- top_up[order(top_up$adj.P.Val), ][1:15, ]

top_up <- top_up %>%
  rownames_to_column(var = "Gene_Symbol")

# Filter and select top 10 downregulated genes
top_down <- res.ild.GPL14550[res.ild.GPL14550$expression == "Downregulated" & res.ild.GPL14550$adj.P.Val < pvalue_threshold, ]
top_down <- top_down[order(top_down$adj.P.Val), ][1:15, ]

top_down <- top_down %>%
  rownames_to_column(var = "Gene_Symbol")

# Combine top genes
top_genes <- rbind(top_up, top_down)

# Create the volcano plot
vplot <- ggplot(res.ild.GPL14550, aes(x = logFC, y = -log10(adj.P.Val), color = expression)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(name = "Regulation",
                     values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not Significant" = "gray")) +
  theme_bw() +
  geom_text_repel(data = top_genes,
                  aes(x = logFC, y = -log10(adj.P.Val), label = Gene_Symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_ILD vs Control_GPL14550",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_ILDvsCtrl_GPL14550.png", vplot, width = 12, height = 10)




##### HP vs IPF

# Create a customized volcano plot
# Set significance thresholds
logFC_threshold <- 1
pvalue_threshold <- 0.05

# Create the expression classification column
res.hp.ipf$expression <- "Not Significant"
res.hp.ipf$expression[res.hp.ipf$logFC > logFC_threshold & res.hp.ipf$adj.P.Val < pvalue_threshold] <- "Upregulated"
res.hp.ipf$expression[res.hp.ipf$logFC < -logFC_threshold & res.hp.ipf$adj.P.Val < pvalue_threshold] <- "Downregulated"

# Filter and select upregulated genes
top_up <- res.hp.ipf[res.hp.ipf$expression == "Upregulated" & res.hp.ipf$adj.P.Val < pvalue_threshold, ]
top_up <- top_up[order(top_up$adj.P.Val), ][1:15, ]
top_up <- top_up %>%
  rownames_to_column(var = "Gene_Symbol")

# Filter and select downregulated genes
top_down <- res.hp.ipf[res.hp.ipf$expression == "Downregulated" & res.hp.ipf$adj.P.Val < pvalue_threshold, ]
top_down <- top_down[order(top_down$adj.P.Val), ][1:15, ]
top_down <- top_down %>%
  rownames_to_column(var = "Gene_Symbol")

# Combine top genes
top_genes <- rbind(top_up, top_down)

# Create the volcano plot
vplot <- ggplot(res.hp.ipf, aes(x = logFC, y = -log10(adj.P.Val), color = expression)) +
  geom_point(alpha = 0.8, size = 3) +
  scale_color_manual(name = "Regulation",
                     values = c("Upregulated" = "red3", "Downregulated" = "blue3", "Not Significant" = "gray")) +
  theme_bw() +
  geom_text_repel(data = top_genes,
                  aes(x = logFC, y = -log10(adj.P.Val), label = Gene_Symbol),
                  inherit.aes = FALSE,
                  size = 4,
                  max.overlaps = Inf,
                  box.padding = 0.3,
                  point.padding = 0.3,
                  segment.color = NA) +
  labs(title = "Volcano Plot_HP vs IPF_GPL14550",
       x = "Log Fold Change (logFC)",
       y = "-Log10 P-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black")
ggsave("volcano_plot_HPvsIPF_GPL14550.png", vplot, width = 12, height = 10)
