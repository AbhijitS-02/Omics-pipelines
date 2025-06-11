#######################################
########### META ANALYSIS #############
#######################################
##### FUSED INVERSE NORMAL METHOD #####
#######################################


# Set working directory
setwd("C:/IIT KGP")


## RNA-seq datasets and sample size of each dataset.
datasets <- c("GSE184316", "GSE150910_biopsy", "GSE150910_explant", "GSE47460_GPL14550")          
n_samp <- c(76, 62, 123, 143)
sig_cutoff <- 1e-16       # the lower limit of system cutoff for a number.


# Import all the required datasets
GSE184316.degs <- read.csv("C:/IIT KGP/GSE184316/HPvsIPF/HPvsIPF_deseq2_results.csv", header = TRUE)

GSE150910.biopsy.degs <- read.csv("C:/IIT KGP/GSE150910/Biopsy samples/CHP vs IPF/CHPvsIPF(Biopsy)_deseq2_results.csv",
                                  header = TRUE)

GSE150910.explant.degs <- read.csv("C:/IIT KGP/GSE150910/Explant samples/CHP vs IPF/CHPvsIPF(Explant)_deseq2_results.csv",
                                   header = TRUE)

GSE47460.GPL14550.degs <- read.csv("C:/IIT KGP/GSE47460/GPL14450/HPvsIPF/topTable result_GPL14550_HP_vs_IPF.csv",
                                   header = TRUE)



####### Define function for fused inverse normal method 

metaFIN <- function(study_1, study_2, study_3, study_4, n_samp, datasets){
  
  ## Determine direction of expression for the differentially expressed genes
  study_1$dir <- sign(study_1$logFC)
  study_2$dir <- sign(study_2$logFC)
  study_3$dir <- sign(study_3$logFC)
  study_4$dir <- sign(study_4$logFC)
  
  # Get all unique genes from the 4 studies
  unique_genes <- Reduce(union, list(study_1$Gene_symbol, 
                                     study_2$Gene_symbol, 
                                     study_3$Gene_symbol, 
                                     study_4$Gene_symbol))
  
  print(unique_genes)
  
  # Get direction status of expression of the unique genes to determine if they are conflicting or not
  sign <- matrix(data = 0, nrow = length(unique_genes), ncol = length(n_samp)+1)
  row.names(sign) <- unique_genes
  
  sign[,1] <- study_1$dir[match(unique_genes, study_1$Gene_symbol)]
  sign[,2] <- study_2$dir[match(unique_genes, study_2$Gene_symbol)]
  sign[,3] <- study_3$dir[match(unique_genes, study_3$Gene_symbol)]
  sign[,4] <- study_4$dir[match(unique_genes, study_4$Gene_symbol)]
  
  for(l in 1:length(unique_genes))
  {
    if (1 %in% sign[l, c(1:length(n_samp))] & -1 %in% sign[l, c(1:length(n_samp))])
    {
      sign[l, (length(n_samp)+1)] <- 1
    }
  }
  
  write.csv(as.data.frame(sign), "Sign.csv")
  ### Calculation of Ng terms
  
  # Step 1 : Estimation of weights
  
  ## Initialize an empty weights matrix with the number of rows equal to the number of unique genes and
  ## columns equal to the number of studies
  weights <- matrix(0, nrow = length(unique_genes), ncol = length(n_samp))
  
  # each element in weight matrix corresponds to a unique gene 
  # numerator terms of w_s
  weights[which((unique_genes %in% study_1$Gene_symbol) == TRUE), 1] <- n_samp[1]
  weights[which((unique_genes %in% study_2$Gene_symbol) == TRUE), 2] <- n_samp[2]
  weights[which((unique_genes %in% study_3$Gene_symbol) == TRUE), 3] <- n_samp[3]
  weights[which((unique_genes %in% study_4$Gene_symbol) == TRUE), 4] <- n_samp[4]
  
  # denominator terms of w_s
  denom <- apply(weights, 1, sum)
  
  # divide numerator by denominator and square root to get final weights
  weights <- weights/denom
  weights <- sqrt(weights)
  row.names(weights) <- as.character(unique_genes)
  colnames(weights) <- c("study_1", "study_2", "study_3", "study_4")
  
  weights <- as.data.frame(weights, stringsAsFactors = FALSE)
  write.csv(weights, file = "weights.csv") #save weights as csv - optional
  
  
  # Step 2 : Calculation of each term of Ng for a gene g
  ng_terms <- matrix(0, nrow = nrow(weights), ncol = ncol(weights))
  
  for(j in 1:nrow(ng_terms))
  {
    #check if a gene has same direction of expression across studies
    #if yes, then use first case definition for Ng
    
    if (sign[j, ncol(sign)] == 0) #for non conflicting genes
    {
      
      #dataset 1
      if (unique_genes[j] %in% study_1$Gene_symbol)
      {
        k = which(study_1$Gene_symbol == unique_genes[j])
        p_val = min(max(study_1$Pvalue[k],1e-16), 1 - 1e-16)
        ng_terms[j,1] <- weights$study_1[j] * qnorm((1 - p_val), mean = 0, sd = 1)
      }
      
      #dataset 2
      if (unique_genes[j] %in% study_2$Gene_symbol)
      {
        k = which(study_2$Gene_symbol == unique_genes[j])
        p_val = min(max(study_2$Pvalue[k],1e-16), 1 - 1e-16)
        ng_terms[j,2] <- weights$study_2[j] * qnorm((1 - p_val), mean = 0, sd = 1)
      }
      
      #dataset 3
      if (unique_genes[j] %in% study_3$Gene_symbol)
      {
        k = which(study_3$Gene_symbol == unique_genes[j])
        p_val = min(max(study_3$Pvalue[k],1e-16), 1 - 1e-16)
        ng_terms[j,3] <- weights$study_3[j] * qnorm((1 - p_val), mean = 0, sd = 1)
      }
      
      #dataset 1
      if (unique_genes[j] %in% study_4$Gene_symbol)
      {
        k = which(study_4$Gene_symbol == unique_genes[j])
        p_val = min(max(study_4$Pvalue[k],1e-16), 1 - 1e-16)
        ng_terms[j,4] <- weights$study_4[j] * qnorm((1 - p_val), mean = 0, sd = 1)
      }
    }
    
    #check if a gene has conflicting direction of expression across studies 
    #if yes, use second case definition of Ng
    
    if(sign[j, ncol(sign)] == 1)
    {
      
      #dataset 1
      if(unique_genes[j] %in% study_1$Gene_symbol)
      {
        k = which(study_1$Gene_symbol == unique_genes[j])
        p_val = min(max(study_1$Pvalue[k], sig_cutoff), 1 - sig_cutoff)
        ber.rand.var = sign(study_1$logFC[j]) #Bernoulli random variable, only takes 1 and -1 as values
        ng_terms[j,1] <- weights$study_1[j] * ber.rand.var * abs(qnorm((1 - p_val), mean = 0, sd = 1))
      }
      
      #dataset 2
      if(unique_genes[j] %in% study_2$Gene_symbol)
      {
        k = which(study_2$Gene_symbol == unique_genes[j])
        p_val = min(max(study_2$Pvalue[k], sig_cutoff), 1 - sig_cutoff)
        ber.rand.var = sign(study_2$logFC[j]) #Bernoulli random variable, only takes 1 and -1 as values
        ng_terms[j,2] <- weights$study_2[j] * ber.rand.var * abs(qnorm((1 - p_val), mean = 0, sd = 1))
      }
      
      #dataset 3
      if(unique_genes[j] %in% study_3$Gene_symbol)
      {
        k = which(study_3$Gene_symbol == unique_genes[j])
        p_val = min(max(study_3$Pvalue[k], sig_cutoff), 1 - sig_cutoff)
        ber.rand.var = sign(study_3$logFC[j]) #Bernoulli random variable, only takes 1 and -1 as values
        ng_terms[j,3] <- weights$study_3[j] * ber.rand.var * abs(qnorm((1 - p_val), mean = 0, sd = 1))
      }
      
      #dataset 4
      if(unique_genes[j] %in% study_4$Gene_symbol)
      {
        k = which(study_4$Gene_symbol == unique_genes[j])
        p_val = min(max(study_4$Pvalue[k], sig_cutoff), 1 - sig_cutoff)
        ber.rand.var = sign(study_4$logFC[j]) #Bernoulli random variable, only takes 1 and -1 as values
        ng_terms[j,4] <- weights$study_4[j] * ber.rand.var * abs(qnorm((1 - p_val), mean = 0, sd = 1))
      }
    }
  }
  
  colnames(ng_terms) <- datasets
  ng_terms <- as.data.frame(ng_terms, stringsAsFactors = FALSE)
  row.names(ng_terms) <- row.names(weights)
  
  
  # Step 3 : Sum all the Ng terms row-wise
  
  ng <- as.data.frame(rowSums(ng_terms))
  colnames(ng) <- c("ng")
  row.names(ng) <- row.names(ng_terms)
  
  # Step 4 : Hypothesis testing
  
  #One sided test on right-hand tail of the distribution performed for genes with same direction of expression
  #Two sided test performed for genes with conflicting direction of expression
  
  ## First do one sided test for all genes, then replace with two sided test for conflicting genes
  ng$metaFIN_pval <- 1 - pnorm(ng$ng)
  
  conflict_index <- which(sign[, ncol(sign)] == 1) #index of conflicting genes
  ng$metaFIN_pval[conflict_index] <- 2 * (1 - pnorm(abs(ng$ng[conflict_index]))) #replaced by two-sided test
  
  # FDR correction for multiple hypothesis testing using Benjamini-Hochberg method
  ng$metaFIN_padj <- p.adjust(ng$metaFIN_pval, method = "BH", n = length(ng$metaFIN_pval))
  
  return(ng)
} 


# Prepare the datasets for proper input format
## Input format:  Gene_symbol | logFC | Pvalue

# Dataset 1
head(GSE184316.degs)

study.1_GSE184316 <- GSE184316.degs %>%
  select(Gene_Symbol, log2FoldChange, pvalue)

colnames(study.1_GSE184316) <- c("Gene_symbol", "logFC", "Pvalue")




#Dataset 2
head(GSE150910.biopsy.degs)

study.2_GSE150910.biopsy <- GSE150910.biopsy.degs %>%
  select(Gene_symbol, log2FoldChange, pvalue)

colnames(study.2_GSE150910.biopsy) <- c("Gene_symbol", "logFC", "Pvalue")



#Dataset 3
head(GSE150910.explant.degs)

study.3_GSE150910.explant <- GSE150910.explant.degs %>%
  select(Gene_symbol, log2FoldChange, pvalue)

colnames(study.3_GSE150910.explant) <- c("Gene_symbol", "logFC", "Pvalue")



#Dataset 4
head(GSE47460.GPL14550.degs)

study.4_GSE47460 <- GSE47460.GPL14550.degs %>%
  select(Gene_Symbol, logFC, P.Value)

colnames(study.4_GSE47460) <- c("Gene_symbol", "logFC", "Pvalue")

##Datasets are prepared

### Call the metaFIN function to run Fused Inverse Normal Method

meta_FIN_results <- metaFIN(study_1 = study.1_GSE184316,
                            study_2 = study.2_GSE150910.biopsy,
                            study_3 = study.3_GSE150910.explant,
                            study_4 = study.4_GSE47460,
                            n_samp = n_samp,
                            datasets = datasets)


meta_FIN_results <- meta_FIN_results %>%
  rownames_to_column(var = "Gene_symbol")

write.csv(meta_FIN_results, "FINM_Meta Analysis Results.csv", row.names = FALSE)


### Calculate average logFC for all the unique genes

head(meta_FIN_results)


## Dataset 1 logFC values
# Find matching indices in both data frames
gene_index <- match(meta_FIN_results$Gene_symbol, study.1_GSE184316$Gene_symbol)

# Update the GSE184316 column with corresponding logFC values where matches exist
meta_FIN_results$GSE184316 <- study.1_GSE184316$logFC[gene_index]




## Dataset 2 logFC
# Find matching indices in both data frames
gene_index <- match(meta_FIN_results$Gene_symbol, study.2_GSE150910.biopsy$Gene_symbol)

# Update the GSE184316 column with corresponding logFC values where matches exist
meta_FIN_results$GSE150910.biopsy <- study.2_GSE150910.biopsy$logFC[gene_index]





## Dataset 3 logFC
# Find matching indices in both data frames
gene_index <- match(meta_FIN_results$Gene_symbol, study.3_GSE150910.explant$Gene_symbol)

# Update the GSE184316 column with corresponding logFC values where matches exist
meta_FIN_results$GSE150910.explant <- study.3_GSE150910.explant$logFC[gene_index]




## Dataset 4 logFC
# Find matching indices in both data frames
gene_index <- match(meta_FIN_results$Gene_symbol, study.4_GSE47460$Gene_symbol)

# Update the GSE184316 column with corresponding logFC values where matches exist
meta_FIN_results$GSE47460 <- study.4_GSE47460$logFC[gene_index]



# Calculate mean absolute log fold change

for(l in 1:nrow(meta_FIN_results)) #iterate through all the genes
{
  logfc.sum = 0
  counter <- 0 #counter to keep track of non-NaN values to calculate average as 
               #dividing by NaN value will give NaN as output
  
  ## eg: if 2 values are NaN and 2 are not, sum will be divided only by 2 and not by 4 to get average
  
  if (is.na(meta_FIN_results$GSE184316[l]) == FALSE) #NaN values will be skipped through
  {
    counter = counter + 1
    logfc.sum = logfc.sum + meta_FIN_results$GSE184316[l]
  }
  if(is.na(meta_FIN_results$GSE150910.biopsy[l]) == FALSE) #NaN values will be skipped through
  {
    counter = counter + 1
    logfc.sum = logfc.sum + meta_FIN_results$GSE150910.biopsy[l]
  }
  if(is.na(meta_FIN_results$GSE150910.explant[l]) == FALSE) #NaN values will be skipped through
  {
    counter = counter + 1
    logfc.sum = logfc.sum + meta_FIN_results$GSE150910.explant[l]
  }
  if(is.na(meta_FIN_results$GSE47460[l]) == FALSE) #NaN values will be skipped through
  {
    counter = counter + 1
    logfc.sum = logfc.sum + meta_FIN_results$GSE47460[l]
  }
  
  meta_FIN_results$mean_abs_logfc[l] <- abs(logfc.sum/counter) #get absolute average mean log fold change
}

# Rearrange the columns

meta_FIN_results <- meta_FIN_results[, c("Gene_symbol", "ng", "mean_abs_logfc", "GSE184316",
                                         "GSE150910.biopsy", "GSE150910.explant", "GSE47460",
                                         "metaFIN_pval", "metaFIN_padj")]

meta_FINM_degs <- meta_FIN_results %>%
  filter(mean_abs_logfc > 1 & metaFIN_padj < 0.05) %>%
  arrange(desc(ng))

write.csv(meta_FINM_degs, "FINM_significant DEGs.csv", row.names = FALSE)


# Assess directionality of expression
direction.expr <- read.csv("C:/IIT KGP/Sign.csv")

