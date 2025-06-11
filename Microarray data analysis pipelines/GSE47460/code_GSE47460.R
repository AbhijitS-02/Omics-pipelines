############################################
# GSE47460 <- TISSUE DATASET
############################################

# Set working directory
setwd("C:/IIT KGP/GSE47460")

#install libraries
BiocManager::install("sva")

# Load the libraries
library(GEOquery)
library(tidyverse)
library(limma)
library(sva)
library(stringr)

# Get supplementary files
getGEOSuppFiles("GSE47460")

# Untar the files
untar("GSE47460_RAW.tar", exdir = "data/")

# Get phenotypic data
gse <- getGEO("GSE47460", GSEMatrix = TRUE)

pheno.data.GPL14550 <- gse$`GSE47460-GPL14550_series_matrix.txt.gz`@phenoData@data
pheno.data.GPL6480 <- gse$`GSE47460-GPL6480_series_matrix.txt.gz`@phenoData@data






######### Filter the phenodata

### GPL14550
filtered.pheno.data.GPL14550 <- pheno.data.GPL14550 %>%
  filter(`disease state:ch1` == "Control" | `disease state:ch1` == "Interstitial lung disease")

# create a new column to separate IPF, HP and Control
filtered.pheno.data.GPL14550 <- filtered.pheno.data.GPL14550 %>%
  mutate(condition = case_when(
    grepl("IPF", `ild subtype:ch1`) == TRUE ~ "IPF",
    grepl("HP", `ild subtype:ch1`) == TRUE ~ "HP",
    is.na(`ild subtype:ch1`) ~ "Control",
    TRUE ~ `ild subtype:ch1`
  ))

# filter again
filtered.pheno.data.GPL14550 <- filtered.pheno.data.GPL14550 %>%
  filter(condition == "IPF" | condition == "HP" | condition == "Control")











### Repeat for GPL6480
filtered.pheno.data.GPL6480 <- pheno.data.GPL6480 %>%
  filter(`disease state:ch1` == "Control" | `disease state:ch1` == "Interstitial lung disease")

# create a new column to separate IPF, HP and Control
filtered.pheno.data.GPL6480 <- filtered.pheno.data.GPL6480 %>%
  mutate(condition = case_when(
    grepl("IPF", `ild subtype:ch1`) == TRUE ~ "IPF",
    grepl("HP", `ild subtype:ch1`) == TRUE ~ "HP",
    is.na(`ild subtype:ch1`) ~ "Control",
    TRUE ~ `ild subtype:ch1`
  ))

# filter again
filtered.pheno.data.GPL6480 <- filtered.pheno.data.GPL6480 %>%
  filter(condition == "IPF" | condition == "HP" | condition == "Control")








##### Separate the files according to Platform ID

# Filter GSMs with "ILD" or "CTRL" in the title
gsm_ids_GPL14550<- rownames(pheno.data.GPL14550)[grepl("ILD|CTRL", pheno.data.GPL14550$title)]

for (gsm in gsm_ids_GPL14550){
  
  # get the title to complete whole filename
  title <- pheno.data.GPL14550[gsm, "title"]
  
  filename <- paste0(gsm, "_", title, ".txt.gz")
  
  # Paths
  from <- file.path("data/", filename)
  to <- file.path("data/GPL14550", filename)
  
  if (file.exists(from)) {
    file.rename(from, to) #if file exists, copy to GPL14550 folder
  }else {
    paste0(filename, "is missing")
  } 
}


gsm_ids_GPL6480 <- rownames(pheno.data.GPL6480)[grepl("ILD|CTRL", pheno.data.GPL6480$title)]

for (gsm in gsm_ids_GPL6480){
  
  # get the title to complete whole filename
  title <- pheno.data.GPL6480[gsm, "title"]
  
  filename <- paste0(gsm, "_", title, ".txt.gz")
  
  # Paths
  from <- file.path("data/", filename)
  to <- file.path("data/GPL6480", filename)
  
  if (file.exists(from)) {
    file.rename(from, to) #if file exists, copy to GPL14550 folder
  }else {
    paste0(filename, "is missing")
  } 
}
########## All files are now separated based on platform and ILD and CONTROL




################ Filter the files to have only IPF, HP and Control patient data

### GPL14550
# Get all files in the directory
files <- list.files("data/GPL14550", full.names = TRUE)

# Extract GSM IDs from filenames (before first underscore)
file_gsm_ids <- gsub("_.*", "", basename(files))

# Get GSM IDs from phenotype data
gsm_ids <- rownames(filtered.pheno.data.GPL14550)

# Identify files to delete 
files_to_delete <- files[!file_gsm_ids %in% gsm_ids]

# Delete the files
file.remove(files_to_delete)



### GPL6480
# Get all files in the directory
files <- list.files("data/GPL6480", full.names = TRUE)

# Extract GSM IDs from filenames (before first underscore)
file_gsm_ids <- gsub("_.*", "", basename(files))

# Get GSM IDs from phenotype data
gsm_ids <- rownames(filtered.pheno.data.GPL6480)

# Identify files to delete 
files_to_delete <- files[!file_gsm_ids %in% gsm_ids]

# Delete the files
file.remove(files_to_delete)

##########################################



### Read the datasets

# For GPL14550
files.GPL14550 <- list.files("data/GPL14550", full.names = TRUE)
raw.GPL14550 <- read.maimages(files = files.GPL14550, source = "agilent", 
                              green.only = TRUE, other.columns = "gIsWellAboveBG")
dim(raw.GPL14550)

# For GPL6480
files.GPL6480 <- list.files("data/GPL6480", full.names = TRUE)
raw.GPL6480 <- read.maimages(files = files.GPL6480, source = "agilent", 
                             green.only = TRUE, other.columns = "gIsWellAboveBG")
dim(raw.GPL6480)

#################################################


#### Workflow to get Expression data of Platform -> GPL14550

# Gene annotation
annotation.GPL14550 <- read.delim("GPL14550_annotation.txt", sep = "\t", header = TRUE)

colnames(annotation.GPL14550)[colnames(annotation.GPL14550) == "ID"] <- "ProbeName"

# merge
raw.GPL14550$genes <- merge(annotation.GPL14550, raw.GPL14550$genes, by= "ProbeName")
raw.GPL14550$genes <- raw.GPL14550$genes %>%
  select(Row, Col, ControlType, ProbeName, SystematicName, GENE_SYMBOL, GENE_NAME)

## Background correction
raw.GPL14550.bg.corrected <- backgroundCorrect(raw.GPL14550, method = "normexp", offset = 40)

## Normalize
raw.GPL14550.bg.corrected <- normalizeBetweenArrays(raw.GPL14550.bg.corrected, 
                                                    method = "cyclicloess", 
                                                    cyclic.method = "pairs")

## Gene filtering

raw.GPL14550.bg.corrected$genes[raw.GPL14550.bg.corrected$genes == ""] <- NA

# Filter out Control probes
Control <- raw.GPL14550.bg.corrected$genes$ControlType == 1L

# Filter out probes with no symbol
NoSymbol <- is.na(raw.GPL14550.bg.corrected$genes$GENE_SYMBOL)

# Keep probes having above background expression in atleast 20% of arrays
isExpr <- (rowSums(raw.GPL14550.bg.corrected$other$gIsWellAboveBG > 0)) >= (0.2 * ncol(raw.GPL14550.bg.corrected$other$gIsWellAboveBG))


filtered.raw.GPL14550.bg.corrected <- raw.GPL14550.bg.corrected[!Control & !NoSymbol & isExpr,]
dim(filtered.raw.GPL14550.bg.corrected)

# Keep required columns
filtered.raw.GPL14550.bg.corrected$genes <- filtered.raw.GPL14550.bg.corrected$genes[, c("ProbeName", "GENE_SYMBOL", "SystematicName")]
head(filtered.raw.GPL14550.bg.corrected$genes)


## Probe summarisation
exprs.GPL14550 <- avereps(filtered.raw.GPL14550.bg.corrected$E, ID = filtered.raw.GPL14550.bg.corrected$genes$GENE_SYMBOL)


#####################################################################



#### Workflow to get Expression data of Platform -> GPL6480

# Gene annotation
annotation.GPL6480 <- read.delim("GPL6480_annotation.txt", sep = "\t", header = TRUE)

colnames(annotation.GPL6480)[colnames(annotation.GPL6480) == "ID"] <- "ProbeName"

# merge
raw.GPL6480$genes <- merge(annotation.GPL6480, raw.GPL6480$genes, by= "ProbeName")
raw.GPL6480$genes <- raw.GPL6480$genes %>%
  select(Row, Col, ControlType, ProbeName, SystematicName, GENE_SYMBOL, GENE_NAME)

## Background correction
raw.GPL6480.bg.corrected <- backgroundCorrect(raw.GPL6480, method = "normexp", offset = 35)

## Normalize
raw.GPL6480.bg.corrected <- normalizeBetweenArrays(raw.GPL6480.bg.corrected, 
                                                    method = "cyclicloess", 
                                                    cyclic.method = "pairs")

## Gene filtering

raw.GPL6480.bg.corrected$genes[raw.GPL6480.bg.corrected$genes == ""] <- NA

# Filter out Control probes
Control <- raw.GPL6480.bg.corrected$genes$ControlType == 1L

# Filter out probes with no symbol
NoSymbol <- is.na(raw.GPL6480.bg.corrected$genes$GENE_SYMBOL)

# Keep probes having above background expression in atleast 20% of arrays
isExpr <- (rowSums(raw.GPL6480.bg.corrected$other$gIsWellAboveBG > 0)) >= (0.2 * ncol(raw.GPL6480.bg.corrected$other$gIsWellAboveBG))


filtered.raw.GPL6480.bg.corrected <- raw.GPL6480.bg.corrected[!Control & !NoSymbol & isExpr,]
dim(filtered.raw.GPL6480.bg.corrected)

# Keep required columns
filtered.raw.GPL6480.bg.corrected$genes <- filtered.raw.GPL6480.bg.corrected$genes[, c("ProbeName", "GENE_SYMBOL", "SystematicName")]
head(filtered.raw.GPL6480.bg.corrected$genes)


## Probe summarisation
exprs.GPL6480 <- avereps(filtered.raw.GPL6480.bg.corrected$E, ID = filtered.raw.GPL6480.bg.corrected$genes$GENE_SYMBOL)



############################################
#For GPL14550
colnames(exprs.GPL14550) <- str_extract(colnames(exprs.GPL14550), "GSM[0-9]+")

# For GPL6480
colnames(exprs.GPL6480) <- str_extract(colnames(exprs.GPL6480), "GSM[0-9]+")


## Create sample columns

# GPL14550
sample.info.GPL14550 <- filtered.pheno.data.GPL14550 %>%
  select(geo_accession, condition)

# GPL6480
sample.info.GPL6480 <- filtered.pheno.data.GPL6480 %>%
  select(geo_accession, condition)



#####################################################
############### WORKFLOW 1 - MERGED #################
#####################################################

# Find common genes and merge the expression
common_genes <- intersect(rownames(exprs.GPL14550), rownames(exprs.GPL6480))

combined_expr <- cbind(exprs.GPL14550[common_genes,],
                       exprs.GPL6480[common_genes, ])

# Format column names
colnames(combined_expr) <- str_extract(colnames(combined_expr), "GSM[0-9]+")

# Create platform info for batch effect removal
platform_info <- c(rep("GPL14550", ncol(exprs.GPL14550)), rep("GPL6480", ncol(exprs.GPL6480)))

# Create combined sample info
combined.sample.info <- rbind(sample.info.GPL14550, sample.info.GPL6480)

# Check if order is correct
all(rownames(combined.sample.info) == colnames(combined_expr))

# Create design matrix
condition <- factor(combined.sample.info$condition, levels = c("Control", "IPF", "HP"))
design <- model.matrix(~ condition)

# Correct Platform effects
combat_exprs <- ComBat(dat = combined_expr, batch = platform_info, mod = design, par.prior = TRUE, prior.plots = TRUE)

# Fit linear model
fit <- lmFit(combat_exprs, design = design)
fitted.ebayes <- eBayes(fit, trend = TRUE, robust = TRUE)
summary(decideTests(fitted.ebayes[,-1], p.value = 0.05, lfc = 1))

#### 1) HP vs Control
res.HP <- topTable(fitted.ebayes, coef = "conditionHP", adjust.method = "BH", n = Inf)
degs.HP <- res.HP %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
  rownames_to_column(var = "Gene_Symbol")

write.csv(degs.HP, "DEGs_HP_vs_Ctrl.csv", row.names = FALSE)


#### 2) IPF vs Control
res.IPF <- topTable(fitted.ebayes, coef = "conditionIPF", adjust.method = "BH", n = Inf)
degs.IPF <- res.IPF %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
  rownames_to_column(var = "Gene_Symbol")

write.csv(degs.IPF, "DEGs_IPF_vs_Ctrl.csv", row.names = FALSE)



#### 3) ILD vs Control
sample.info <- merged.pheno.data %>%
  select(geo_accession, `disease state:ch1`)

colnames(sample.info)[colnames(sample.info) == "disease state:ch1"] <- "condition"

# Check if order is correct
all(rownames(sample.info) == colnames(combat_exprs))

# Create design matrix
condition <- factor(sample.info$condition, levels = c("Control", "Interstitial lung disease"))
design <- model.matrix(~ condition)

# Fit linear model
fit.ild <- lmFit(combat_exprs, design = design)
fit.ild <- eBayes(fit.ild, trend = TRUE, robust = TRUE)
summary(decideTests(fit.ild[,-1], p.value = 0.05, lfc = 1))

res.ild <- topTable(fit.ild, coef = "conditionInterstitial lung disease", adjust.method = "BH", n = Inf)
degs.ild <- res.ild %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
  rownames_to_column(var = "Gene_Symbol")

write.csv(degs.ild, "DEGs_ild_vs_Ctrl.csv", row.names = FALSE)




#### 4) HP vs IPF
sample.info.hp.ipf <- merged.pheno.data %>%
  filter(condition == "HP" | condition == "IPF") %>%
  select(geo_accession, condition) 
  
# Check if order is correct
combat_exprs <- combat_exprs[, colnames(combat_exprs) %in% rownames(sample.info.hp.ipf)]
all(rownames(sample.info.hp.ipf) == colnames(combat_exprs))

# Create design matrix
condition <- factor(sample.info.hp.ipf$condition, levels = c("IPF", "HP"))
design <- model.matrix(~ condition)

# Fit linear model
fit.hp.ipf <- lmFit(combat_exprs, design = design)
fit.hp.ipf <- eBayes(fit.hp.ipf, trend = TRUE, robust = TRUE)
summary(decideTests(fit.hp.ipf[,-1], p.value = 0.05, lfc = 1))

res.hp.ipf <- topTable(fit.hp.ipf, coef = "conditionHP", adjust.method = "BH", n = Inf)
degs.hp.ipf <- res.hp.ipf %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
  rownames_to_column(var = "Gene_Symbol")

write.csv(degs.hp.ipf, "DEGs_HP_vs_IPF.csv", row.names = FALSE)










#######################################################
############### WORKFLOW 2 - GPL14550 #################
#######################################################

# GPL14550 sample info
sample.info.GPL14550 <- filtered.pheno.data.GPL14550 %>%
  select(geo_accession, condition)

# Create design matrix
condition.GPL14550 <- factor(sample.info.GPL14550$condition, levels = c("Control", "IPF", "HP"))
design.GPL14550 <- model.matrix(~ condition.GPL14550)

# Fit linear model
fit.GPL14550 <- lmFit(exprs.GPL14550, design = design.GPL14550)
fitted.ebayes.GPL14550 <- eBayes(fit.GPL14550, trend = TRUE, robust = TRUE)
summary(decideTests(fitted.ebayes.GPL14550[,-1], p.value = 0.05, lfc = 1))

#### 1) HP vs Control
res.HP.GPL14550 <- topTable(fitted.ebayes.GPL14550, coef = "condition.GPL14550HP", 
                            adjust.method = "BH", n = Inf)
degs.HP.GPL14550 <- res.HP.GPL14550 %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
  rownames_to_column(var = "Gene_Symbol")

write.csv(degs.HP.GPL14550, "DEGs_HP_vs_Ctrl_GPL14550.csv", row.names = FALSE)


#### 2) IPF vs Control
res.IPF.GPL14550 <- topTable(fitted.ebayes.GPL14550, coef = "condition.GPL14550IPF", adjust.method = "BH", n = Inf)
degs.IPF.GPL14550 <- res.IPF.GPL14550 %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
  rownames_to_column(var = "Gene_Symbol")

write.csv(degs.IPF.GPL14550, "DEGs_IPF_vs_Ctrl_GPL14550.csv", row.names = FALSE)



#### 3) ILD vs Control
sample.info.GPL14550 <- sample.info.GPL14550 %>%
  mutate(disease = ifelse(condition == "Control", "Control", "ILD"))

# Check if order is correct
all(rownames(sample.info.GPL14550) == colnames(exprs.GPL14550))

# Create design matrix
condition <- factor(sample.info.GPL14550$disease, levels = c("Control", "ILD"))
design <- model.matrix(~ condition)

# Fit linear model
fit.ild.GPL14550 <- lmFit(exprs.GPL14550, design = design)
fit.ild.GPL14550 <- eBayes(fit.ild.GPL14550, trend = TRUE, robust = TRUE)
summary(decideTests(fit.ild.GPL14550[,-1], p.value = 0.05, lfc = 1))

res.ild.GPL14550 <- topTable(fit.ild.GPL14550, coef = "conditionILD", adjust.method = "BH", n = Inf)
degs.ild.GPL14550 <- res.ild.GPL14550 %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
  rownames_to_column(var = "Gene_Symbol")

write.csv(degs.ild.GPL14550, "DEGs_ILD_vs_Ctrl_GPL14550.csv", row.names = FALSE)




#### 4) HP vs IPF
sample.info.hp.ipf <- sample.info.GPL14550 %>%
  filter(condition == "HP" | condition == "IPF") %>%
  select(geo_accession, condition) 

# Check if order is correct
exprs.GPL14550 <- exprs.GPL14550[, colnames(exprs.GPL14550) %in% rownames(sample.info.hp.ipf)]
all(rownames(sample.info.hp.ipf) == colnames(exprs.GPL14550))

# Create design matrix
condition <- factor(sample.info.hp.ipf$condition, levels = c("IPF", "HP"))
design <- model.matrix(~ condition)

# Fit linear model
fit.hp.ipf <- lmFit(exprs.GPL14550, design = design)
fit.hp.ipf <- eBayes(fit.hp.ipf, trend = TRUE, robust = TRUE)
summary(decideTests(fit.hp.ipf[,-1], p.value = 0.05, lfc = 1))

res.hp.ipf <- topTable(fit.hp.ipf, coef = "conditionHP", adjust.method = "BH", n = Inf)

res.hp.ipf %>%
  rownames_to_column(var = "Gene_Symbol") %>%
  write.csv(., "topTable result_GPL14550_HP_vs_IPF.csv", row.names = FALSE)

degs.hp.ipf <- res.hp.ipf %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
  rownames_to_column(var = "Gene_Symbol")

write.csv(degs.hp.ipf, "DEGs_HP_vs_IPF_GPL14550.csv", row.names = FALSE)


######################################################
############### WORKFLOW 3 - GPL6480 #################
######################################################

# Create design matrix
condition.GPL6480 <- factor(sample.info.GPL6480$condition, levels = c("Control", "IPF", "HP"))
design.GPL6480 <- model.matrix(~ condition.GPL6480)

# Fit linear model
fit.GPL6480 <- lmFit(exprs.GPL6480, design = design.GPL6480)
fitted.ebayes.GPL6480 <- eBayes(fit.GPL6480, trend = TRUE, robust = TRUE)
summary(decideTests(fitted.ebayes.GPL6480[,-1], p.value = 0.05, lfc = 1))



#### 1) HP vs Control
res.HP.GPL6480 <- topTable(fitted.ebayes.GPL6480, coef = "condition.GPL6480HP", 
                            adjust.method = "BH", n = Inf)
degs.HP.GPL6480 <- res.HP.GPL6480 %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
  rownames_to_column(var = "Gene_Symbol")

write.csv(degs.HP.GPL6480, "DEGs_HP_vs_Ctrl_GPL6480.csv", row.names = FALSE)


#### 2) IPF vs Control
res.IPF.GPL6480 <- topTable(fitted.ebayes.GPL6480, coef = "condition.GPL6480IPF", 
                            adjust.method = "BH", n = Inf)
degs.IPF.GPL6480 <- res.IPF.GPL6480 %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
  rownames_to_column(var = "Gene_Symbol")

write.csv(degs.IPF.GPL6480, "DEGs_IPF_vs_Ctrl_GPL6480.csv", row.names = FALSE)



#### 3) ILD vs Control
sample.info.GPL6480 <- sample.info.GPL6480 %>%
  mutate(disease = ifelse(condition == "Control", "Control", "ILD"))

# Check if order is correct
all(rownames(sample.info.GPL6480) == colnames(exprs.GPL6480))

# Create design matrix
condition <- factor(sample.info.GPL6480$disease, levels = c("Control", "ILD"))
design <- model.matrix(~ condition)

# Fit linear model
fit.ild.GPL6480 <- lmFit(exprs.GPL6480, design = design)
fit.ild.GPL6480 <- eBayes(fit.ild.GPL6480, trend = TRUE, robust = TRUE)
summary(decideTests(fit.ild.GPL6480[,-1], p.value = 0.05, lfc = 1))

res.ild.GPL6480 <- topTable(fit.ild.GPL6480, coef = "conditionILD", adjust.method = "BH", n = Inf)
degs.ild.GPL6480 <- res.ild.GPL6480 %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
  rownames_to_column(var = "Gene_Symbol")

write.csv(degs.ild.GPL6480, "DEGs_ILD_vs_Ctrl_GPL6480.csv", row.names = FALSE)




#### 4) HP vs IPF
sample.info.hp.ipf <- sample.info.GPL6480 %>%
  filter(condition == "HP" | condition == "IPF") %>%
  select(geo_accession, condition) 

# Check if order is correct
exprs.GPL6480 <- exprs.GPL6480[, colnames(exprs.GPL6480) %in% rownames(sample.info.hp.ipf)]
all(rownames(sample.info.hp.ipf) == colnames(exprs.GPL6480))

# Create design matrix
condition <- factor(sample.info.hp.ipf$condition, levels = c("IPF", "HP"))
design <- model.matrix(~ condition)

# Fit linear model
fit.hp.ipf <- lmFit(exprs.GPL6480, design = design)
fit.hp.ipf <- eBayes(fit.hp.ipf, trend = TRUE, robust = TRUE)
summary(decideTests(fit.hp.ipf[,-1], p.value = 0.05, lfc = 1))

#OUTPUT - conditionHP
#Down             0
#NotSig       17525
#Up               0

# No further analysis done with this as no DEGs
