library(dplyr)
library(limma)
library(ggplot2)
library(EnhancedVolcano)

###############################################################################
# Read Data
## set working directory
dir <- "/Users/mamourie/Library/CloudStorage/Box-Box/"
dir_ampAD <- "/Users/mamourie/Library/CloudStorage/Box-Box/Elizabeth-Project-Data/AMP-AD/DiversityCohort/"
# dir <- "/Users/elizabethmamourian/Library/CloudStorage/Box-Box/"
# dir_ampAD <- "/Users/elizabethmamourian/Library/CloudStorage/Box-Box/Elizabeth-Project-Data/AMP-AD/DiversityCohort/"

## load QC'ed proetomics data
prot.orig <- read.csv(paste0(dir_ampAD, "syn59611693/QC/normAbundances_post-FP_frontal.csv"), stringsAsFactors = FALSE)

## load clinical data w/ Subject ID
clin.orig <- read.csv(paste0(dir_ampAD, "syn51732482/Data/Metadata/AMP-AD_DiverseCohorts_individual_metadata.csv"), stringsAsFactors = FALSE)

## load mapping between sample ID & Subject ID
id_map <- read.csv(paste0(dir_ampAD, "syn51732482/Data/Metadata/AMP-AD_DiverseCohorts_biospecimen_metadata.csv"), stringsAsFactors = FALSE)

################################################################################
# Prepare data format, primary keys; review demographics
## create subject metadata file
clin_id <- merge(clin.orig, id_map, by="individualID") %>%
  dplyr::filter(specimenID %in% names(prot.orig)) # selects sample IDs with available proteomics data
table(clin_id$ADoutcome, exclude=NULL)
clin_id[clin_id$ageDeath=="90+", "ageDeath"] <- 90

### determine demographics of each diagnosis group ####
clin_ad <- clin_id %>%
  dplyr::filter(ADoutcome== "AD")
clin_ctrl <- clin_id %>%
  dplyr::filter(ADoutcome== "Control")
clin_other <- clin_id %>%
  dplyr::filter(ADoutcome== "Other")
clin_missing <- clin_id %>%
  dplyr::filter(ADoutcome== "missing or unknown")
# sex
table(clin_ad$sex) 
table(clin_ctrl$sex) 
table(clin_other$sex) 
table(clin_missing$sex) 
# age
mean(as.numeric(clin_ad$ageDeath), na.rm = TRUE)
sd(as.numeric(clin_ad$ageDeath), na.rm = TRUE)
mean(as.numeric(clin_ctrl$ageDeath))
sd(as.numeric(clin_ctrl$ageDeath))
mean(as.numeric(clin_other$ageDeath))
sd(as.numeric(clin_other$ageDeath))
mean(as.numeric(clin_missing$ageDeath))
sd(as.numeric(clin_missing$ageDeath))
# remove temp data frames
rm(clin_ad, clin_ctrl, clin_other, clin_missing)
### end demographics ####

## match genes list to SDoH genes of interest
prot_ids <- as.data.frame(stringr::str_extract(prot.orig$X, "[^|]+")) %>%
  dplyr::rename(Gene= "stringr::str_extract(prot.orig$X, \"[^|]+\")")
prot_ids$X <- prot.orig$X

## select samples with TRIM11 measurement
prot_sel <- as.data.frame(t(prot.orig))
names(prot_sel) <- prot_sel["X",] # name columns with gene names
prot_sel_named <- prot_sel[-1,] # remove gene names from data contents
prot_sel_named <- prot_sel_named %>%
  mutate_all(~ as.numeric(as.character(.))) # make data numeric
TRIM11_df <- prot_sel_named %>%
  dplyr::filter(!is.na(`TRIM11|Q96F44`)) # selects 228 observations of TRIM11 protein

## create order of subjects/samples to apply to metadata
prot_order <- as.data.frame(row.names(TRIM11_df))
names(prot_order) <- "specimenID"
prot_order$order <- 1:length(prot_order$specimenID)

## apply order of samples in "TRIM11_df" to subject metadata
clin_order <- merge(clin_id, prot_order, all.y=TRUE) %>%
  dplyr::arrange(order)

## check order of metadata against order of proteomics data
clin_order$specimenID[1:10]
row.names(TRIM11_df)[1:10]

### try with & without log normalization
## Add pseudo count to avoid log(0) 
# TRIM11_df_orig <- TRIM11_df
# TRIM11_df <- log2(TRIM11_df_orig + 1)  

# summarize distribution of TRIM11 observations
summary(TRIM11_df$`TRIM11|Q96F44`)
hist(TRIM11_df$`TRIM11|Q96F44`)


###############################################################################
# Use linear regression to find association of other measured proteins with TRIM11
## prepare response, predictors, and design matrix
response <- TRIM11_df[["TRIM11|Q96F44"]] # extract response
predictors <- TRIM11_df[, setdiff(colnames(TRIM11_df), "TRIM11|Q96F44")] # prepare predictors (exclude the response variable)
predictors_t <- t(predictors) # transpose predictors: features as rows, samples as columns
design <- model.matrix(~ response) # write design matrix

## Fit linear model using limma
fit <- lmFit(predictors_t, design)
fit <- eBayes(fit)

## Get all associations
results <- topTable(fit, coef="response", number=Inf)

## Filter significant features
sig_results <- results[results$adj.P.Val < 0.05, ]

## Write output
write.csv(results, paste0(dir,"mamourie/ShenLab/TRIM11/TRIM11_differentialAbundance_AMP-AD.csv"), row.names = TRUE)
write.csv(sig_results, paste0(dir,"mamourie/ShenLab/TRIM11/TRIM11_differentialAbundance_AMP-AD_significant.csv"), row.names = TRUE)

###############################################################################

table(clin_order$sex)
table(clin_order$ADoutcome)

###############################################################################