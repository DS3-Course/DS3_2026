############################################################
# version: 3.1 # Do not edit except for new version
# title: DESeq2 differential expression analysis # Do not edit except for new version
# project: HTP example for DS3
# author:
#   - Matthew Galbraith
# email:
#   - matthew.galbraith@cuanschutz.edu
# affiliation:
#   - Linda Crnic Institute and Department of Pharmacology
#   - University of Colorado Anschutz Medical Campus
############################################################
# Script original author: Matthew Galbraith                
# version: 3.0  Date: 01_13_2024                           
############################################################
# Change Log:
# V3.0 010524
# NEW MAJOR VERSION
# V3.1.DS3 051424
# Cleaning up and streamlining workflow for DS3
#


### Summary:  
# Testing for differential expression in HTP PAXgene Whole Blood RNA samples.
# Paired-end PolyA+ Globin-depleted libraries generated and sequenced by Novagene.
# See PMID 37379383.
# 

### Data type(s):
# A. Bulk RNAseq counts
#      Where/who did this data come from?
#      What is the source of the original data and where is it stored?
# B. Bulk RNAseq RPKMs data
#      Where/who did this data come from?
#      What is the source of the original data and where is it stored?
# C. HTP meta data
#      Where/who did this data come from?
#      What is the source of the original data and where is it stored?
#

### Workflow:
# 1.  Read in meta, counts, rpkms data for all samples and QC  
# 2.  Filter by minimum cpm  
# 3.  Comparison groups and covariates setup  
# 4.  Generate DESeqDataSet(s) and assess models
# 5.  Run DESeq2 analysis  
# 6.  QC checks and overall sample groupings
# 7.  Get results for comparison(s) of interest  
# 8.  Export results
# 9.  Extra plots
# 10. Volcano plot(s)
# 11. Manhattan plot(s)
# 12. Individual gene plots
# 13. Heatmap(s)
# 14. GSEA
#

### Comments:
# Any further relevant details?  
#   


# 0 General Setup -----
# RUN THIS FIRST TIME - Initialize and install packages with renv:
# renv::init(bioconductor = TRUE)
## 0.1 Load required libraries ----
library("DESeq2") # differential expression analysis
library("edgeR") # for cpm() function (can also be used for differential expression analysis)
library("limma") # for removeBatchEffect() function (can also be used for differential expression analysis)
library("BiocParallel") # enables mutli-cpu for some of DEseq2 functions
ncores <- parallel::detectCores() - 1
register(MulticoreParam(workers = ncores)) # enables mutli-cpu for some of DEseq2 functions; can set number of workers - default is all cores; 
library("apeglm") # used with DESeq2 to 'moderate' fold-changes
# library("sva") # for Surrogate Vector Analysis (optional)
library("readxl") # reading Excel files
library("openxlsx") # for exporting results as Excel workbooks
#
library("gplots") # required for current sample-sample distance heatmap
library("genefilter") # used for function rowVars?
library("graphics") # used for dendrograms?
library("dendextend") # used for coloring dendograms
#
library("tidyverse") # required for ggplot2, dplyr etc
library("ggforce") # used for sina plots
library("ggrastr") # required for rasterizing some layers of plots
library("ggrepel") # required for using geom_text and geom_text_repel() to make sample labels for PCA plot
library("RColorBrewer") # color palettes
library("circlize") # color scale generation
library("tidyHeatmap") # tidy heatmaps
# library("factoextra") # extraction and visualization for PCA
library("conflicted") # force all conflicts to become errors
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("count", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("paste", "base")
conflict_prefer("rowVars", "matrixStats")
library("here")
# detach("package:here", unload=TRUE) # run this to reset here()
source(here("helper_functions_DESeq.R")) # load helper functions
#

## 0.2 renv setup ----
# see vignette("renv")
# The general workflow when working with renv is:
#  1 Call renv::init() to initialize a new project-local environment with a private R library,
#  2 Work in the project as normal, installing and removing new R packages as they are needed in the project,
#    recommend using renv::install() as this will automatically update lockfile
#  3 Call renv::snapshot() to save the state of the project library to the lockfile (called renv.lock),
#  4 Continue working on your project, installing and updating R packages as needed.
#  5 Call renv::snapshot() again to save the state of your project library if
# your attempts to update R packages were successful, or call renv::restore() to
# revert to the previous state as encoded in the lockfile if your attempts to
# update packages introduced some new problems.
# 
# IMPORTANT
# To initialize renv in a project using Bioconductor
# # use the latest-available Bioconductor release
# renv::init(bioconductor = TRUE)
# to use a specific version of Bioconductor
# renv::init(bioconductor = "3.14")
#

## 0.3 Set required parameters ----
# Input data files
gene_anno_file <- here("data", "gencode.v33.basic.annotation_gtf_anno.txt.gz") # HUMAN
counts_file <- here("data", "HTP_WholeBlood_RNAseq_GencodeV33_Counts.txt.gz")
rpkms_file <- here("data", "HTP_WholeBlood_RNAseq_GencodeV33_RPKMs.txt.gz")
htp_meta_data_file <- here("data", "HTP_Metadata_v0.5_RNAseq.txt.gz")
#
# Settings
min_cpm <- 0.5 # used for low count filtering; default is 0.5
min_samples <- "auto" # used for low count filtering; use a number, "all", or "auto" (sets to half number of samples)
standard_colors <- c("Control" = "grey30", "T21" = "#009b4e") # these should be named
#
out_file_prefix <- "HTP_RNAseq_DESeq2_v3.1.DS3_" # should match this script title
# End required parameters ###


# 1. Read in and inspect data ----
#
## 1.0 Load gene names and other annotation ----
gene_anno <- gene_anno_file %>%
  read_tsv()
gene_anno
#


## 1.1 Read in counts and rpkms data ----
counts_data <- counts_file |> 
  read_tsv() |> 
  rename(Geneid = EnsemblID)
counts_data
#
rpkms_data <- rpkms_file |> 
  read_tsv() |> 
  rename(Geneid = EnsemblID)
rpkms_data
#


## 1.2. Read in meta data ----
#
meta_data <- htp_meta_data_file |> 
  read_tsv() |> 
  mutate(
    Karyotype = fct_relevel(Karyotype, c("Control", "T21")), # convert to factor and set order
    Sex = fct_relevel(Sex, "Female"), # convert to factor and set order
    Sample_source = as_factor(Sample_source_code), # convert to factor - default is numerical order
    Sample_source_code = NULL # drop
  ) |> 
  rename(Sampleid = LabID)
######## CAN SUBSET HERE TO SAVE TIME ########
# meta_data <- meta_data |> 
#   group_by(Karyotype, Sex) |> 
#   sample_n(10) |> 
#   ungroup()
######## 
# inspect
meta_data
meta_data |> skimr::skim()
#
here("data/HTP_Metadata_v0.5_dictionary.txt") |> read_tsv()
#


## 1.3 Convert unfiltered counts to tidy/long format and join with meta data ----
# (one row per observation, one column per variable)
#
# Subselect counts data and rpkms to match meta_data
# (ensures same samples and order as in meta_data)
counts_data <- counts_data %>% select(Geneid, meta_data %>% pull(Sampleid))
rpkms_data <- rpkms_data %>% select(Geneid, meta_data %>% pull(Sampleid))
#
# Gather to long format
counts_data_long <- counts_data %>% 
  pivot_longer(-Geneid, names_to = "Sampleid", values_to = "raw_count")
counts_data_long
#
rpkms_data_long <- rpkms_data %>%
  pivot_longer(-Geneid, names_to = "Sampleid", values_to = "RPKM")
rpkms_data_long
#
# Join with metadata
counts_data_long <- counts_data_long %>% 
  inner_join(meta_data)
# Check number of samples before and after join
n_samples <- counts_data %>% select(-Geneid) %>% colnames() %>% length()
n_samples_join <- counts_data_long %>% distinct(Sampleid) %>% nrow()
cat("Number of samples in counts_data:", n_samples,
    "\nNumber of samples after join with meta data:", n_samples_join)
if(n_samples != n_samples_join) warning("Possible problem joining counts_data_long with meta data - different number of samples after join\n (ignore if subsetting by Karyotype and/or comorbidities etc)")
#

## 1.4 Check raw read count distributions across samples (not normalized) ----
counts_data_long %>%
  ######## ONLY USING A SUBSET HERE TO SAVE TIME ########
  filter(Sampleid %in% (meta_data |> group_by(Karyotype) |> sample_n(5) |> pull(Sampleid))) |> 
  mutate(Sampleid = fct_relevel(Sampleid, meta_data |> arrange(Karyotype) %>% pull(Sampleid) %>% as.character)) %>% 
  ggplot(aes(Sampleid, log2(raw_count + 0.1), color = Karyotype)) + # CUSTOMIZE AS NEEDED
  geom_sina(size = 0.01) +
  scale_color_manual(values = standard_colors) +
  labs(title = "Raw read count distributions across samples") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_counts_unfiltered.png")), width = 10, height = 5, units = "in")
#

## 1.5 Check globin genes (Optional) -----
# Given the presence of RBCs in the PAX RNA preps, what is the proportion of reads taken up by Globin mRNAs?
# Up to 70% of the mRNA (by mass) in whole blood total RNA are globin transcripts
# Globin depletion was carried out using:
#
# 'Expecting' ~95% depeletion: .05 x .7 = ≤ 3.5% of reads
gene_percents <- counts_data_long %>%
  ######## ONLY USING A SUBSET HERE TO SAVE TIME ########
  filter(Sampleid %in% (meta_data |> group_by(Karyotype) |> sample_n(5) |> pull(Sampleid))) |> 
  group_by(Sampleid) %>%
  mutate(total = sum(raw_count)) %>%
  group_by(Sampleid, Geneid) %>%
  mutate(percent_total = raw_count / total * 100) %>%
  group_by(Sampleid) %>%
  mutate(rank = min_rank(-percent_total)) %>%
  ungroup() %>%
  distinct(Geneid, Sampleid, total, percent_total, rank)
globins <- gene_anno %>%
  filter(str_detect(gene_name, "^HBA1$|^HBA2$|^HBB$|^HBG1$|^HBG2$|^HBD$"))
gene_percents |>
  inner_join(globins) %>%
  inner_join(meta_data) |>
  mutate(Sampleid = fct_relevel(Sampleid, (gene_percents |> filter(Geneid == "ENSG00000244734.4") |> arrange(-percent_total))$Sampleid)) |>
  ggplot(aes(Sampleid, percent_total, fill = gene_name)) +
  geom_col() +
  facet_wrap(~ Karyotype, scale = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Globin mRNAs - percentage of total reads")
#


# 2. Filter/remove genes with low expression ----
#
## 2.1 Summary of total reads per sample ----
counts_data %>% 
  select(-Geneid) %>% 
  colSums() %>% 
  summary()
counts_data %>% 
  select(-Geneid) %>% 
  colSums() %>% 
  enframe(name = "Sampleid", value = "Total_reads") |> 
  arrange(Total_reads)
#

## 2.2 Filter by minimum counts per million ----
# Keep only rows (transcripts / genes) with greater than `min_cpm` cpm in `min_samples`
# NOTE: 10 counts=0.5 cpm for 20 million reads,  15 counts=0.5 cpm for 30 million reads...
# from Michael Love: https://support.bioconductor.org/p/95840/
# The independent filtering is designed only to filter out low count genes to
# the extent that they are not enriched with small p-values. Here the problem is
# not independent filtering, but that these two genes get a small p-value rather
# than being filtered or having an insignificant p-value. Datasets can be
# different in many ways, and for whatever reason, these two genes survive the
# filtering and get a counterintuitive small p-value. I'd recommend you just use
# a more strict filter in the very beginning, e.g. at least three samples with
# counts greater than 10
#
# Check min_samples and calculate if needed (ie if number is not supplied)
if (min_samples == "all") {
  min_samples=ncol(counts_data) - 1
} else if (min_samples == "auto") {
  min_samples=(ncol(counts_data) - 1) / 2
}
before <- counts_data %>% 
  transmute(
    Geneid = Geneid, 
    row_sum = rowSums(select(., -Geneid))
  ) %>% 
  filter(row_sum > 0) %>% 
  nrow()
cpm_data <- counts_data %>% 
  column_to_rownames("Geneid") %>% 
  cpm()
keep <- cpm_data |> 
  as_tibble(rownames = "Geneid") |> 
  pivot_longer(-Geneid, names_to = "Sampleid", values_to = "cpm") |> 
  mutate(cpm > min_cpm) |> # check against min_cpm
  filter(`cpm > min_cpm` == TRUE) |> # and filter
  dplyr::count(Geneid) |> # count samples remaining per Geneid
  filter(n >= min_samples) # filter against min_samples
counts_filtered <- counts_data |> 
  filter(Geneid %in% keep$Geneid)

## 2.3 Summarize rows before and after filtering ----
cat("Total number of rows: ", counts_data %>% nrow())
cat("Number of rows with non-zero total read counts before filtering: ", before)
cat("Number of rows after filtering: ", counts_filtered %>% nrow())
# Gather to long format and join with meta data
counts_filtered_long <- counts_filtered %>% 
  pivot_longer(-Geneid, names_to = "Sampleid", values_to = "raw_count") |> 
  inner_join(meta_data)
# 

## 2.4 Check filtered read count distributions for each sample (not normalized) ----
counts_filtered_long %>%
  ######## ONLY PLOTTING A SUBSET HERE TO SAVE TIME ########
  filter(Sampleid %in% (meta_data |> group_by(Karyotype) |> sample_n(5) |> pull(Sampleid))) |> 
  mutate(Sampleid = fct_relevel(Sampleid, meta_data |> arrange(Karyotype) %>% pull(Sampleid) %>% as.character)) %>% 
  ggplot(aes(Sampleid, log2(raw_count + 0.1), color = Karyotype)) + # CUSTOMIZE AS NEEDED
  geom_sina(size = 0.01) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.2, color = "black", fill = "transparent", size = 0.75) +
  scale_color_manual(values = standard_colors) +
  labs(title = "CPM-filtered read count distributions across samples", x = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_counts_filtered.png")), width = 10, height = 5, units = "in")
#


# 3. Groups and/or Covariates setup ----
groups <- meta_data %>% 
  # CUSTOMIZE
  select(Sampleid, Karyotype, Sex, Age, Sample_source)
  # mutate(
  #   Age = (Age - mean(Age)) / sd(Age), # convert to Z-scores to scale and center
  # )
#


# 4. Generate DESeqDataSet object(s) ----
# Creates *DESeqDataSet* object and populates with count data and experimental design  
#
## 4.1 simple model ----
simple_formula <- as.formula(paste0("~", "Karyotype")) # REQUIRES CUSTOMIZATION
#
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered %>% 
    select(Geneid, groups %>% pull(Sampleid)) %>%  # ensures correct order of columns
    column_to_rownames("Geneid"), # must be converted to data frame from tibble
  colData = groups,
  design = simple_formula
)
#
# Check meta data read in to DESeqDataSet (simple version):
colData(dds)
#

## 4.2 multivariable model(s) ----
# Put the variable of interest at the end of the formula. Thus the results
# function will by default pull the condition results unless contrast or name
# arguments are specified.
# see http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#
multivar_formula <- as.formula(paste0("~", "Sex +", "Karyotype")) # REQUIRES CUSTOMIZATION
multivar_formula2 <- as.formula(paste0("~ ", "Age +", "Karyotype")) # REQUIRES CUSTOMIZATION
multivar_formula3 <- as.formula(paste0("~ ", "Sample_source +", "Karyotype")) # REQUIRES CUSTOMIZATION
multivar_formula4 <- as.formula(paste0("~ ", "Sex + Age + Sample_source +", "Karyotype")) # REQUIRES CUSTOMIZATION
#
dds_multi <- DESeqDataSetFromMatrix(
  countData = counts_filtered %>% 
    select(Geneid, groups %>% pull(Sampleid)) %>%  # ensures correct order of columns
    column_to_rownames("Geneid"), # must be converted to data frame from tibble
  colData = groups,
  design = multivar_formula
)
# Check meta data read in to DESeqDataSet (multivariable version):
colData(dds_multi)
#
dds_multi2 <- DESeqDataSetFromMatrix(
  countData = counts_filtered %>%
    select(Geneid, groups %>% pull(Sampleid)) %>%  # ensures correct order of columns
    column_to_rownames("Geneid"), # must be converted to data frame from tibble
  colData = groups,
  design = multivar_formula2
)
# 
dds_multi3 <- DESeqDataSetFromMatrix(
  countData = counts_filtered %>%
    select(Geneid, groups %>% pull(Sampleid)) %>%  # ensures correct order of columns
    column_to_rownames("Geneid"), # must be converted to data frame from tibble
  colData = groups,
  design = multivar_formula3
)
# 
dds_multi4 <- DESeqDataSetFromMatrix(
  countData = counts_filtered %>%
    select(Geneid, groups %>% pull(Sampleid)) %>%  # ensures correct order of columns
    column_to_rownames("Geneid"), # must be converted to data frame from tibble
  colData = groups,
  design = multivar_formula4
)
# 
## 4.3 Likelihood ratio test for multivariable model(s) ----
# For how many genes does adding additional terms give `better` fit?
dds_multi_lrt <- DESeq(dds_multi, parallel = TRUE, test = "LRT", reduced = simple_formula)
dds_multi_lrt2 <- DESeq(dds_multi2, parallel = TRUE, test = "LRT", reduced = simple_formula)
dds_multi_lrt3 <- DESeq(dds_multi3, parallel = TRUE, test = "LRT", reduced = simple_formula)
dds_multi_lrt4 <- DESeq(dds_multi4, parallel = TRUE, test = "LRT", reduced = simple_formula)
#
dds_multi_lrt %>% results() %>% elementMetadata() %>% as_tibble() %>% filter(str_detect(description, "LRT p-value")) %>% pull()
dds_multi_lrt %>% results(.) %>% as_tibble(rownames = "Geneid") %>% count(padj < 0.1)
# 
dds_multi_lrt2 %>% results() %>% elementMetadata() %>% as_tibble() %>% filter(str_detect(description, "LRT p-value")) %>% pull()
dds_multi_lrt2 %>% results(.) %>% as_tibble(rownames = "Geneid") %>% count(padj < 0.1)
#
dds_multi_lrt3 %>% results() %>% elementMetadata() %>% as_tibble() %>% filter(str_detect(description, "LRT p-value")) %>% pull()
dds_multi_lrt3 %>% results(.) %>% as_tibble(rownames = "Geneid") %>% count(padj < 0.1)
#
dds_multi_lrt4 %>% results() %>% elementMetadata() %>% as_tibble() %>% filter(str_detect(description, "LRT p-value")) %>% pull()
dds_multi_lrt4 %>% results(.) %>% as_tibble(rownames = "Geneid") %>% count(padj < 0.1)
#


## 4.4 SVA model (optional) -------
# The goal of the sva is to remove all unwanted sources of variation while
# protecting the contrasts due to the primary variables included in mod. This
# leads to the identification of features that are consistently different
# between groups, removing all common sources of latent variation.
# dat  <- counts(DESeq(dds, parallel=TRUE), normalized = TRUE)
# # best to ensure genes with low counts have already been removed
# mod  <- model.matrix(~ Group, colData(dds)) # designate variable of interest to "protect"
# mod0 <- model.matrix(~ 1, colData(dds))
# #
# svseq <- svaseq(dat, mod, mod0, n.sv = NULL) # if NULL, number of factors will be estimated for you 
# #
# groups_sva <- svseq %>% 
#   # extract + add additional coefficients to the estimated surrogate variables to allow comparison
#   biobroom::tidy_sva(addVar = colData(DESeq(dds, parallel=TRUE))) # need to run DESeq to get sizeFactors
# #

### Compare SVs with known categorical variables ----
# groups_sva %>%
#   pivot_longer(contains("sv"), names_to = "SV", values_to = "estimate") %>%
#   ggplot(aes(Sex, estimate)) +
#   geom_sina() +
#   facet_wrap(~SV)
# groups_sva %>% 
#   pivot_longer(contains("sv"), names_to = "SV", values_to = "estimate") %>% 
#   ggplot(aes(Group, estimate)) + # Group should be preserved ie not correlated to any SV
#   geom_sina() +
#   facet_wrap(~SV)

### Compare SVs with known continuous variables
# groups_sva %>% 
#   pivot_longer(contains("sv"), names_to = "SV", values_to = "estimate") %>% 
#   ggplot(aes(sizeFactor, estimate)) +
#   geom_point() +
#   facet_wrap(~SV)
# groups_sva %>% 
#   pivot_longer(contains("sv"), names_to = "SV", values_to = "estimate") %>% 
#   ggplot(aes(Age, estimate)) +
#   geom_point() +
#   facet_wrap(~SV)
# groups_sva %>%
#   pivot_longer(contains("sv"), names_to = "SV", values_to = "estimate") %>%
#   ggplot(aes(RIN, estimate)) +
#   geom_point() +
#   facet_wrap(~SV)

### Build DESeq model with SVs -----
# sva_formula <- as.formula(
#   paste0("~", groups_sva %>% select(matches("^sv")) %>% colnames() %>% paste(collapse = " + "),
#          "+ Group" # CUSTOMIZE
#   ))
# #
# dds_sva <- DESeqDataSetFromMatrix(
#   countData = counts_filtered %>% 
#     select(Geneid, groups_sva %>% pull(Sampleid)) %>%  # ensures correct order of columns
#     column_to_rownames("Geneid"), # must be converted to data frame from tibble
#   colData = groups_sva,
#   design = sva_formula
# )
# # Check metadata read in to DESeqDataSet (sva version):
# colData(dds_sva)
# # 

### Likelihood ratio test for SVA model(s) ----
# # For how many genes does adding surrogate vectors give `better` fit?
# dds_sva_lrt <- DESeq(dds_sva, parallel = TRUE, test = "LRT", reduced = simple_formula)
# dds_sva_lrt %>% results() %>% elementMetadata() %>% as_tibble() %>% filter(str_detect(description, "LRT p-value")) %>% pull()
# dds_sva_lrt %>% results(.) %>% as_tibble(rownames = "Geneid") %>% count(padj < 0.1) # in this case = 10376
# #


# 5. Run DESeq2 analysis ----
# Default analysis runs the following steps:  
# 1. estimation of size factors  
# 2. estimation of dispersion  
# 3. Negative Binomial GLM fitting and Wald statistics  
# note: since ~v1.16, shrinkage of log2foldChange is not run by default - this has moved to lfcShrink function
#
## 5.1 Run simple model ----
dds <- DESeq(dds, parallel=TRUE)
#
## 5.2 Run multivariable model(s) ----
dds_multi4 <- DESeq(dds_multi4, parallel=TRUE)
#
## 5.3 Run SVA model(s) (optional) ----
# dds_sva <- DESeq(dds_sva, parallel=TRUE)
#



# 6. QC checks and overall sample grouping(s) ----  
#
## 6.1 Check Size Factors used for normalization ----
get_size_fcts(dds) |> 
  arrange(SizeFactor) |> 
  summary()
#
## 6.2 Get normalized counts ----
nc <- dds %>% counts(normalized = TRUE) %>% as_tibble(rownames="Geneid")
nc
#
### Generate covariate-adjusted normalized counts ----
nc_SexAgeSource_adj <- nc %>%
  mutate_at(2:ncol(.), ~ log2(.)) %>%  # log2 transformation of counts
  column_to_rownames(var = "Geneid") %>% # convert to dataframe to preserve Geneid
  limma::removeBatchEffect(
    batch = colData(dds_multi) %>% as_tibble() %>% pull(Sex), # REQUIRES CUSOMIZATION; DO NOT INCLUDE PREDICTOR OF INTEREST
    batch2 = colData(dds_multi) %>% as_tibble() %>% pull(Sample_source), # REQUIRES CUSOMIZATION; DO NOT INCLUDE PREDICTOR OF INTEREST
    covariates = colData(dds_multi) %>% as_tibble() %>% select(Sampleid, Age) |> column_to_rownames(var = "Sampleid"),
    design = colData(dds_multi) %>% as_tibble() %>% model.matrix(simple_formula, data = .) # INCLUDE ONLY PREDICTOR OF INTEREST
  ) %>%
  as_tibble(rownames="Geneid") %>% # convert back to tibble
  mutate_at(2:ncol(.), ~ 2^(.)) # remove log2 transformation
nc_SexAgeSource_adj
#
### Generate sv-adjusted normalized counts (optional) ----
# nc_sva_adj <- nc %>%
#   mutate_at(2:ncol(.), ~ log2(.)) %>%  # log2 transformation of counts
#   column_to_rownames(var = "Geneid") %>% # convert to dataframe to preserve Geneid
#   limma::removeBatchEffect(
#     covariates = colData(dds_sva) %>% as_tibble() %>% select(matches("^sv")),
#     design = colData(dds_sva) %>% as_tibble() %>% model.matrix(simple_formula, data = .) # INCLUDE ONLY PREDICTOR OF INTEREST
#   ) %>%
#   as_tibble(rownames="Geneid") %>% # convert back to tibble
#   mutate_at(2:ncol(.), ~ 2^(.)) # remove log2 transformation
# #


## 6.3 Get vst transformed values ----
# Variance-stablizing transformation and normalization ± covariate correction  
vst_mat <- assay(vst(dds))
vst_mat |> as_tibble(rownames = "Geneid")
#
### Generate covariate-adjusted vst values ----
vst_mat_SexAgeSource_adj <-
  vst_mat %>%
  limma::removeBatchEffect(
    batch = colData(dds_multi4) %>% as_tibble() %>% pull(Sex), # REQUIRES CUSOMIZATION; DO NOT INCLUDE PREDICTOR OF INTEREST
    batch2 = colData(dds_multi4) %>% as_tibble() %>% pull(Sample_source), # REQUIRES CUSOMIZATION; DO NOT INCLUDE PREDICTOR OF INTEREST
    covariates = colData(dds_multi4) %>% as_tibble() %>% select(Sampleid, Age) |> column_to_rownames(var = "Sampleid"),
    design = colData(dds_multi4) %>% as_tibble() %>% model.matrix(simple_formula, data = .) # INCLUDE ONLY PREDICTOR OF INTEREST
  )
vst_mat_SexAgeSource_adj |> as_tibble(rownames = "Geneid")
#
### Generate sv-adjusted vst values (optional) ----
# vst_mat_sva_adj <-
#   vst_mat %>%
#   limma::removeBatchEffect(
#     covariates = colData(dds_sva) %>% as_tibble() %>% select(matches("^sv")),
#     design = colData(dds_sva) %>% as_tibble() %>% model.matrix(simple_formula, data = .) # INCLUDE ONLY PREDICTOR OF INTEREST
#   )
# #


## 6.4 Dendrogram and hierarchical clustering of sample-sample distances -----
#
plotDendClust2(vst_mat, adjustment = "Unadjusted", color_var = "Karyotype")
#
plotDendClust2(vst_mat_SexAgeSource_adj, adjustment = "SexAgeSource_adj", color_var = "Karyotype")
#



## 6.5 PCA plot(s) of normalized and variance-transformed count data ± covariate correction ----
#
### Overall groupings by PCA ----
plotPCA_custom2(vst_mat, 
                PCA_by = "variance", 
                save_PCs = TRUE, # available as "PC_loadings"
                plot_title = "PCA plot", 
                subtitle = "unadjusted", 
                color_var = "Karyotype",
                shapes = "Sex", 
                labels = FALSE, 
                x_lower_lim = -85, # CUSTOMIZE LIMITS
                x_upper_lim = 85, 
                y_lower_lim = -85, 
                y_upper_lim = 85) 
#
plotPCA_custom2(x = vst_mat_SexAgeSource_adj, 
                PCA_by = "variance", 
                save_PCs = TRUE, # available as "PC_loadings"
                plot_title="PCA plot", 
                subtitle = "SexAgeSource-adjusted",
                color_var = "Karyotype",
                shapes = "Sex", 
                labels = FALSE, 
                x_lower_lim = -85, # CUSTOMIZE LIMITS
                x_upper_lim = 85, 
                y_lower_lim = -85, 
                y_upper_lim = 85)
#
### Verify Karyotype groupings by PCA ----
# plotPCA_custom2(x = vst_mat[rownames(vst_mat) %in% (gene_anno %>% filter(chr == "chr16" & start >= 75540514 & end <= 97962622))$Geneid,], # MOUSE
plotPCA_custom2(x = vst_mat[rownames(vst_mat) %in% (gene_anno %>% filter(chr == "chr21"))$Geneid,], # HUMAN
                PCA_by = "variance",
                save_PCs = TRUE, # available as "PC_loadings"
                plot_title="PCA plot - chr21 genes", # HUMAN
                # plot_title="PCA plot - chr16 triplicated genes", # MOUSE
                subtitle = "unadjusted",
                color_var = "Karyotype",
                shapes = "Karyotype", 
                labels = FALSE,
                x_lower_lim = -10,
                x_upper_lim = 10,
                y_lower_lim = -10,
                y_upper_lim = 10)
#
# plotPCA_custom2(x = vst_mat[rownames(vst_mat) %in% (gene_anno %>% filter(chr == "chr16" & start >= 75540514 & end <= 97962622))$Geneid,], # MOUSE
plotPCA_custom2(x = vst_mat_SexAgeSource_adj[rownames(vst_mat_SexAgeSource_adj) %in% (gene_anno %>% filter(chr == "chr21"))$Geneid,], # HUMAN
                PCA_by = "variance",
                save_PCs = TRUE, # available as "PC_loadings"
                plot_title="PCA plot - chr21 genes", # HUMAN
                # plot_title="PCA plot - chr16 triplicated genes", # MOUSE
                subtitle = "SexAgeSource-adjusted",
                color_var = "Karyotype",
                shapes = "Karyotype", 
                labels = FALSE,
                x_lower_lim = -10,
                x_upper_lim = 10,
                y_lower_lim = -10,
                y_upper_lim = 10)
#
### Verify sex groupings by PCA ----
plotPCA_custom2(x = vst_mat[rownames(vst_mat) %in% (gene_anno %>% filter(chr == "chrY" | chr == "chrX"))$Geneid,],
                PCA_by = "variance",
                save_PCs = TRUE, # available as "PC_loadings"
                plot_title="PCA plot - chrY+X",
                subtitle = "unadjusted",
                color_var = "Sex",
                shapes = "Sex",
                labels = FALSE,
                x_lower_lim = -20,
                x_upper_lim = 20,
                y_lower_lim = -20,
                y_upper_lim = 20)
#
### Check meta data associations with PCs -------
c("Sex", "Age", "Sample_source", "Karyotype") %>% 
  paste("loadings ~", .) |> 
  set_names() %>% 
  map(., ~ pca_lm_function(unadjusted_PC_loadings, unadjusted_PC_percVar, .))
#
c("Sex", "Age", "Sample_source", "Karyotype") %>% 
  paste("loadings ~", .) |> 
  set_names() %>% 
  map(., ~ pca_lm_function(`SexAgeSource-adjusted_PC_loadings`, `SexAgeSource-adjusted_PC_percVar`, .))
#


# 7. Get DESeq2 results ---- 
# see http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# To see available results: dds_batch %>% resultsNames() list of results
# available does not show all combinations - only comparisons against reference
# level of last variable in the design formula are available by default
#
# The results function without any arguments will automatically perform a
# contrast of the last level of the last variable in the design formula over the
# first level, based on levels of thisdesign$condition and thus may not be
# useful; design could be releveled to report the main comparison of interest
# but can specify the comparison of interest when calling results()
#
# Usually need to use contrast argument to results() but some types of
# comparison such as continuous variables or interaction terms will require the
# alternative names argument.
#
## NOTE cpm filtering does not make much difference with independent filtering
## on, but will get rid of some odd cases that get a small p-value despite low
## counts and may speed up analysis time
#
#


## 7.1 Define comparisons of interest ----
# CONTRASTS VERSION:
# REQUIRES CUSTOMIZATION (one row for each comparison)
# ensure that levels are in desired order so fold-change is calculated in correct direction
comparisons <- list( 
  c("Karyotype", "T21", "Control") 
)
#

## 7.2 Results summaries ----
#
for (i in comparisons) { # CONTRASTS VERSION
  dds %>% get_results_sum(i, show_ind_filt_off=FALSE)
}
for (i in comparisons) { # CONTRASTS VERSION
  dds_multi4 %>% get_results_sum(i, show_ind_filt_off=FALSE)
}
#

## 7.3 Assemble DESeq2 results table(s) ----
# CONTRASTS VERSION (spaces or dashes in variables will cause problems here)
#
# Initialize empty vector to store names of results objects for later reference
comparisons_results <- character() 
### Simple model ----
for (comparison in comparisons) {
  name <- paste("res_simple", comparison[2], "vs", comparison[3], sep = "_") # CUSTOMIZE name of results table
  comparisons_results <- c(comparisons_results, name)
  res_temp <- dds |> # SET to correct dds object
    get_results_tbl(
      contrast = comparison,
      shrink_type = "apeglm"
    )
  res_temp %>% assign(name, ., pos=1)
}
### Multivariable model(s) ----
for (comparison in comparisons) {
  name <- paste("res_multi_SexAgeSource", comparison[2], "vs", comparison[3], sep="_")  # CUSTOMIZE name of results table
  comparisons_results <- c(comparisons_results, name)
  res_temp <- dds_multi4 |> # SET to correct dds object
    get_results_tbl(
      contrast = comparison,
      shrink_type = "apeglm"
    )
  res_temp %>% assign(name, ., pos=1)
}
### SVA model(s) (Optional) ----
# for (comparison in comparisons) {
#   name <- paste("res_sva", comparison[2], "vs", comparison[3], sep="_") # CUSTOMIZE name of results table
#   comparisons_results <- c(comparisons_results, name)
#   res_temp <- dds_sva |> # SET to correct dds object
#     get_results_tbl(
#       contrast = comparison,
#       shrink_type = "apeglm"
#     )
#   res_temp %>% assign(name, ., pos=1)
# }
#
# List of results tables:
comparisons_results
# Preview results tables:
res_simple_T21_vs_Control
res_multi_SexAgeSource_T21_vs_Control


## 7.4 MA plot(s) ----
#
res_simple_T21_vs_Control |> 
  plotDEgg(
    sig = 0.1,
    title = "T21 vs. Control (simple)", # CUSTOMIZE
    subtitle = c("Model: ", as.character(simple_formula)) %>% paste(collapse = "")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "simple_T21_vs_Control", "_MA.png")), width = 5, height = 5, units = "in")
#
res_multi_SexAgeSource_T21_vs_Control |> 
  plotDEgg(
    sig = 0.1,
    title = "T21 vs. Control (multi)", # CUSTOMIZE
    subtitle = c("Model: ", as.character(multivar_formula4)) %>% paste(collapse = "")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "multi_SexAgeSource_T21_vs_Control", "_MA.png")), width = 5, height = 5, units = "in")
#



## 7.5 Standard volcano plot(s) - no gene labels ----
#
res_simple_T21_vs_Control |> 
  plotVolcano(
    sig = 0.1,
    title = "T21 vs. Control (simple)", # CUSTOMIZE
    subtitle = paste("Model:", paste(simple_formula, collapse = ""))
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "simple_T21_vs_Control", "_Volcano.png")), width = 5, height = 5, units = "in")
#
res_multi_SexAgeSource_T21_vs_Control |> 
  plotVolcano(
    sig = 0.1,
    title = "T21 vs. Control (multi)", # CUSTOMIZE
    subtitle = paste("Model:", paste(multivar_formula4, collapse = ""))
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "multi_SexAgeSource_T21_vs_Control", "_Volcano.png")), width = 5, height = 5, units = "in")
#


# 8. Export results ----
#
comparisons_results |> export_res()
# May want only preferred model results to avoid later confusion


# 9. Extra plots ----
#
## 9.1 Plot dispersion estimates ----
#
dds |> ggplotDispEsts()
dds_multi4 |> ggplotDispEsts()
#
## 9.2 Plot p-value distributons ----
#
res_simple_T21_vs_Control |> plotPvals()
res_multi_SexAgeSource_T21_vs_Control |> plotPvals()
#


# 10. labelled volcano plots ----
#
## 10.1 label top significant genes ----
# Compare volcano plots across various models OR across various comparisons
# with labeling of top differential genes OR selected genes of interest
# using 'patchwork' to assemble multiple plots
#
# get y limits 
max(
  res_simple_T21_vs_Control %>% pull(padj) %>% -log10(.) %>% max(na.rm = TRUE),
  res_multi_SexAgeSource_T21_vs_Control %>% pull(padj) %>% -log10(.) %>% max(na.rm = TRUE)
)
# get x limits (nb expand_limits x does not seem to work if xlim is already set?! but y does?!)
max(
  res_simple_T21_vs_Control %>% pull(log2FoldChange_adj) %>% abs(.) %>% max(na.rm = TRUE) %>% ceiling(),
  res_multi_SexAgeSource_T21_vs_Control %>% pull(log2FoldChange_adj) %>% abs(.) %>% max(na.rm = TRUE) %>% ceiling()
)
#
v1 <- res_simple_T21_vs_Control %>%
  volcano_plot_lab(
    title = "T21 vs. Control (simple)",
    subtitle = paste0(
      paste("Model:", paste(simple_formula, collapse = ""), "\n"), # CUSTOMIZE formula
      "[Down: ", (.) %>% filter(padj < 0.1 & FoldChange_adj <1) %>% nrow(), "; Up: ", (.) %>% filter(padj < 0.1 & FoldChange_adj >1) %>% nrow(), "]"
      ),
    labels = TRUE,
    n_labels = 3,
    raster = TRUE
  ) +
  expand_limits(y = c(0, 220.3075), x = c(-1, 1) * 4) # usually set to same limits across all plots
#
v2 <- res_multi_SexAgeSource_T21_vs_Control %>%
  volcano_plot_lab(
    title = "T21 vs. Control (multi)",
    subtitle = paste0(
      paste("Model:", paste(multivar_formula4, collapse = ""), "\n"), # CUSTOMIZE formula
      "[Down: ", (.) %>% filter(padj < 0.1 & FoldChange_adj <1) %>% nrow(), "; Up: ", (.) %>% filter(padj < 0.1 & FoldChange_adj >1) %>% nrow(), "]"
    ),
    labels = TRUE,
    n_labels = 3,
    raster = TRUE
  ) +
  expand_limits(y = c(0, 220.3075), x = c(-1, 1) * 4) # usually set to same limits across all plots
#
v1 + v2 +
  patchwork::plot_layout(guides = 'collect', nrow = 1)
  # patchwork::plot_annotation( # use for single title/subtitle
  #   title = "T21 vs. Control",
  #   subtitle = ""
  # )
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_models_combined", ".pdf")), device = cairo_pdf, width = 25, height = 5, units = "in")
#

## 10.2 label triplicated genes ----
# May want to highlight Chr21 or triplicated genes in mouse models
res_multi_SexAgeSource_T21_vs_Control %>%
  volcano_plot_chr21(
    title = "T21 vs. Control (multi)",
    subtitle = paste0(
      paste("Model:", paste(multivar_formula4, collapse = ""), "\n"), # CUSTOMIZE formula
      "[Down: ", (.) %>% filter(padj < 0.1 & FoldChange_adj <1) %>% nrow(), "; Up: ", (.) %>% filter(padj < 0.1 & FoldChange_adj >1) %>% nrow(), "]"
    ),
    raster = FALSE
  ) +
  expand_limits(y = c(0, 220.3075), x = c(-1, 1) * 4) # usually set to same limits across all plots
ggsave(filename = here("plots", paste0(out_file_prefix, "multi_SexAgeSource_T21_vs_Control", "_Chr21_Volcano.png")), width = 5, height = 5, units = "in")
#


# 11. Manhattan plot(s) (Optional) ----
#
all_detected_genes <- res_multi_SexAgeSource_T21_vs_Control %>% nrow()
chr21_detected_genes <- res_multi_SexAgeSource_T21_vs_Control %>% filter(chr == "chr21") %>% nrow()
chr21_up_genes <- res_multi_SexAgeSource_T21_vs_Control %>% filter(chr == "chr21") %>% filter(padj < 0.1 & FoldChange_adj > 1) %>% nrow()
chr21_dn_genes <- res_multi_SexAgeSource_T21_vs_Control %>% filter(chr == "chr21") %>% filter(padj < 0.1 & FoldChange_adj < 1) %>% nrow()
all_up_genes <- res_multi_SexAgeSource_T21_vs_Control %>% filter(padj < 0.1 & FoldChange_adj > 1) %>% nrow()
all_dn_genes <- res_multi_SexAgeSource_T21_vs_Control %>% filter(padj < 0.1 & FoldChange_adj < 1) %>% nrow()
cat("Proportion of detected Chr21 genes upregulated: ", chr21_up_genes, "/", chr21_detected_genes, " (", (chr21_up_genes / chr21_detected_genes * 100) %>% round(1), "%)", sep = "")
cat("Proportion of detected Chr21 genes downregulated: ", chr21_dn_genes, "/", chr21_detected_genes, " (", (chr21_dn_genes / chr21_detected_genes * 100) %>% round(1), "%)", sep = "")
cat("Proportion of detected genes upregulated: ", all_up_genes, "/", all_detected_genes, " (", (all_up_genes / all_detected_genes * 100) %>% round(1), "%)", sep = "")
cat("Proportion of detected genes downregulated: ", all_dn_genes, "/", all_detected_genes, " (", (all_dn_genes / all_detected_genes * 100) %>% round(1), "%)", sep = "")
cat("Proportion of DE up genes NOT on Chr21: ", ((all_up_genes - chr21_up_genes) / all_up_genes * 100) %>% round(1), "%", sep = "")
cat("Proportion of DE down genes NOT on Chr21: ", ((all_dn_genes - chr21_dn_genes) / all_dn_genes * 100) %>% round(1), "%", sep = "")
#
m1 <- res_multi_SexAgeSource_T21_vs_Control %>%
  # filter(chr == "chr21") %>%
  mutate(
    color = case_when( # color by significance
      padj < 0.1 ~ "padj < 0.1",
      .default = "n.s."
    )
  )|> 
  ggplot(aes(start, log2(FoldChange_adj), color = color)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  geom_point(data = . %>% filter(color == "n.s."), size = 0.3, ) +
  geom_point(data = . %>% filter(color == "padj < 0.1"), size = 0.3, ) +
  scale_color_manual(values = c("padj < 0.1" = "red", "n.s." = "black")) +
  facet_wrap(~chr, scales = "free_x", nrow = 3) +
  # facet_grid(~chr, scales = "free_x", space = "free_x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.9, 0.15)
  ) +
  labs(
    title = "T21 vs. Control (multi)",
    subtitle = paste("Model:", paste(multivar_formula4, collapse = "")),
    x = "Chromosome position"
  )
m1
ggsave(filename = here("plots", paste0(out_file_prefix, "manhattan", ".png")), width = 8, height = 5, units = "in")
# May want to rasterize before saving as pdf
m1r <- ggrastr::rasterize(m1, layers='Point', dpi = 600, dev = "ragg_png") 
ggsave(m1r, filename = here("plots", paste0(out_file_prefix, "manhattan", ".pdf")), device = cairo_pdf, width = 8, height = 5, units = "in")
# 


# 12. Individual gene plots ----
#
## Generate covariate-adjusted RPKM values ----
# Can plot normalized counts or raw RPKMs but usually prefer RPKMs adjusted with same covariables as preferred model
rpkm_SexAgeSource_adj <- rpkms_data %>%
  select(Geneid, (colData(dds_multi4) %>% as_tibble() %>% pull(Sampleid))) %>%
  mutate_at(2:ncol(.), ~ log2(.)) %>%  # log2 transformation of counts
  column_to_rownames(var = "Geneid") %>% # convert to dataframe to preserve Geneid
  limma::removeBatchEffect(
    batch = colData(dds_multi) %>% as_tibble() %>% pull(Sex), # REQUIRES CUSOMIZATION; DO NOT INCLUDE PREDICTOR OF INTEREST
    batch2 = colData(dds_multi) %>% as_tibble() %>% pull(Sample_source), # REQUIRES CUSOMIZATION; DO NOT INCLUDE PREDICTOR OF INTEREST
    covariates = colData(dds_multi) %>% as_tibble() %>% select(Sampleid, Age) |> column_to_rownames(var = "Sampleid"),
    design = colData(dds_multi) %>% as_tibble() %>% model.matrix(simple_formula, data = .) # INCLUDE ONLY PREDICTOR OF INTEREST
  ) %>%
  as_tibble(rownames="Geneid") %>% # convert back to tibble
  mutate_at(2:ncol(.), ~ 2^(.)) # remove log2 transformation
# # Export adjusted RPKMs
# rpkms_data_long %>%
#   inner_join(
#     rpkm_SexAgeSource_adj %>%
#       pivot_longer(-Geneid, names_to = "Sampleid", values_to = "RPKM_adj")
#   ) %>%
#   mutate(adjustment = "SexAgeSource") %>%
#   inner_join(gene_anno %>% select(Geneid, gene_name, chr, gene_type)) %>%
#   write_tsv(file = (rpkms_file |> str_replace("RPKMs.txt.gz", "RPKMs_adj.txt.gz")))
# #

## Get genes of interest ----
# usually run only for preferred model
top_signif <- res_multi_SexAgeSource_T21_vs_Control %>%
  filter(padj < 0.1) %>%
  slice_min(pvalue, n = 10) %>%
  arrange(-FoldChange_adj) %>%
  .[1:10,] %>% # ensure only 10 as slice_ will sometimes have ties
  select(Geneid, Gene_name)
#
## Sina plots ----
s1 <- rpkm_SexAgeSource_adj %>%
  pivot_longer(-Geneid, names_to = "Sampleid", values_to = "RPKM") |> 
  inner_join(meta_data) %>%
  inner_join(top_signif) %>%
  mutate(Gene_name = fct_relevel(Gene_name, top_signif %>% pull(Gene_name))) %>% # control plotting order
  group_by(Geneid) %>%
  mutate(extreme = rstatix::is_extreme((log2(RPKM)))) %>%
  filter(extreme != TRUE) %>%
  ungroup() %>%
  ggplot(aes(Karyotype, log2(RPKM), color = Karyotype)) + # CUSTOMIZE
  geom_sina(maxwidth = 0.5) +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.5, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Gene_name, scales = "free", nrow = 2) +
  scale_color_manual(values = standard_colors) +
  theme(aspect.ratio = 1.3) + # set fixed aspect ratio
  labs(
    title = "HTP Whole Blood RNAseq: top significant genes",
    subtitle = "SexAgeSource-adjusted; Extreme outliers removed",
    x = NULL
  )
s1
ggsave(s1, filename = here("plots", paste0(out_file_prefix, "top_signif_RPKM_sina", ".png")), width = 15, height = 5, units = "in")
# May want to rasterize before saving as pdf
s1r <- ggrastr::rasterize(s1, layers='Point', dpi = 600, dev = "ragg_png") 
ggsave(s1r, filename = here("plots", paste0(out_file_prefix, "top_signif_RPKM_sina", ".pdf")), device = cairo_pdf, width = 15, height = 5, units = "in")
#


# 13. Heatmaps ----
#
## Generate gene-wise Z-scores ----
# usually need to transform in some way to compensate for wide range of expression levels
zscores_SexAgeSource_adj <- rpkm_SexAgeSource_adj %>%
  pivot_longer(-Geneid, names_to = "Sampleid", values_to = "RPKM") |> 
  group_by(Geneid) |> 
  mutate(
    zscore = (log2(RPKM) - mean(log2(RPKM), na.rm = TRUE)) / sd(log2(RPKM), na.rm = TRUE)
  ) |> 
  ungroup()
#
## Plot genes of interest as heatmap
# collect data for heatmap
hm_dat <- zscores_SexAgeSource_adj |> 
  inner_join(gene_anno) |> 
  inner_join(meta_data) %>%
  inner_join(top_signif) %>%
  mutate(gene_name = fct_relevel(Gene_name, top_signif %>% pull(Gene_name))) # control plotting order
# Generate centered color scale
hm_lim <- hm_dat |> 
  pull(zscore) |> 
  abs() |> 
  max() |> 
  round(2)
breaks <- seq(-hm_lim, hm_lim, length.out = 11)
hm_palette <- circlize::colorRamp2(
  breaks, 
  RColorBrewer::brewer.pal(11, "RdBu") |> rev()
  )
# plot heatmap
hm_dat |> 
  group_by(Karyotype) |> # to split heatmap
  tidyHeatmap::heatmap(
    gene_name, # be careful of non-unique
    Sampleid,
    zscore,
    palette_value = hm_palette,
    heatmap_legend_param = list(color_bar = "continuous", at = seq(-hm_lim, hm_lim, length.out = 5)),
    cluster_rows = FALSE,
    row_title = NULL,
    show_column_names = FALSE,
    column_title = NULL,
    border = TRUE
  ) |> 
  wrap_heatmap() +
  labs(
    title = "HTP Whole Blood RNAseq: top significant genes",
    subtitle = "SexAgeSource-adjusted Z-scores",
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "top_signif_zcore_heatmap", ".pdf")), device = cairo_pdf, width = 10, height = 5, units = "in")
#

# 14. GSEA Hallmarks analysis ----
#
# Human
hallmarks <- here("data/GSEA/human/h.all.v7.4.symbols.gmt") %>%
  fgsea::gmtPathways(gmt.file = .)
# # Mouse
# hallmarks <- here("data/GSEA/mouse/mh.all.v2022.1.Mm.symbols.gmt") %>%
#   fgsea::gmtPathways(gmt.file = .)
#
## Generate ranks ----
ranks_multi_SexAgeSource_T21_vs_Control <- res_multi_SexAgeSource_T21_vs_Control %>%
  filter(!is.na(log2FoldChange_adj)) %>% # need to remove NA rows that will break plotEnrichment2()
  select(ID = Gene_name, t = log2FoldChange_adj) %>%
  arrange(-abs(t)) %>% # to keep strongest of any duplicates
  distinct(ID, .keep_all = TRUE) %>% # to avoid duplicates
  tibble::deframe() # convert to named numerical vector
#
## Run fgsea ----
hallmarks_multi_SexAgeSource_T21_vs_Control <- run_fgsea2(geneset = hallmarks, ranks = ranks_multi_SexAgeSource_T21_vs_Control, weighted = FALSE)
#
weighted_hallmarks_multi_SexAgeSource_T21_vs_Control <- run_fgsea2(geneset = hallmarks, ranks = ranks_multi_SexAgeSource_T21_vs_Control, weighted = TRUE)
#
## Export GSEA results ----
list(
  "Unweighted" = hallmarks_multi_SexAgeSource_T21_vs_Control %>% unnest(leadingEdge) %>% group_by(pathway, pval, padj, log2err, ES, NES, size) %>% summarize(leadingEdge = paste(leadingEdge, collapse = ",")) %>% arrange(padj, -abs(NES)),
  "Weighted" = weighted_hallmarks_multi_SexAgeSource_T21_vs_Control %>% unnest(leadingEdge) %>% group_by(pathway, pval, padj, log2err, ES, NES, size) %>% summarize(leadingEdge = paste(leadingEdge, collapse = ",")) %>% arrange(padj, -abs(NES))
) |> 
  export_excel(filename = "GSEA_Hallmarks")
#
## GSEA barplot(s) ----
hallmarks_multi_SexAgeSource_T21_vs_Control %>%
  filter(padj < 0.1) %>%
  slice_max(order_by = abs(NES), n = 20) %>%
  arrange(NES) %>%
  mutate(
    pathway = str_remove(pathway, "^HALLMARK_") %>% str_replace_all("_", " ") %>% str_to_title(),
    pathway = fct_inorder(pathway)
  ) %>%
  ggplot(aes(-log10(padj), pathway, fill = NES)) +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_col(color = "black") +
  scale_fill_gradient2(
    low = "#542788",
    mid = "#f7f7f7",
    high = "#b35806",
    midpoint = 0,
    guide = "colourbar",
  ) +
  labs(
    title = "HTP Whole Blood RNAseq: T21 vs. Control\nGSEA Hallmarks (Unweighted, Top 20 q < 0.1)",
    y = NULL
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "barplot_Hallmarks_sig_top20", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#
## GSEA enrichment plot(s) ----
plotEnrichment2(
  pathway = hallmarks$HALLMARK_INTERFERON_GAMMA_RESPONSE,
  stats = ranks_multi_SexAgeSource_T21_vs_Control,
  res = hallmarks_multi_SexAgeSource_T21_vs_Control,
  title = "HTP Whole Blood RNAseq: T21 vs. Control\nInterferon Gamma Response (Unweighted)"
)
#



#+ save_workspace, warning = FALSE, message = FALSE, collapse = TRUE
# ##################
# # Save workspace #
# ##################
save.image(file = here("rdata", paste0(out_file_prefix, ".RData")), compress = TRUE, safe = TRUE) # saves entire workspace (can be slow)
# ################
# # RESTART HERE #
# ################
# load(here("rdata", paste0(out_file_prefix, ".RData")))


#' ***
#' ### Session info 
#+ session_info, collapse=TRUE
# Session info  ---------------------------------------------------
# Report generated at:
date()
cat("\nSessionInfo:")
sessionInfo()


