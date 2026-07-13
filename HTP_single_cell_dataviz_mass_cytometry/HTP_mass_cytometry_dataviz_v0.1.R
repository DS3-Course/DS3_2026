################################################
# Title: HTP Mass Cytometry Dataviz
# Project: Analysis of mass cytometry data in R
# Author(s):
#   - author 1
# email(s):
#   - author1@institution.edu
# affiliation(s):
#   - Department of XXX
#   - University of XXX
################################################
# Script original author: Matthew Galbraith
# version: 0.1  Date: 05_24_2024
################################################
# Change Log:
# v0.1
# Initial version
#

### Summary:  
# Data visualization of CD45+CD66lo-gated mass cytometry (CyTOF) data from 388 HTP participants.
# See PMID 37379383.
# See also:
# https://www.nature.com/articles/s41467-019-13055-y
# https://distill.pub/2016/misread-tsne/
# https://pair-code.github.io/understanding-umap/


### Data type(s):
#   A. HTP sample meta data
#      Where/who did this data come from?
#      What is the source of the original data and where is it stored?
#   B. HTP CyTOF tSNE data
#.     tSNE run on 500 cells per sample, perplexity set at 440
#      Where/who did this data come from?
#      What is the source of the original data and where is it stored?
#

### Workflow:
#   Step 1 - Read in and inspect sample meta data + CyTOF data
#   Step 2 - Data exploration
#   Step 3 - 
#   Step 4 - 
#   Step 5 - 
#  

## Comments:
#  Any further relevant details?  
#  


# 0 General Setup -----
# RUN FIRST TIME
# renv::init()
## 0.1 Load required libraries ----
library("readxl") # used to read .xlsx files
library("openxlsx") # used for data export as Excel workbooks
library("tidyverse") # data wrangling and ggplot2
# library("rstatix") # pipe- and tidy-friendly statistical tests
library("ggrepel") # for labelling genes
library("ggforce") # for sina plots
library("tictoc") # timer
library("skimr") # data summary
library("janitor") # data cleaning
library("patchwork") # assembling multiple plots
library("ggrastr") # required for rasterizing some layers of plots
library("conflicted")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("count", "dplyr")
library("here") # generates path to current project directory
# detach("package:here", unload=TRUE) # run this to reset here()
source(here("helper_functions.R")) # load helper functions
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

## 0.3 Set required parameters ----
# Input data files
htp_meta_data_file <- here("data", "HTP_Metadata_v0.5_Synapse.txt") # comments/notes on this file?
#
CD45_cells_clusters_data_file <- here("data", "HTP_CyTOF_CD45_FlowSOMv1_clusters.txt.gz") # FlowSOM clustering assignments for 500 cells per sample
CD45_cells_cluster_24_info_file <- here("data", "HTP_CyTOF_CD45_FlowSOMv1_cell_cluster_24_info.txt") # named clusters metadata
CD45_cells_tSNE_data_file <- here("data", "HTP_CyTOF_CD45_FlowSOMv1_tSNE.txt.gz") # 2D tSNE coordinates for 500 cells per sample (perplexity = 440)
CD45_cells_marker_expression_data_file <- here("data", "HTP_CyTOF_CD45_FlowSOMv1_markers.txt.gz") # Marker expression levels for 500 cells per sample ± arcsinh transformation
CD45_betreg_results_file <- here("data", "DS3_CyTOF_betareg_v0.1_betareg_multi_SexAgeSource_results.txt")
#
# Other parameters
# standard_colors <- c("Group1" = "#F8766D", "Group2" = "#00BFC4")
standard_colors <- c("Control" = "gray60", "T21" = "#009b4e")
out_file_prefix <- "HTP_mass_cytometry_dataviz_v0.1_"
# End required parameters ###


# 1 Read in and inspect data ----
## 1.1 Read in meta data ----
htp_meta_data <- htp_meta_data_file |> 
  read_tsv() |> 
  mutate(
    Karyotype = fct_relevel(Karyotype, c("Control", "T21")), # convert to factor and set order
    Sex = fct_relevel(Sex, "Female"), # convert to factor and set order
    Sample_source_code = as_factor(Sample_source_code) # convert to factor - default is numerical order
  )
# inspect
htp_meta_data
htp_meta_data |> skimr::skim()
#
here("data/HTP_Metadata_v0.5_dictionary.txt") |> read_tsv()
#

## 1.2 Read in cluster assignments ----
CD45_cells_clusters_data <- CD45_cells_clusters_data_file |> 
  read_tsv() |> 
  mutate(
    cell = as_factor(cell), # to prevent being treated as a number
    som100_cluster = as_factor(som100_cluster),
    meta20_cluster = as_factor(meta20_cluster),
    meta30_cluster = as_factor(meta30_cluster)
  )
# inspect
CD45_cells_clusters_data # 194,000 rows = 388 LabIDs x 500 cells each
CD45_cells_clusters_data |> skimr::skim()
CD45_cells_clusters_data |> distinct(LabID) # 388 LabIDs
CD45_cells_clusters_data |> group_by(LabID) |> summarize("n_cells" = n_distinct(cell)) # 500 cells per sample
CD45_cells_clusters_data |> group_by(cell) # unique
CD45_cells_clusters_data |> group_by(Time, Event_length) # unique
CD45_cells_clusters_data |> distinct(som100_cluster) # 100 SOM clusters
CD45_cells_clusters_data |> distinct(meta20_cluster) # 20 meta clusters
CD45_cells_clusters_data |> distinct(meta30_cluster) # 30 meta clusters
CD45_cells_clusters_data |> distinct(Cell_cluster_name) # 24 cell type labels
CD45_cells_clusters_data |> filter(!str_detect(Cell_cluster_name, "EXCLUDE")) |> distinct(Cell_cluster_name) # 20 cell type labels
# look at merging of clusters
CD45_cells_clusters_data |> distinct(meta20_cluster, meta30_cluster, Cell_cluster_name) |> arrange(meta20_cluster)
#

### 1.2.1 Read in cluster_24_info
CD45_cells_cluster_24_info <- CD45_cells_cluster_24_info_file |> read_tsv()
#


## 1.3 Read in tSNE data ----
CD45_cells_tSNE_data <- CD45_cells_tSNE_data_file |> 
  read_tsv() |> 
  mutate(
    cell = as_factor(cell) # to prevent being treated as a number
  )
# inspect
CD45_cells_tSNE_data # 194,000 rows = 388 LabIDs x 500 cells each
CD45_cells_tSNE_data |> skimr::skim()
CD45_cells_tSNE_data |> distinct(LabID) # 388 LabIDs
CD45_cells_tSNE_data |> group_by(LabID) |> summarize("n_cells" = n_distinct(cell)) # 500 cells per sample
#


## 1.3 Read in marker expression data ----
CD45_cells_marker_expression_data <- CD45_cells_marker_expression_data_file |> 
  read_tsv() |> 
  mutate(
    cell = as_factor(cell) # to prevent being treated as a number
  )
# inspect
CD45_cells_marker_expression_data # 6,984,000 rows = 388 LabIDs x 500 cells each x 36 markers (features)
CD45_cells_marker_expression_data |> skimr::skim()
CD45_cells_marker_expression_data |> distinct(LabID) # 388 LabIDs
CD45_cells_marker_expression_data |> distinct(feature) # 36 markers
CD45_cells_marker_expression_data |> group_by(LabID) |> summarize("n_cells" = n_distinct(cell)) # 500 cells per sample
#


## 1.4 Read in T21 vs. Control fold-changes ----
CD45_betreg_results <- CD45_betreg_results_file |> 
  read_tsv()
CD45_betreg_results
#


## 1.5 Check joining against meta data ----
CD45_cells_clusters_data |> 
  inner_join(htp_meta_data)
CD45_cells_tSNE_data |> 
  inner_join(htp_meta_data)
CD45_cells_marker_expression_data |> 
  inner_join(htp_meta_data)
# check number of rows returned !!!
CD45_cells_clusters_data |> 
  inner_join(htp_meta_data) |> 
  distinct(LabID, Karyotype, Sex) |> 
  count(Karyotype)
CD45_cells_clusters_data |> 
  inner_join(htp_meta_data) |> 
  distinct(LabID, Karyotype, Sex) |> 
  count(Sex)
CD45_cells_clusters_data |> 
  inner_join(htp_meta_data) |> 
  distinct(LabID, Karyotype, Sex) |> 
  count(Karyotype, Sex)
#



# 2 tSNE plots  ----

## 2.1 Dealing with overploting ----
#
# High amount of overplotting
CD45_cells_tSNE_data |> 
  ggplot(aes(TSNE1, TSNE2)) +
  geom_point() +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "overplotting of points"
  )
# reduce point size
CD45_cells_tSNE_data |> 
  ggplot(aes(TSNE1, TSNE2)) +
  geom_point(size = 0.1) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "reduced point size"
  )
# sub-sampling
CD45_cells_tSNE_data |> 
  group_by(LabID) |> 
  slice_sample( n = 100) |> 
  ungroup() |> 
  ggplot(aes(TSNE1, TSNE2)) +
  geom_point(size = 0.1) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "reduced point size + 100 cells per sample"
  )
# transparancy (alpha)
# loses low density regions and/or still too much over plotting
# also looks fuzzy
CD45_cells_tSNE_data |> 
  ggplot(aes(TSNE1, TSNE2)) +
  geom_point(size = 0.1, alpha = 0.1) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "reduced point size + alpha = 0.1"
  )
CD45_cells_tSNE_data |> 
  ggplot(aes(TSNE1, TSNE2)) +
  geom_point(size = 0.1, alpha = 0.01) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "reduced point size + alpha = 0.01"
  )
# Color by density
renv::install("LKremer/ggpointdensity") 
library("ggpointdensity")
CD45_cells_tSNE_data |> 
  ggplot(aes(TSNE1, TSNE2)) +
  geom_pointdensity(size = 0.1) + # SLOW with all data points
  scale_color_viridis_c() +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "reduced point size + geom_pointdensity"
  )
# Custom density function
CD45_cells_tSNE_data |> 
  mutate(
    density = getDenCols(TSNE1, TSNE2, transform = FALSE) 
  ) |> 
  arrange(density) |> 
  ggplot(aes(TSNE1, TSNE2, color = density)) +
  geom_point(size = 0.1) +
  scale_color_viridis_c() +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "reduced point size + custom density function"
  )
#


## 2.2 Comparing two groups ----
# be careful with plotting order
CD45_cells_tSNE_data |> 
  inner_join(htp_meta_data) |> 
  ggplot(aes(TSNE1, TSNE2, color = Karyotype)) +
  geom_point(size = 0.1) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "Levels: Control T21"
  )
CD45_cells_tSNE_data |> 
  inner_join(htp_meta_data) |> 
  mutate(Karyotype = fct_relevel(Karyotype, c("T21", "Control"))) |> 
  ggplot(aes(TSNE1, TSNE2, color = Karyotype)) +
  geom_point(size = 0.1) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "Levels: T21 Control"
  )
# better to facet
# but how to compare when n is different?
CD45_cells_tSNE_data |> 
  inner_join(htp_meta_data) |> 
  distinct(LabID, Karyotype) |> 
  count(Karyotype)
#
CD45_cells_tSNE_data |> 
  inner_join(htp_meta_data) |> 
  ggplot(aes(TSNE1, TSNE2, color = Karyotype)) +
  geom_point(size = 0.1) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = ""
  ) +
  facet_wrap(~ Karyotype) +
  scale_color_manual(values = standard_colors)
#
# Color by density + split by Karyotype (but same color scale?)
CD45_cells_tSNE_data |> 
  inner_join(htp_meta_data) |> 
  ggplot(aes(TSNE1, TSNE2)) +
  geom_pointdensity(size = 0.1) +
  scale_color_viridis_c() +
  facet_wrap(~ Karyotype) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "reduced point size + geom_pointdensity"
  )
# Split by Karyotype (~ distinct color scales)
CD45_cells_tSNE_data |> 
  inner_join(htp_meta_data) |> 
  group_by(Karyotype) |> 
  mutate(
    density = getDenCols(TSNE1, TSNE2, transform = FALSE) 
  ) |> 
  arrange(density) |> # note sorting by density
  ggplot(aes(TSNE1, TSNE2, color = density)) +
  geom_point(size = 0.1) +
  scale_color_viridis_c() +
  facet_wrap(~ Karyotype) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "reduced point size + custom density function"
  )
# Normalized density per facet
# currently requires dev version of ggpoint density
CD45_cells_tSNE_data |> 
  inner_join(htp_meta_data) |> 
  ggplot(aes(TSNE1, TSNE2)) +
  # geom_pointdensity(size = 0.1) +
  ggpointdensity::stat_pointdensity(aes(col = after_stat(ndensity)), size = 0.1) + 
  scale_color_viridis_c() +
  facet_wrap(~ Karyotype) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "reduced point size + normalized density"
  )
# Still better if sorted by density
CD45_cells_tSNE_data |> 
  inner_join(htp_meta_data) |> 
  group_by(Karyotype) |> 
  mutate(
    density = getDenCols(TSNE1, TSNE2, transform = FALSE) 
  ) |> 
  arrange(density) |> # sorting by density
  ggplot(aes(TSNE1, TSNE2)) +
  # geom_pointdensity(size = 0.1) +
  ggpointdensity::stat_pointdensity(aes(col = stat(ndensity)), size = 0.1) + 
  scale_color_viridis_c() +
  facet_wrap(~ Karyotype) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "reduced point size + normalized density"
  )
# Alternative: plot equal sample numbers per group
CD45_tSNE_data_sample_subset <- CD45_cells_tSNE_data |> 
  inner_join(htp_meta_data) |> 
  distinct(LabID, Karyotype, Sex) |> 
  group_by(Karyotype, Sex) |> 
  slice_sample( n = 43) |> 
  ungroup()
CD45_tSNE_data_sample_subset |> 
  count(Karyotype, Sex)
#
CD45_cells_tSNE_data |> 
  inner_join(htp_meta_data) |> 
  filter(LabID %in% CD45_tSNE_data_sample_subset$LabID) |>
  ggplot(aes(TSNE1, TSNE2)) +
  geom_pointdensity(size = 0.1) +
  scale_color_viridis_c() +
  facet_wrap(~ Karyotype) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "equal n per group + geom_pointdensity"
  )
#


## 2.3 Coloring by clusters ----
# too many clusters
CD45_cells_tSNE_data |> 
  inner_join(CD45_cells_clusters_data) |> 
  group_by(LabID) |> 
  slice_sample( n = 100) |> 
  ggplot(aes(TSNE1, TSNE2, color = som100_cluster)) +
  geom_point(size = 0.1) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "FlowSOM 100 clusters"
  )
# more sensible: 30 meta clusters
# still too many colors for default ggplot color scale
CD45_cells_tSNE_data |> 
  inner_join(CD45_cells_clusters_data) |> 
  group_by(LabID) |> 
  slice_sample( n = 100) |> 
  ggplot(aes(TSNE1, TSNE2, color = meta30_cluster)) +
  geom_point(size = 0.1) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "FlowSOM 30 meta clusters"
  )
# cell cluster names
# still too many colors for default ggplot color scale
CD45_cells_tSNE_data |> 
  inner_join(CD45_cells_clusters_data) |> 
  filter(!str_detect(Cell_cluster_name, "EXCLUDE")) |> # remove excluded clusters
  group_by(LabID) |> 
  slice_sample( n = 100) |> 
  ggplot(aes(TSNE1, TSNE2, color = Cell_cluster_name)) +
  geom_point(size = 0.1) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "FlowSOM named clusters"
  )
#

### 2.3.1 Generating color schemes using I Want Hue ----
# https://medialab.github.io/iwanthue/
# see https://cran.r-project.org/web/packages/hues/index.html
library("hues")
iwanthue(5)
iwanthue(5, plot = TRUE)
# see also https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3
# see also https://cran.r-project.org/web/packages/RColorBrewer/index.html
#
CD45_cells_tSNE_data |> 
  inner_join(CD45_cells_clusters_data) |> 
  filter(!str_detect(Cell_cluster_name, "EXCLUDE")) |> # remove excluded clusters
  group_by(LabID) |> 
  slice_sample( n = 100) |> 
  ggplot(aes(TSNE1, TSNE2, color = Cell_cluster_name)) +
  geom_point(size = 0.1) +
  scale_color_iwanthue() +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "improved color palette"
  )
# see also iwanthue() function in https://rdocumentation.org/packages/SC3/versions/1.1.4
#
### 2.3.2 Using a named color scheme ----
# better control if we save and name colors
hues::iwanthue(20)
cell_cluster_colors <- c("#636ad8","#60b646","#af5dcb","#a9b539","#d14396","#5bbe7d","#d34258","#4bb7a7","#d0522d","#5faadc","#d69b3a","#5778be","#757327","#9b72b8","#447d41","#de87bb","#b2ad6a","#9f4866","#a26431","#de8373", "grey", "grey", "grey", "grey")
names(cell_cluster_colors) <- CD45_cells_cluster_24_info |> arrange(num) |> pull(label)
cell_cluster_colors
#
t20 <- CD45_cells_tSNE_data |> 
  inner_join(htp_meta_data) |> 
  inner_join(CD45_cells_clusters_data) |> 
  inner_join(CD45_cells_cluster_24_info) |> 
  # To put less abundant clusters on top:
  arrange(num) |> # sort by overall cluster percentage
  mutate(label = fct_inorder(label)) |>  # to control plotting order
  filter(!str_detect(Cell_cluster_name, "EXCLUDE")) |> # remove excluded clusters
  group_by(LabID) |> 
  slice_sample( n = 100) |> 
  ggplot(aes(TSNE1, TSNE2, color = label)) +
  geom_point(size = 0.1) +
  scale_color_manual(values = cell_cluster_colors) +
  theme(aspect.ratio = 1) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "named color palette"
  )
t20
#
# tidy up color legend
t20 <- t20 +
  guides(
    colour = guide_legend(
      override.aes = list(size=2), 
      ncol = 2,
      title = "Cell cluster"
      )
    )
t20
#
### 2.3.3 Adding cluster labels ----
cell_cluster_tSNE_labels <- CD45_cells_tSNE_data |> 
  inner_join(htp_meta_data) |> 
  inner_join(CD45_cells_clusters_data) |> 
  group_by(Cell_cluster_name) |> 
  summarize(
    TSNE1 = median(TSNE1), # set median xy coordinates of each cluster as label position
    TSNE2 = median(TSNE2)
  ) |> 
  inner_join(CD45_cells_cluster_24_info) # numbered by overall percentage
cell_cluster_tSNE_labels
#
t20 <- t20 +
  geom_text(
    data = cell_cluster_tSNE_labels |> filter(!str_detect(Cell_cluster_name, "EXCLUDE")), # remove excluded clusters
    aes(label = num), 
    color = "black"
    ) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "named color palette + cluster labels"
  )
t20
#
t20 + facet_wrap(~ Karyotype)
#


## 2.4 Coloring by expression values ----
t_expr <- CD45_cells_tSNE_data |> 
  inner_join(htp_meta_data) |> 
  inner_join(CD45_cells_clusters_data) |> 
  inner_join(CD45_cells_cluster_24_info) |> 
  filter(!str_detect(Cell_cluster_name, "EXCLUDE")) |> # remove excluded clusters
  inner_join(CD45_cells_marker_expression_data) |> 
  filter(feature %in% c("CD3", "CD4", "CD8a", "CD19")) |> 
  group_by(feature) %>% 
  mutate(
    trimmed = scales::oob_squish(asinh5_exprs, range = c(quantile(asinh5_exprs, 0.0025), quantile(asinh5_exprs, 0.9975))),
    scaled = scales::rescale(asinh5_exprs),
    scaled_trimmed = scales::rescale(trimmed)
  ) %>% 
  arrange(scaled) %>% 
  ungroup() %>% 
  ggplot(aes(TSNE1, TSNE2, color = scaled_trimmed)) +
  geom_point(size = 0.05) +
  scale_color_viridis_c() +
  facet_wrap(~ feature, nrow = 1) +
  theme(
    aspect.ratio = 1,
    legend.key=element_blank(),
    axis.text.x=element_blank(), axis.text.y=element_blank()
  ) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "Scaled marker expression; 0.25%ile trim"
  )
t_expr
#


## 2.5 Coloring by fold-change ----
#
# Get color scale limits
fc_limits <- CD45_betreg_results |> 
  filter(BHadj_pval < 0.1) |> 
  pull(FoldChange) |> 
  log2() |> 
  abs() |> 
  max() |> 
  round(2) * c(-1, 1) # ensures color scale is centered on 0
#
t_betareg_fc <- CD45_cells_tSNE_data |> 
  inner_join(CD45_cells_clusters_data) |> 
  inner_join(CD45_cells_cluster_24_info) |> 
  filter(!str_detect(Cell_cluster_name, "EXCLUDE")) |> # remove excluded clusters
  inner_join(CD45_betreg_results) |> 
  ggplot(aes(TSNE1, TSNE2, color)) +
  geom_point( # plot n.s. clusters in grey
    data = . %>% filter(is.na(BHadj_pval) | BHadj_pval >= 0.1),
    size = 0.05, color = "grey"
  ) +
  geom_point( # plot significant cluster colored by fold-change
    data = . %>% filter(BHadj_pval < 0.1),
    aes(color = log2(FoldChange)),
    size = 0.05
  ) +
  # customize color scale to center on 0
  scale_color_distiller(palette = "RdBu", limit = fc_limits, oob = scales::squish) +
  theme(
    aspect.ratio = 1,
    legend.key=element_blank(),
    axis.text.x=element_blank(), axis.text.y=element_blank()
  ) +
  labs(
    title = "HTP Mass Cytometry: CD45+CD66lo",
    subtitle = "Significant clusters colored by fold-change in T21 vs. D21"
  )
t_betareg_fc
#

## 2.6 Exporting plots with rasterized points ----
# Editing plots outside R generally required for finalizing figures (eg labels, font sizes)
# Plots with 100k+ points likely to cause problems and usually do not require editing
# Exporting as PDF with rasterized points layers preserves editability of labels, axes etc
#
t20_rast <- ggrastr::rasterize(t20, layers='Point', dpi = 600, dev = "ragg_png")
ggsave(t20_rast, filename = here("plots", paste0(out_file_prefix, "tSNE_plot_clusters", ".pdf")), device = cairo_pdf, width = 8, height = 5, units = "in")
#
t_expr_rast <- ggrastr::rasterize(t_expr, layers='Point', dpi = 600, dev = "ragg_png")
ggsave(t_expr_rast, filename = here("plots", paste0(out_file_prefix, "tSNE_plot_markers", ".pdf")), device = cairo_pdf, width = 16, height = 5, units = "in")
#
t_betareg_fc_rast <- ggrastr::rasterize(t_betareg_fc, layers='Point', dpi = 600, dev = "ragg_png")
ggsave(t_betareg_fc_rast, filename = here("plots", paste0(out_file_prefix, "tSNE_plot_foldchange", ".pdf")), device = cairo_pdf, width = 8, height = 5, units = "in")
#


################################################
# save workspace ----
save.image(file = here("rdata", paste0(out_file_prefix, ".RData")), compress = TRUE, safe = TRUE) # saves entire workspace (can be slow)
# To reload previously saved workspace:
# load(here("rdata", paste0(out_file_prefix, ".RData")))

# session_info ----
date()
sessionInfo()
################################################
