################################################
# Title: Rproject analysis template
# Project: Rprojects Introduction
# Author(s):
#   - author 1
# email(s):
#   - author1@institution.edu
# affiliation(s):
#   - Department of XXX
#   - University of XXX
################################################
# Script original author: Matthew Galbraith
# version: 0.1  Date: 05_08_2024
################################################
# Change Log:
# v0.1
# Initial version
#

### Summary:  
# Description of data wrangling and/or analysis being performed.
#  

### Data type(s):
#   A. Meta data
#      Where/who did this data come from?
#      What is the source of the original data and where is it stored?
#   B. Data type 1
#      Where/who did this data come from?
#      What is the source of the original data and where is it stored?
#   C. Data type 2
#      Where/who did this data come from?
#      What is the source of the original data and where is it stored?
#

### Workflow:
#   1. Step 1 description
#   2. Step 2 description
#   3. Step 3 description
#   4. Step 4 description
#

## Comments:
#  Any further relevant details?  
#   


# 0 General Setup -----
# Initialize and install packages with renv
# renv::init()
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

## 0.1 Load required libraries ----
library("tidyverse")
library("readxl") # read .xlsx files
library("openxlsx") # data export to Excel workbooks
library("skimr") # data summary and validation
library("janitor") # data cleaning etc
# library("ggrepel") # labelling points in plots
# library("ggforce") # sina plots etc 
# library("patchwork") # arranging plots
# library("tidyHeatmap") # tidy interface to ComplexHeatmap
# library("plotly") # generating interactive plots
# library("tictoc") # timer
library("conflicted") # force all conflicts to become errors
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("count", "dplyr")
library("here") # generate path to current project directory
#
source(here("helper_functions.R")) # load helper functions
#


## 0.2 Set required parameters ----
# Input data files
meta_data_file <- here("data", "meta_data_file.txt") # comments/notes on this file?
data_type1_file <- here("data", "data_type1_file.txt") # comments/notes on this file?
data_type2__file <- here("data", "data_type1_file.txt")  # comments/notes on this file?
# Other parameters
standard_colors <- c("Group1" = "#F8766D", "Group2" = "#00BFC4")
out_file_prefix <- "analysis_script_template_v0.1_"
# End required parameters ###


# 1 Read in and inspect data ----
## 1.1 Read in meta data ----
meta_data <- meta_data_file |> 
  read_tsv() |> 
  janitor::clean_names(case = "none")
# inspect
meta_data
meta_data |> skimr::skim()
#

## 1.2 Read in data type 1 ----
data_type1 <- data_type1_file |> 
  read_tsv() |> 
  janitor::clean_names(case = "none")
# inspect
data_type1
data_type1 |> skimr::skim()
#

## 1.3 Read in data type 2 ----
data_type2 <- data_type2_file |> 
  read_tsv() |> 
  janitor::clean_names(case = "none")
# inspect
data_type2
data_type2 |> skimr::skim()
#

## 1.3 Join meta data with data type 1 and data type 2 ----
combined <- meta_data |> 
  inner_join(data_type1) |> 
  inner_join(data_type2)
# check number of rows returned


# 2 Data exploration / summary ----
# check data distribution(s), outliers etc


# 3 Analysis ----
# statistical testing
# assemble and export results

# 4 Plot results ----
# plot results summaries


# 5 Plot individual features ----
# plot interesting/significant features


################################################
# save workspace ----
save.image(file = here("rdata", paste0(out_file_prefix, ".RData")), compress = TRUE, safe = TRUE) # saves entire workspace (can be slow)
# To reload previously saved workspace:
# load(here("rdata", paste0(out_file_prefix, ".RData")))

# session_info ----
date()
sessionInfo()
################################################
