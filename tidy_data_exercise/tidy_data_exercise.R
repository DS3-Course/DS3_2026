################################################
# Title: Tidy data exercise using HTP examples
# Project: Introduction to tidy data
################################################
# Script original author: Matthew Galbraith
# version: 0.1  Date: 05_13_2024
################################################
# Change Log:
# v0.1
# Initial version
#


# 0 General Setup -----
# RUN FIRST TIME
# renv::init()
## 0.1 Load required libraries ----
library("tidyverse") # data wrangling and ggplot2
library("here")
#


# Inspired by:
# https://tidyr.tidyverse.org/articles/tidy-data.html
# https://vita.had.co.nz/papers/tidy-data.html

# Data semantics
# - A dataset is a collection of values, usually either numbers (if quantitative) or strings (if qualitative).
# - Every value belongs to a variable and an observation.

# 1 The common 'wide' data format ---- 

## 1.1 - 1 column per sample, 1 row per feature ----
# PROBLEM: hard to store data about samples
htp_cytokines_wide <- here("data/htp_cytokines_data_wide.txt.gz") |> 
  read_tsv()
htp_cytokines_wide

## 1.2 - 1 column per feature, 1 row per sample ----
# Same data as above, but the rows and columns have been transposed
# PROBLEM: hard to store data about features
htp_cytokines_wide2 <- here("data/htp_cytokines_data_wide2.txt.gz") |> 
  read_tsv()
htp_cytokines_wide2

## 1.3 Tidy version ----
htp_cytokines_long <- htp_cytokines_wide |> 
  pivot_longer(-Analyte, names_to = "LabID", values_to = "Value")
htp_cytokines_long
# This makes the values, variables, and observations more clear.
# Every combination of name and assessment is a single measured observation
# Also allows for additional data to be associated (although see 2.3 below)

# Tidy data
# 1. Each variable is a column; each column is a variable.
# 2. Each observation is a row; each row is an observation.
# 3. Each value is a cell; each cell is a single value.


# 2 Common messy data problems ----

## 2.1 Column headers are values, not variable names ----
htp_Age_groups <- here("data/HTP_values_in_colnames.txt") |> 
  read_tsv()
htp_Age_groups
# Tidied
htp_Age_groups %>% 
  pivot_longer(-c(Karyotype, Sex), names_to = "Age_group", values_to = "n")
#

## 2.2 Multiple variables stored in one column ----
htp_demo_groups <- here("data/HTP_multiple_vars_one_col.txt") |> 
  read_tsv()
htp_demo_groups
# Tidied
htp_demo_groups |> 
  separate(Demographic_group, into = c("Age_group", "Sex", "Karyotype"))
#

## 2.3 Multiple data types in one table ----
# Problem: some values are repeated
# Problem: could be unclear which variables apply to which entity (Samples vs. Features)
htp_meta_cytokines <- here("data/HTP_meta_cytokines_data.txt.gz") |> 
  read_tsv()
htp_meta_cytokines
# Tidied sample data
htp_meta_cytokines |> 
  select(LabID, ParticipantID:Ethnicity) |> 
  distinct()
# Tidied feature data
htp_meta_cytokines |> 
  select(LabID, Sample_type:N_imputed_wells) |> 
  distinct()

## 2.4 One data type in multiple tables ----
# variables are the same in both files
here("data/HTP_meta_Controls.txt.gz") |> read_tsv()
here("data/HTP_meta_T21.txt.gz") |> read_tsv()
#
# method 1
combined <- bind_rows(
  here("data/HTP_meta_Controls.txt.gz") |> read_tsv(),
  here("data/HTP_meta_T21.txt.gz") |> read_tsv()
)
combined
#
# method 2 - useful when many files
combined2 <- read_tsv(
  list.files(path = here("data"), pattern = "HTP_meta_[Controls|T21]", full.names = TRUE),
  id = "path" # keep track of file paths
)
combined2
#







    