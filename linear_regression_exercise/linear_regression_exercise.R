################################################
# Title: Linear regression exercise using simple built-in data
# Project: Linear regression
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
library("ggforce")
library("skimr") # data summary
library("conflicted")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("count", "dplyr")
#

# 1 Inspect data ----
# Fuel economy data from 1999 to 2008 for 38 popular models of cars
mpg
mpg |> skim()
#
# Create new variable
mpg <- mpg |> 
  # mutate(trans_type = str_extract(trans, "\\w+(?=\\()")) %>% 
  mutate(
    trans_type = case_when(
      str_detect(trans, "auto") ~ "auto",
      str_detect(trans, "manual") ~ "manual",
    )
  )
#



# 2 Data exploration  ----
mpg |> 
  ggplot(aes(displ, hwy)) +
  geom_point() +
  theme(aspect.ratio = 1) +
  labs(title = "Displacement vs. Highway mpg (hwy)")

mpg |> 
  ggplot(aes(displ, hwy)) +
  geom_point() +
  theme(aspect.ratio = 1) +
  labs(title = "Displacement vs. Highway mpg (hwy)") +
  geom_smooth(method = "lm")


# 3 Linear regression modelling (continuous dependent and predictor) ----
# Create and inspect model object: hwy vs. displ
fit <- lm(hwy ~ displ, data = mpg)
fit
fit |> class()
fit |> str()

# Get results summary
fit |> summary()

# Get tidy model information
fit |> broom::tidy()
fit |> broom::glance()
fit |> broom::augment()

# Compare to model with transformed data: hwy vs. displ
fit2 <- lm(log2(hwy) ~ log2(displ), data = mpg)
fit2 |> summary()
fit2 |> broom::tidy()
fit2 |> broom::glance()


# 4 Linear regression modelling (categorical dependent and continuous predictor) ----
mpg %>% 
  ggplot(aes(trans_type, hwy)) +
  geom_point()
mpg %>% 
  ggplot(aes(trans_type, hwy, color = trans_type)) +
  geom_sina()
mpg %>% 
  ggplot(aes(trans_type, hwy, color = trans_type)) +
  geom_sina() +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75)
mpg %>% 
  ggplot(aes(trans_type, hwy, color = trans_type)) +
  geom_sina() +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  geom_smooth(method = "lm", aes(group = "1"))

# Create and inspect model object: hwy vs. htrans_typewy
fit3 <- lm(hwy ~ trans_type, data = mpg)
fit3 |> summary()
fit3 |> broom::tidy()
fit3 |> broom::glance()

# Compare to model with transformed data
fit4 <- lm(log2(hwy) ~ trans_type, data = mpg)
fit4 |> summary()
fit4 |> broom::tidy()
fit4 |> broom::glance()


# Compare to t test
ttest <- t.test(hwy ~ trans_type, data = mpg)
ttest |> class()
ttest
ttest |> broom::tidy() # 3.38e-05
fit3 |> broom::tidy() # 1.89e-5



