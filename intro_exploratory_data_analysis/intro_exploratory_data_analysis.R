################################################
# Title: Exploratory Exploratory Data Analysis exercise
# Project: Introduction to Exploratory Data Analysis
################################################
# Script original author: Matthew Galbraith
# version: 0.1  Date: 07_10_2025
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
library("ggforce") # sina plots etc 
library("skimr") # data summary and validation
library("here")
#
out_file_prefix <- "intro_exploratory_data_analysis_v0.1_"
#


# 1 Built-in datasets for testing / practicing ----
?mtcars
mtcars
# what type of object is mtcars?

?iris
iris
# what type of object is iris?

?trees
trees
# what type of object is iris?

?mpg
mpg
# what type of object is mpg?


# 2 Convert between data.frames and tibbles ----
## 2.1 data.frame to tibble ± rownames ----
mtcars |> as_tibble()
mtcars |> as_tibble(rownames = "car")
mtcars_tbl <- rownames_to_column(mtcars, var = "car") %>% as_tibble()

iris |> as_tibble()

## 2.2 tibble to data.frame ± rownames ----
mtcars_tbl |> as.data.frame()
mtcars_tbl |> column_to_rownames(var = "car") %>% head()

mpg |> rowid_to_column() |> head()
#


# 3 Review syntax for selection ----
# see also https://rstudio.github.io/cheatsheets/syntax.pdf

## 3.1 Base R 'dollar' syntax ----
mtcars$mpg

## 3.2 Base R 'square backets' extraction operator ----
mtcars[1]
mtcars["mpg"]

# [rows, cols]
mtcars[,1]
# using a range
mtcars[1:10,1]
# using name
mtcars[,"mpg"]

# subset/filter by a condition
mtcars[mtcars$mpg > 30, ]
mtcars[mtcars$cyl == 4, ]


## 3.3 Tidyverse syntax for selection ----
mtcars |> pull(1)
mtcars |> pull(mpg)

mtcars_tbl |> pull(2)
mtcars_tbl |> pull(mpg)

mtcars |> select(mpg) |> head(n = 10)
mtcars_tbl |> select(car, mpg)

# subset/filter by a condition
mtcars_tbl |> filter(mpg > 30)
mtcars_tbl |> filter(cyl == 4)


# 4 Exploratory data analysis ----

## 4.1 Summary statistics using Base R vs. Tidyverse ----

# one continuous variable
mean(mtcars$mpg)
mtcars_tbl |> dplyr::summarize(mean(mpg))
mtcars_tbl |> dplyr::summarize(mean_mpg = mean(mpg)) # better to name new variables

# one categorical variable
table(mtcars$cyl) # no labels!
mtcars_tbl |> 
  dplyr::group_by(cyl) |> 
  dplyr::summarize(n = n())
mtcars_tbl |> 
  dplyr::count(cyl) # even simpler

# two categorical variables
table(mtcars$cyl, mtcars$am)
mtcars_tbl |> 
  dplyr::group_by(cyl, am) |> 
  dplyr::summarize(n = n())
mtcars_tbl |> 
  dplyr::count(cyl, am)

# one continuous, one categorical
mean(mtcars$mpg[mtcars$cyl==4])
mean(mtcars$mpg[mtcars$cyl==6])
mean(mtcars$mpg[mtcars$cyl==8])

mtcars_tbl |> 
  dplyr::group_by(cyl) |> 
  dplyr::summarize(mean_mpg = mean(mpg))

# summarize multiple variables all at once
mtcars_tbl |> 
  dplyr::group_by(cyl) |> 
  dplyr::summarize(
    mean_mpg = mean(mpg),
    mean_disp = mean(disp)
    )

mtcars_tbl |> summary()

mtcars_tbl |> skimr::skim()

# add groups
mtcars_tbl |> 
  dplyr::group_by(cyl) |> 
  skimr::skim()


## 4.2 Data visualization using Base R vs. Tidyverse ----

### 4.2.1 Two continuous variables - scatter plot ----
# Base R
plot(mpg$displ, mpg$hwy)

# ggplot
mpg |> 
  ggplot(aes(displ, hwy)) +
  geom_point()

# add a third variable
mpg |> 
  ggplot(aes(displ, hwy, color = class)) +
  geom_point()

# add a fourth variable
mpg |> 
  ggplot(aes(displ, hwy, color = class, shape = drv)) +
  geom_point()

# add a title / subtitle
mpg |> 
  ggplot(aes(displ, hwy, color = class, shape = drv)) +
  geom_point() +
  labs(
    title = "Highway mpg (hwy) vs. displacement",
    subtitle = "subtitle"
  )

# make the plot look nicer
mpg |> 
  ggplot(aes(displ, hwy, color = class, shape = drv)) +
  geom_point() +
  labs(
    title = "Highway mpg (hwy) vs. displacement",
    subtitle = "subtitle"
  ) +
  theme(aspect.ratio = 1)

mpg |> 
  ggplot(aes(displ, hwy, color = class, shape = drv)) +
  geom_point() +
  labs(
    title = "Highway mpg (hwy) vs. displacement",
    subtitle = "subtitle"
  ) +
  theme(
    aspect.ratio = 1,
    panel.border=element_rect(colour="black", fill="transparent"),
    plot.title=element_text(face="bold", hjust=0),
    axis.text=element_text(color="black", size=14),
    axis.text.x=element_text(angle=0, hjust=0.5),
    axis.ticks = element_line(color = "black"), # make sure tick marks are black
    panel.background=element_blank(),
    panel.grid=element_blank(),
    plot.background=element_blank(),
    strip.background = element_blank(), # facet label borders
    legend.key=element_blank(), legend.background=element_blank() # remove grey bg from legend
  )

# save the plot as PNG
ggsave(filename = here("plots", paste0(out_file_prefix, "mpg_scatter", ".png")), width = 7, height = 5, units = "in")
# save the plot as PDF
ggsave(filename = here("plots", paste0(out_file_prefix, "mpg_scatter", ".pdf")), device = cairo_pdf, width = 7, height = 5, units = "in")


### 4.2.2 Categorical vs. continuous variable - box + sina plots ----
# Base R
plot(as.factor(mpg$class), mpg$hwy)

# ggplot
# Setting and modifying default theme for plots
theme_set(theme_gray(base_size=12, base_family="Arial") +
            theme(
              panel.border=element_rect(colour="black", fill="transparent"),
              plot.title=element_text(face="bold", hjust=0),
              axis.text=element_text(color="black", size=14),
              axis.text.x=element_text(angle=0, hjust=0.5),
              axis.ticks = element_line(color = "black"), # make sure tick marks are black
              panel.background=element_blank(),
              panel.grid=element_blank(),
              plot.background=element_blank(),
              strip.background = element_blank(), # facet label borders
              legend.key=element_blank(), legend.background=element_blank() # remove grey bg from legend
            )
)

# Box plot
mpg %>% 
  ggplot(aes(class, hwy)) +
  geom_boxplot() +
  labs(title = "Highway mpg (hwy) vs. class")

# Sina + boxes
mpg %>% 
  ggplot(aes(class, hwy)) +
  ggforce::geom_sina() +
  geom_boxplot(
    notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE,  width=0.3, color="black",  fill="transparent", size=0.75
  ) +
  labs(title = "Highway mpg (hwy) vs. class")

# add color
mpg %>% 
  ggplot(aes(class, hwy, color = class)) +
  ggforce::geom_sina() +
  geom_boxplot(
    notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE,  width=0.3, color="black",  fill="transparent", size=0.75
  ) +
  labs(title = "Highway mpg (hwy) vs. class")

# arrange by means
mpg %>% 
  mutate(mean = mean(hwy), .by = class) %>% 
  arrange(mean) %>%
  mutate(class = fct_inorder(class)) %>% # control plotting order
  ggplot(aes(class, hwy, color = class)) +
  ggforce::geom_sina() +
  geom_boxplot(
    notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE,  width=0.3, color="black",  fill="transparent", size=0.75
  ) +
  labs(title = "Highway mpg (hwy) vs. class")

# save the plot
ggsave(filename = here("plots", paste0(out_file_prefix, "mpg_sina_plot", ".pdf")), device = cairo_pdf, width = 10, height = 5, units = "in")




# 5 Defining functions ----
# think about moving these functions to a 'helper_functions.R' script (see https://github.com/DS3-Course/DS3_2026/Rproject_template for example)
# function_name <- function(arg1, arg2) { use arg inputs to do something here }

# basic scatter plot function
scatter_plot <- function(tbl) {
  tbl |> 
    ggplot(aes(displ, hwy, color = class, shape = drv)) +
    geom_point() +
    theme(
      aspect.ratio = 1,
      panel.border=element_rect(colour="black", fill="transparent"),
      plot.title=element_text(face="bold", hjust=0),
      axis.text=element_text(color="black", size=14),
      axis.text.x=element_text(angle=0, hjust=0.5),
      axis.ticks = element_line(color = "black"), # make sure tick marks are black
      panel.background=element_blank(),
      panel.grid=element_blank(),
      plot.background=element_blank(),
      strip.background = element_blank(), # facet label borders
      legend.key=element_blank(), legend.background=element_blank() # remove grey bg from legend
    )
}

# use the function to make the plot
scatter_plot(mpg)
mpg |> scatter_plot()


# improved scatter plot function
scatter_plot2 <- function(tbl, x_var, y_var, color_var, shape_var) {
  tbl |> 
    ggplot(aes(
     x = !!enquo(x_var),
     y = !!enquo(y_var),,
     color = !!enquo(color_var), 
     shape = !!enquo(shape_var)
    )) +
    geom_point() +
    theme(
      aspect.ratio = 1,
      panel.border=element_rect(colour="black", fill="transparent"),
      plot.title=element_text(face="bold", hjust=0),
      axis.text=element_text(color="black", size=14),
      axis.text.x=element_text(angle=0, hjust=0.5),
      axis.ticks = element_line(color = "black"), # make sure tick marks are black
      panel.background=element_blank(),
      panel.grid=element_blank(),
      plot.background=element_blank(),
      strip.background = element_blank(), # facet label borders
      legend.key=element_blank(), legend.background=element_blank() # remove grey bg from legend
    )
}

# use the function to make the plot
mpg |> 
  scatter_plot2(
    x_var = displ,
    y_var = hwy,
    color_var = class,
    shape_var = drv
  ) +
  labs(title = "Highway mpg (hwy) vs. displacement")

# change the variables
mpg |> 
  scatter_plot2(
    x_var = hwy,
    y_var = displ,
    color_var = class,
    shape_var = drv
  ) +
  labs(title = "Displacement vs Highway mpg (hwy)")

