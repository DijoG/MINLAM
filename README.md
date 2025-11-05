# MINLAM 🏔️

<img src="https://latex.codecogs.com/svg.latex?\mathcal{N}(\mu,\sigma^2)" title="Normal Distribution" align="right" height="40">

**MINLAM** is an R package for Bayesian probability estimation in categorical multimodal data. It performs subpopulation detection and probability assignment using Metropolitan-Hastings sampling built upon the `INLA` (Integrated Nested Laplace Approximation) framework.

## Overview

The package's core function `fuss_PARALLEL()`:

  - Uses the Metrolpolitan-Hastings sampling algorithm written by Virgilio Gómez-Rubio

  - Enables Bayesian probability estimation for multimodal categorical data

  - Performs subpopulation detection and classification

  - Supports parallel processing for efficient computation

  - Provides probability assignment to identified subpopulations.

## Installation
```r
devtools::install_github("DijoG/MINLAM")
```
## Dependencies

`INLA`, `tidyverse`, `multimode`, `furrr`

## Quick Start

### Basic Usage
```r
require(MINLAM)

# Prepare data and run analysis
result <- fuss_PARALLEL(
  data = df_GROUPS,
  varCLASS = "Category", 
  varY = "Value", 
  method = "dpi", 
  within = 1, 
  maxNGROUP = 5, 
  df_prob = FALSE, 
  out_dir = ".../output", 
  n_workers = parallelly::availableCores()
)
```

## Detailed Example

### Data Generation

For dummy data creation the **truncnorm** package is needed.

```r
library(tidyverse);library(truncnorm)

# Set seed for reproducibility
set.seed(5)

# Define nine categories with three subpopulations
categories <- rep(LETTERS[1:9], each = 75)
subpopulations <- rep(rep(c("Group 1", "Group 2", "Group 3"), each = 25), times = 9)

# Generate data with single-peaked distributions within each subgroup
values <- c(
  rtruncnorm(25, a = 5, b = 10, mean = 6, sd = 0.4), rtruncnorm(25, a = 5, b = 10, mean = 7.5, sd = 0.4), rtruncnorm(25, a = 5, b = 10, mean = 9, sd = 0.4),
  rtruncnorm(25, a = 5, b = 10, mean = 6.2, sd = 0.5), rtruncnorm(25, a = 5, b = 10, mean = 7.7, sd = 0.5), rtruncnorm(25, a = 5, b = 10, mean = 9.2, sd = 0.5),
  rtruncnorm(25, a = 5, b = 10, mean = 5.8, sd = 0.6), rtruncnorm(25, a = 5, b = 10, mean = 7.4, sd = 0.6), rtruncnorm(25, a = 5, b = 10, mean = 8.9, sd = 0.6),
  rtruncnorm(25, a = 5, b = 10, mean = 6.1, sd = 0.4), rtruncnorm(25, a = 5, b = 10, mean = 7.8, sd = 0.4), rtruncnorm(25, a = 5, b = 10, mean = 9.3, sd = 0.4),
  rtruncnorm(25, a = 5, b = 10, mean = 6.3, sd = 0.5), rtruncnorm(25, a = 5, b = 10, mean = 7.8, sd = 0.5), rtruncnorm(25, a = 5, b = 10, mean = 9.4, sd = 0.5),
  rtruncnorm(25, a = 5, b = 10, mean = 5.9, sd = 0.6), rtruncnorm(25, a = 5, b = 10, mean = 7.5, sd = 0.6), rtruncnorm(25, a = 5, b = 10, mean = 9.2, sd = 0.6),
  rtruncnorm(25, a = 5, b = 10, mean = 6.4, sd = 0.4), rtruncnorm(25, a = 5, b = 10, mean = 7.9, sd = 0.4), rtruncnorm(25, a = 5, b = 10, mean = 9.5, sd = 0.4),
  rtruncnorm(25, a = 5, b = 10, mean = 6.0, sd = 0.5), rtruncnorm(25, a = 5, b = 10, mean = 7.6, sd = 0.5), rtruncnorm(25, a = 5, b = 10, mean = 9.3, sd = 0.5),
  rtruncnorm(25, a = 5, b = 10, mean = 6.2, sd = 0.6), rtruncnorm(25, a = 5, b = 10, mean = 7.8, sd = 0.6), rtruncnorm(25, a = 5, b = 10, mean = 9.6, sd = 0.6)
)

# Create data frame
df <- data.frame(Category = categories, Subpopulation = subpopulations, Value = values)
```
### Data Visualization
```r
# Plot 01 ~ subpopulations/subgroups not shown
ggplot(df, aes(x = Value)) +
  geom_density(color = NA, fill = "grey98", adjust = .8) +
  facet_wrap(~Category) +
  theme_dark() +
  labs(title = "Multimodal Data ~ Density", 
       x = "Value", y = "Density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme(legend.position = "top",
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5))
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/MM_01.png">

```r
# Plot 02 ~ subgroups shown
ggplot(df, aes(x = Value, fill = Subpopulation)) +
  geom_density(alpha = 0.5, color = NA) +
  scale_fill_manual(values = c("firebrick2", "forestgreen", "cyan3"), 
                     name = "Subgroups") +
  facet_wrap(~Category) +
  theme_dark() +
  labs(title = "Multimodal Data ~ Density with Subgroups", 
       x = "Value", y = "Density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme(legend.position = "top",
        legend.key = element_rect(fill = "transparent", color = NA),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5)) +
  guides(fill = guide_legend(override.aes = list(alpha = .6)))
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/MM_02.png">


### Parallel Processing Setup 
```r
# Configure parallel processing
cores <- length(unique(df$Subpopulation))   # cores = 3
num_classes <- length(unique(df$Category))
num_groups <- ceiling(num_classes / cores)

df <- 
  df %>%
  mutate(GROUP = as.numeric(factor(Category, levels = unique(Category))) %% num_groups + 1)

df_GROUPS <- 
  df %>%
  group_split(GROUP)
```
### Running Analysis
```r
library(furrr)

# Analysis with probability output
tictoc::tic()
MINLAM::fuss_PARALLEL(data = df_GROUPS,
                      varCLASS = "Category", 
                      varY = "Value", 
                      method = "dpi", 
                      within = 1, 
                      maxNGROUP = 5, 
                      df_prob = FALSE, 
                      out_dir = ".../test_MINLAM", 
                      n_workers = cores)
tictoc::toc()
# Processing time: ~16 minutes (3 cores)
```
### Output 

The function generates:
  - **Data CSV** files: Original data with assigned groups and probabities
  - **Weighted CSV** files: Probabilty weights for each subpopulation.
  
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/MM_05.png">

A **Data CSV** file consists of the following fields:
  - `y`: Original/observed value
  - `Group`: Original/observed subgroup
  - `Group_1`: Predicted belonging probability to
  - `Group_2`: Predicted belonging probability to
  - `Group_3`: Predicted belonging probability to
  - `Assigned_Group`: Assigned/predicted subgroup
  - `Min_Assigned`: Minimum value of the assigned/predicted range
  - `Max_Assigned`: Maximum value of the assigned/predicted range
  - `Mean_Assigned`: Mean value of the assigned/predicted range
  - `Mode_Assigned`: Mode of the assigned/predicted range
  - `Main_Class`: Category/class 

### Validation
```r
# Validate subgroup assignments
FIL <- list.files(".../test_MINLAM", pattern = "^df_", full.names = TRUE) 

predicted <- map_dfr(FIL, ~ read_csv2(.x, show_col_types = FALSE)) %>%
  as.data.frame() %>%
  mutate(Main_Class = factor(as.character(Main_Class))) %>%
  arrange(y)

predicted$Main_Class <- fct_recode(predicted$Main_Class, "F" = "FALSE")

# Prepare observed data
df <- df %>% arrange(Value)
df$Subpopulation <- as.numeric(str_remove(df$Subpopulation, "Group "))

# Calculate matching accuracy (each subgroup has 75 records)
matching_indices <- which(df$Subpopulation == predicted$Assigned_Group)
main_class_percent <- table(predicted$Main_Class[matching_indices]) / 75 * 100

label_data <- data.frame(
  Main_Class = names(main_class_percent),
  Percent = format(round(as.numeric(main_class_percent), 1), nsmall = 1)
)

# Plot 03 - validation
ggplot(predicted, aes(x = y)) +
  geom_density(col = NA, fill = "grey98", adjust = 0.8) +
  geom_jitter(aes(y = 0.05, color = factor(Assigned_Group)), 
              height = 0.05, size = 2, shape = 16, alpha = .5) + 
  scale_color_manual(values = c("firebrick2", "forestgreen", "cyan3"), 
                     name = "Assigned Groups") +  
  theme_dark() +
  labs(title = "Validation of Subgroup Assignments", 
       x = "Value", y = "Density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  facet_wrap(~ Main_Class, ncol = 3) +  
  geom_text(data = label_data, aes(x = Inf, y = Inf, 
            label = paste0(Percent, "%")), 
            hjust = 1.2, vjust = 1.2, size = 5, fontface = "bold", 
            inherit.aes = FALSE, col = "grey15") +  
  theme(legend.position = "top",
        legend.key = element_rect(fill = "transparent", color = NA),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5)) +
  guides(color = guide_legend(override.aes = list(alpha = .7)))
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/MM_07.png">

**Validation results show accurate subgroup assignment across categories.**

## Useful Links

  - Mixture models by Virgilio Gómez-Rubio: 
    https://becarioprecario.bitbucket.io/inla-gitbook/ch-mixture.html

  - Algorithm for the Metropolitan-Hastings sampling:
    https://rdrr.io/rforge/INLABMA/src/R/INLAMH.R

  - INLA homepage: 
    https://www.r-inla.org/


**Happy multimoda(e)ling!** 🏔️