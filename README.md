The **MINLAM** package provides the *fuss_PARALLEL()* function to address Bayesian probability estimation for categorical multimodal data depending on a prior density estimation and an assumed minimum tri-modality. Main functionality of **MINLAM** is to provide subpopulation detection and data probability assignment to these subpopulations. The core of **MINLAM** is the Metrolpolitan-Hastings sampling algorithm, written by Virgilio Gómez-Rubio, which is built upon the **INLA** (intergated nested Laplace approximation) package.

### Useful Links
Mixture models by Virgilio Gómez-Rubio: 
https://becarioprecario.bitbucket.io/inla-gitbook/ch-mixture.html

Source code for the Metropolitan-Hastings sampling:
https://rdrr.io/rforge/INLABMA/src/R/INLAMH.R

INLA homepage: 
https://www.r-inla.org/

### Dependencies
INLA, tidyverse, multimode

### Installation

```r
devtools::install_github("DijoG/MINLAM")
```
### Example
For dummy data creation the **truncnorm** package is needed.

```r
require(tidyverse);require(truncnorm)

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

```r
# Check available cores and wrangle data accordingly
parallelly::availableCores() 

cores <- parallelly::availableCores() - 1   # cores = 3

num_classes <- length(unique(df$Category))
num_groups <- ceiling(num_classes / cores)

df <- 
  df %>%
  mutate(GROUP = as.numeric(factor(Category, levels = unique(Category))) %% num_groups + 1)

df_GROUPS <- 
  df %>%
  group_split(GROUP)

# Run fuss_PARALLEL() with parameters: 'within' = 1 and 'df_prob' = FALSE
dir.create(".../test_wi1")

require(furrr)

tictoc::tic()
MINLAM::fuss_PARALLEL(data = df_GROUPS,
                      varCLASS = "Category", 
                      varY = "Value", 
                      method = "dpi", 
                      within = 1, 
                      maxNGROUP = 5, 
                      df_prob = FALSE, 
                      out_dir = ".../test_wi1", 
                      n_workers = cores)
tictoc::toc()
```
57 minutes processing time using 3 cores.
The *test_wi1* output directory contains the *weighted* csv files.

<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/MM_03.png">

A weighted, for example *wdf_F.csv* csv file has the following information.

<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/MM_04.png">

```r
# Run fuss_PARALLEL() with parameters: 'within' = 0.5 and 'df_prob' = TRUE
dir.create(".../test_wi05")

tictoc::tic()
MINLAM::fuss_PARALLEL(data = df_GROUPS,
                      varCLASS = "Category", 
                      varY = "Value", 
                      method = "dpi", 
                      within = 0.5, 
                      maxNGROUP = 5, 
                      df_prob = TRUE, 
                      out_dir = ".../test_wi05", 
                      n_workers = cores)
tictoc::toc()
```
59 minutes processing time using 3 cores.
The *test_wi05* output directory contains the *weighted* as well as the *data* csv files.

<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/MM_05.png">

A data, for example *df_F.csv* csv file has the following information (not all records shown).

<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/MM_06.png">

Validation by matching assigned (predicted) group labels to original subgroup labels.

```r
# Read data csv files (there are 3 subgroups predicted in all 9 categories)
FIL <- list.files(".../test_wi05", pattern = "^df_", full.names = TRUE) 

# Predicted subgroups
V <- 
  map_dfr(FIL, ~ read_csv2(.x, show_col_types = FALSE) %>%
    as.data.frame() %>%
    mutate(Main_Class = factor(as.character(Main_Class)))) %>%
    arrange(y)
V$Main_Class <- fct_recode(V$Main_Class, "F" = "FALSE")

# Observed subgroups
df <-
  df %>%
  arrange(Value)
df$Subpopulation <- as.numeric(str_remove(df$Subpopulation, "Group "))

# Compute MATCHING percentages (each subgroup has 75 records)
matching_indices <- which(df$Subpopulation == V$Assigned_Group)
main_class_percent <- table(V$Main_Class[matching_indices]) / 75 * 100

# Put percentages into a data frame
label_data <- 
  data.frame(
  Main_Class = names(main_class_percent),
  Percent = format(round(as.numeric(main_class_percent), 1), nsmall = 1))

# Plot 03 ~ validation
V %>%
  ggplot(aes(x = y)) +
  geom_density(col = NA, fill = "grey98", adjust = 0.8) +
  geom_jitter(aes(y = 0.05, color = factor(Assigned_Group)), height = 0.05, 
              size = 2, shape = 16, alpha = .5) + 
  scale_color_manual(values = c("firebrick2", "forestgreen", "cyan3"), 
                     name = "Assigned Groups") +  
  theme_dark() +
  labs(title = "Multimodal Data ~ Validation", 
       x = "Value", y = "Density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  facet_wrap(~ Main_Class, ncol = 3) +  
  geom_text(data = label_data, aes(x = Inf, y = Inf, 
                                   label = str_c(Percent, "%")), 
            hjust = 1.2, vjust = 1.2, size = 5, fontface = "bold", 
            inherit.aes = FALSE, col = "grey15") +  
  theme(legend.position = "top",
        legend.key = element_rect(fill = "transparent", color = NA),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5)) +
  guides(color = guide_legend(override.aes = list(alpha = .7)))
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/MM_07.png">

Not bad at all!:)
