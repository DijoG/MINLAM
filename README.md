The **MINLAM** package provides the *fuss_PARALLEL()* function to address Bayesian probabilty estimation for categorical multimodal 
data depending on a prior density estimation. Main functionality of **MINLAM** is to provide subpopulation detection and data probabilty
belonging to these subpopulations. It is based on the Metrolpolitan-Hastings sampling written by Virgilio Gómez-Rubio for internal integration 
within the **INLA** package.

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
  geom_density(alpha = 0.7, color = NA, fill = "grey98", adjust = .8) +
  facet_wrap(~Category, scales = "free_y") +
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
  geom_density(alpha = 0.7, color = NA) +
  facet_wrap(~Category, scales = "free_y") +
  theme_dark() +
  labs(title = "Multimodal Data ~ Density with Subgroups", 
       x = "Value", y = "Density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme(legend.position = "top",
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5))
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/MM_02.png">

```r
# Check available cores and wrangle data accordingly
parallelly::availableCores() 

cores <- parallelly::availableCores() - 3   # cores = 3

num_classes <- length(unique(df$Category))
num_groups <- ceiling(num_classes / cores)

df <- 
  df %>%
  mutate(GROUP = as.numeric(factor(Category, levels = unique(Category))) %% num_groups + 1)

df_GROUPS <- 
  df %>%
  group_split(GROUP)

# Run fuss_PARALLEL() with parameter 'within' = 1
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
# Run fuss_PARALLEL() with parameter 'within' = 0.5
dir.create(".../test_wi05")

tictoc::tic()
MINLAM::fuss_PARALLEL(data = df_GROUPS,
                      varCLASS = "Category", 
                      varY = "Value", 
                      method = "dpi", 
                      within = 0.5, 
                      maxNGROUP = 5, 
                      df_prob = FALSE, 
                      out_dir = ".../test_wi05", 
                      n_workers = cores)
tictoc::toc()
```
59 minutes processing time using 3 cores.
The *test_wi05* output directory contains the *weighted* as well as the *data* csv files.

<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/MM_05.png">

A data, for example *df_F.csv* csv file has the following information.

<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/README/MM_06.png">
