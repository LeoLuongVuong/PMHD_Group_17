# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(table1)
library(skimr)
# for PCA
library('corrr')
library(ggcorrplot)
library("FactoMineR")
library("factoextra")


# Question a: binary outcome for non-gaussian data --------------------------

# Import the data
setwd("./final_report/datasets")
non_gaussian_data <- read.csv("count_data_G17.csv")
gaussian_data <- read.csv("gaussian_data_G17.csv")

# Note

# “bushID”: Index of the bush from which the rose was cut - what is this?
# Don't see in the dataset?

# “T_0”,”T_1”,”…”: Width (cm) of the flower on day 0, 1, …



## Some EDA ---------------------------------------------------------------

str(non_gaussian_data)
summary(non_gaussian_data)

# Convert all variables to factor except tot.vase.days

non_gaussian_data$compound <- as.factor(non_gaussian_data$compound)
non_gaussian_data$garden <- as.factor(non_gaussian_data$garden)
non_gaussian_data$species <- as.factor(non_gaussian_data$species)
non_gaussian_data$subplotID <- as.factor(non_gaussian_data$subplotID)
non_gaussian_data$flowerID <- as.factor(non_gaussian_data$flowerID)
non_gaussian_data$rater <- as.factor(non_gaussian_data$rater)

# Group data by subplot and species, and count the number of flowers

strata <- non_gaussian_data %>%
  group_by(subplotID, species, compound) %>%
  summarise(n = n()) 
View(strata)

# Check crossed design - https://lme4.r-forge.r-project.org/book/Ch2.pdf

# Crossed means we have at least one observation for each combination of a level
# of each fator
xtabs(~ species + rater + subplotID, non_gaussian_data)
# unbalanced study design

# Question c: Explore/Visualize T0 up to T20 with a multivariate method -------

# First, convert gaussian_data to long format

gaussian_long <- gaussian_data |> 
  pivot_longer(
    cols = !(Flower_index:Subplot), 
    names_to = c(".value", "Day"), 
    names_sep = "_", 
    values_drop_na = TRUE
  )

# Rename column "T" to "Width"
gaussian_long <- gaussian_long |> 
  rename(Width = T)

View(gaussian_long)

## Some EDA ---------------------------------------------------------------

# Get some insight into the dataset structure
skim(gaussian_long)
summary(gaussian_long)
skim(gaussian_data) # some data missing, but not substantial -> complete case analysis
summary(gaussian_data)

# Convert columns into factors
for (i in 1:6) {
  gaussian_long[[i]] <- as.factor(gaussian_long[[i]])
}

gaussian_long[["Day"]] <- as.numeric(gaussian_long[["Day"]])

for (i in 1:6) {
  gaussian_data[[i]] <- as.factor(gaussian_data[[i]])
}

# Use table 1 to get an idea about the balance of the data
table1(~ Width | Subplot * Type, data = gaussian_long)

table1(~ Compound | Subplot * Type, data = gaussian_data) # very balance design

# Reorder gaussian_long by Flower_index and Day
gaussian_long <- gaussian_long |> 
  arrange(Flower_index, Day)

# Plot the Width of the flowers over Day, color by Flower_index, facet by Compound
EDA_c <- ggplot(gaussian_long, aes(x = Day, y = Width, group = Flower_index)) +
  geom_line() +
  facet_wrap(~ Compound, ncol = 3, labeller = label_both) +
  theme_minimal() +
  ylab("Width (cm)") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "sans"), 
        strip.text.x = element_text(hjust = 0.5, size = 8, family = "sans"),
        axis.title = element_text(size = 7, family = "sans"),
        axis.text = element_text(size = 7, family = "sans"),
        legend.position = "none",
        panel.grid.minor = element_blank())

ggsave("EDA_c.svg", plot = EDA_c, width = 19, height = 9, units = "cm")

## PCA ---------------------------------------------------------------

# checking for null values
colSums(is.na(gaussian_data))

# creates new data with only T_0 to T_20 columns
T_0_20 <- gaussian_data |> 
  select(T_0:T_20) |> 
  na.omit()
# note: removed all NA observations

# normalise the data for PCA
T_0_20_norm <- scale(T_0_20)
head(T_0_20_norm)

# Compute the correlation matrix
T_0_20_norm_corr_matrix <- cor(T_0_20_norm)
ggcorrplot(T_0_20_norm_corr_matrix)
# they are highly correlated

# Applying PCA
# the princomp() computes the PCA, and summary() function shows the result.
PCA_T_0_20 <- princomp(T_0_20_norm_corr_matrix)
summary(PCA_T_0_20)
# Comp.1 and Comp.2 explains about 98% of the total variance

# explore how the first 2 variables relate to each column
PCA_T_0_20$loadings[, 1:2]

# Scree Plot
scree_plot <- fviz_eig(PCA_T_0_20, addlabels = TRUE)
biplot <- fviz_pca_var(PCA_T_0_20, col.var = "black")

# save the plots
ggsave("scree_plot.svg", plot = scree_plot, width = 19, height = 12.5, units = "cm")
ggsave("biplot.svg", plot = biplot, width = 19, height = 12.5, units = "cm")



