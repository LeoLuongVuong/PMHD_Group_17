# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(table1)
library(skimr)
library(geepack) # for gee
library(MuMIn) # for gee
library(ggeffects) # for ggemmeans()
library(contrast) # for contrast - doesn't work in the first gee model
library(multcomp) # contrasting with glht
library(lme4) # for glmer
library(flexplot) # for model.comparison
library(AICcmodavg) # for AICc()
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

## Some EDA and data manipulations -----------------------------------------

str(non_gaussian_data)
summary(non_gaussian_data)

# Convert all variables to factor except tot.vase.days

non_gaussian_data$compound <- as.factor(non_gaussian_data$compound)
non_gaussian_data$garden <- as.factor(non_gaussian_data$garden)
non_gaussian_data$species <- as.factor(non_gaussian_data$species)
non_gaussian_data$subplotID <- as.factor(non_gaussian_data$subplotID)
non_gaussian_data$flowerID <- as.factor(non_gaussian_data$flowerID)
non_gaussian_data$rater <- as.factor(non_gaussian_data$rater)

### Some quick check ------------------------------------------------

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

### Box plot per compound ------------------------------------------------

Box_plot_total <- non_gaussian_data %>%
  mutate(type = ifelse(compound %in% c(6, 14),"Highlighted","Normal")) %>%
  ggplot(aes(x = factor(compound), y = tot.vase.days, fill = type, alpha = type)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#440154", "grey")) +
  scale_alpha_manual(values = c(0.8, 0.1)) +
  theme_minimal() +
  xlab("Compound") + 
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, by = 5), expand = c(0,0)) +
  ylab("Total vase days") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5)

ggsave("Box_plot_total.jpeg", Box_plot_total, width = 13, height = 7, dpi = 300, unit = "cm")

### Add binary response ------------------------------------------------

# binary_outcome equals non_gaussian_data with no NA

binary_outcome <- non_gaussian_data %>%
  na.omit() 

# replicate each observation of each flowerID 30 times

binary_outcome <- binary_outcome[rep(row.names(binary_outcome), 30),]

# sort by flowerID, remove rownames

binary_outcome <- binary_outcome[order(binary_outcome$flowerID),]

rownames(binary_outcome) <- NULL

# add day column, which goes from 1 to 30 per each flowerID

binary_outcome$day <- rep(1:30, nrow(binary_outcome)/30)

# add fresh column, which is 1 on day that is <= tot.vase.days, 0 otherwise

binary_outcome$fresh <- ifelse(binary_outcome$day <= binary_outcome$tot.vase.days, 1, 0)

## GEE ---------------------------------------------------------------

### only clustered within flowerID --------------------------------------------

## Top-down model selection strategy

# full model
gee_full <- geeglm(fresh ~ compound*day + species + subplotID + rater, 
                    data = binary_outcome, 
                    id = flowerID, 
                    family = binomial, 
                    corstr = "exchangeable")
summary(gee_full)
anova(gee_full)

geepack::QIC(gee_full) #36034     

# remove subplotID
gee_no_subplotID <- geeglm(fresh ~ compound*day + species + rater, 
                           data = binary_outcome, 
                           id = flowerID, 
                           family = binomial, 
                           corstr = "exchangeable")
summary(gee_no_subplotID)
anova(gee_no_subplotID)

geepack::QIC(gee_no_subplotID) #40663     

# remove rater
gee_no_rater <- geeglm(fresh ~ compound*day + species + subplotID, 
                       data = binary_outcome, 
                       id = flowerID, 
                       family = binomial, 
                       corstr = "exchangeable")
summary(gee_no_rater)
anova(gee_no_rater)

geepack::QIC(gee_no_rater) #47884

# remove species
gee_no_species <- geeglm(fresh ~ compound*day + rater + subplotID, 
                         data = binary_outcome, 
                         id = flowerID, 
                         family = binomial, 
                         corstr = "exchangeable")
summary(gee_no_species)
anova(gee_no_species)

geepack::QIC(gee_no_species) #36147

#### model selection with QIC --------------------------------------------

QIC <- MuMIn::QIC
## Not run: 
QIC <- function(x) geepack::QIC(x)[1]

model.sel(gee_full, gee_no_subplotID, gee_no_rater, gee_no_species, rank = QIC)

# gee_full is the best model

#### conclusion -----------------------------------------------------------

# compound 14 & 15 have a better time effect compared to compound 1. However,
# the effect is not significant.

# create effect plot
plot(ggemmeans(gee_full, terms = c("day", "compound"),
               conditions = c(species = 2, rater = 2, subplotID = 3))) + 
  ggplot2::ggtitle("GEE Effect plot")
# doesn't work

## Contrasting with contrast function
# contrast compound14:day vs compound15:day

print(
  contrast(
    gee_full, 
    list(compound = "14"),
    list(compound = "15")
  ),
  X = TRUE)

## Contrasting with glht function

# difference between compound14:day and compound15:day of gee_full
K <- matrix(c(rep(0,55), 1, -1), 1)
t <- glht(gee_full, linfct = K)
summary(t)

#compound14:day and compound15:day are not significantly different, and are 
# not significantly better than water

## GLMM for binary outcome --------------------------------------------------

glmm_full <- glmer(fresh ~ compound*day + rater + species + garden + (1|flowerID) + (1|subplotID),
                          family = binomial(link = "logit"),
                          data = binary_outcome, nAGQ = 0)
summary(glmm_full)

AICc(glmm_full) # 15394

# remove garden

glmm_no_garden <- glmer(fresh ~ compound*day + rater + species + (1|flowerID) + (1|subplotID),
                        family = binomial(link = "logit"),
                        data = binary_outcome, nAGQ = 0)
summary(glmm_no_garden)

AICc(glmm_no_garden) # 15392 # can be removed

# remove rater

glmm_no_rater <- glmer(fresh ~ compound*day + species + (1|flowerID) + (1|subplotID),
                      family = binomial(link = "logit"),
                      data = binary_outcome, nAGQ = 0)
summary(glmm_no_rater)

AICc(glmm_no_rater) # 17016 # can't be removed

# remove species

glmm_no_species <- glmer(fresh ~ compound*day + rater + (1|flowerID) + (1|subplotID),
                        family = binomial(link = "logit"),
                        data = binary_outcome, nAGQ = 0)
summary(glmm_no_species)

AICc(glmm_no_species) # 15416 # can't be removed

# glmm_no_garden with random slope for compound:day

glmm_no_garden_slope <- glmer(fresh ~ compound + day + compound:day + rater + species + (1 + compound:day|flowerID) + (1|subplotID),
                              family = binomial(link = "logit"),
                              data = binary_outcome, nAGQ = 0)
summary(glmm_no_garden_slope)

AICc(glmm_no_garden_slope) # 15524 # shouldn't be added

# Conclusion: glmm_no_garden is the best model

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



