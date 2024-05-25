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
library('corrr') # for PCA
library(ggcorrplot) # for PCA
library("FactoMineR") # for PCA
library("factoextra") # for PCA
library(car) # for Anova()
library(gt) # for making tables
library(webshot2)
library(broom)
library(broom.mixed)
library(xtable) # generate output for latex
library(varTestnlme) # check the variance component

# Question a: binary outcome for non-gaussian data --------------------------

# Import the data
setwd("./final_report/data_analyses")
non_gaussian_data <- read.csv("count_data_G17.csv")

## Some EDA and data manipulations -----------------------------------------

str(non_gaussian_data)
summary(non_gaussian_data)

# check percentage of missing tot.vase.days

sum(is.na(non_gaussian_data$tot.vase.days))/nrow(non_gaussian_data)*100
# 5.5% missing -> important for discussion!

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
  scale_x_discrete(labels = c("Distilled Water",
                              "Apathic Acid",
                              "Beerse Brew",
                              "Concentrate of Caducues",
                              "Distilled of Discovery",
                              "Essence of Epiphanea",
                              "Four in December",
                              "Granule of Geheref",
                              "Kar-Hamel Mooh",
                              "Lucifer’s Liquid",
                              "Noospherol",
                              "Oil of John’s Son",
                              "Powder of Perlimpinpin",
                              "Spirit of Scienza",
                              "Zest of Zen")) +
  ylab("Total vase days") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold"),
        axis.title = element_text(size = 11), #adjust text size 240524
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis text
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.5)
Box_plot_total

# export the plot
setwd("./eda_plots")
ggsave("Box_plot_total.jpeg", Box_plot_total, width = 19, height = 9, dpi = 300, unit = "cm")
# return to the previous directory
Path <- getwd()
setwd(dirname(Path))

### Add binary response ------------------------------------------------

# binary_outcome equals non_gaussian_data with no NA

binary_outcome <- non_gaussian_data %>%
  na.omit() # removes 165 observations

# replicate each observation of each flowerID 30 times

binary_outcome <- binary_outcome[rep(row.names(binary_outcome), 30),]

# sort by flowerID, remove rownames

binary_outcome <- binary_outcome[order(binary_outcome$flowerID),]

rownames(binary_outcome) <- NULL

# add day column, which goes from 1 to 30 per each flowerID

binary_outcome$day <- rep(1:30, nrow(binary_outcome)/30)

# add fresh column, which is 1 on day that is <= tot.vase.days, 0 otherwise

binary_outcome$fresh <- ifelse(binary_outcome$day <= binary_outcome$tot.vase.days, 1, 0)

### Visualize Freshness by Compound ---------------------------------------

sum_binary_outcome <- binary_outcome %>% 
  group_by(day, compound) %>% 
  summarise(fresh_sum = sum(fresh))

# create a percentage column, which is the percentage of fresh flowers per day
# at each compound, being 100% per each compound on day 1, this percentage decreases
# from day 2 on ward (e.g. percentage day 2 = fresh_sum day 2 / fresh_sum day 1 * 100)

sum_binary_outcome <- sum_binary_outcome %>% 
  group_by(compound) %>% 
  mutate(percentage = fresh_sum/max(fresh_sum)*100)

# visual percentage per day
compound_labels = c("Distilled Water",
             "Apathic Acid",
             "Beerse Brew",
             "Concentrate of Caducues",
             "Distilled of Discovery",
             "Essence of Epiphanea",
             "Four in December",
             "Granule of Geheref",
             "Kar-Hamel Mooh",
             "Lucifer’s Liquid",
             "Noospherol",
             "Oil of John’s Son",
             "Powder of Perlimpinpin",
             "Spirit of Scienza",
             "Zest of Zen")

sum_binary <- sum_binary_outcome %>%
  ggplot(aes(day, percentage, color = as.factor(compound))) + 
  geom_line(size = 1) +
  theme_minimal() +
  labs(color = "Compound") +
  ylab("Percentage of fresh flower") +
  xlab("Day") +
  scale_color_viridis_d(labels = compound_labels, option = "D") +
  #scale_color_manual() +
  theme(axis.title = element_text(size = 11, family = "sans"),
        axis.text = element_text(size = 10, family = "sans"),
        legend.title = element_text(size = 9, family = "sans"),
        legend.text = element_text(size = 9, family = "sans"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-1,-1,-1,-1))
sum_binary

# export the plot
setwd("./eda_plots")
ggsave("sum_binary.png", plot = sum_binary, width = 19, height = 12, units = "cm")
# go back to the previous wd
Path = getwd()
setwd(dirname(Path))

## GEE ---------------------------------------------------------------
# only clustered within flowerID

### Choose set of covariates --------------------------------------------

## Top-down model selection strategy

#### full model -----------------------------------------------------------
gee_full <- geeglm(fresh ~ compound*day + species + subplotID + rater, 
                    data = binary_outcome, 
                    id = flowerID, 
                    family = binomial, 
                    corstr = "exchangeable")
summary(gee_full)
anova(gee_full)

geepack::QIC(gee_full) #36034  

#### full remove compound ------------------------------------------------
gee_no_compound <- geeglm(fresh ~ day + compound:day + species + subplotID + rater, 
                          data = binary_outcome, 
                          id = flowerID, 
                          family = binomial, 
                          corstr = "exchangeable")

geepack::QIC(gee_no_compound) # 36384 # significantly worse

summary(gee_no_compound) # best models, compound 6 and 14 are significantly 
# better than distilled water
anova(gee_no_compound)

#### no interaction ------------------------------------------------
gee_no_interaction <- geeglm(fresh ~ compound + day + species + subplotID + rater, 
                             data = binary_outcome, 
                             id = flowerID, 
                             family = binomial, 
                             corstr = "exchangeable")

geepack::QIC(gee_no_interaction)  # 39120 # significantly worse

#### only with interaction --------------------------------------------
gee_only_interaction <- geeglm(fresh ~ compound:day + species + subplotID + rater, 
                                                     data = binary_outcome, 
                                                     id = flowerID, 
                                                     family = binomial, 
                                                     corstr = "exchangeable")

geepack::QIC(gee_only_interaction) # 36384 # significantly worse

summary(gee_only_interaction)
anova(gee_only_interaction)

#### no subplotID ------------------------------------------------
gee_no_subplotID <- geeglm(fresh ~ day + compound:day + species + rater, 
                           data = binary_outcome, 
                           id = flowerID, 
                           family = binomial, 
                           corstr = "exchangeable")
summary(gee_no_subplotID)
anova(gee_no_subplotID)

geepack::QIC(gee_no_subplotID) #40891     

#### no rater ------------------------------------------------
gee_no_rater <- geeglm(fresh ~ day + compound:day + species + subplotID, 
                       data = binary_outcome, 
                       id = flowerID, 
                       family = binomial, 
                       corstr = "exchangeable")
summary(gee_no_rater)
anova(gee_no_rater)

geepack::QIC(gee_no_rater) #48000 

#### no species ------------------------------------------------
gee_no_species <- geeglm(fresh ~ day + compound:day + subplotID + rater, 
                         data = binary_outcome, 
                         id = flowerID, 
                         family = binomial, 
                         corstr = "exchangeable")
summary(gee_no_species)
anova(gee_no_species)

geepack::QIC(gee_no_species) #36473

#### model selection with QIC --------------------------------------------

QIC <- MuMIn::QIC
## Not run: 
QIC <- function(x) geepack::QIC(x)[1]

model.sel(gee_full, gee_no_compound, gee_no_interaction, gee_only_interaction, 
          gee_no_subplotID, gee_no_rater, gee_no_species, rank = QIC)

# gee_no_compound is the best model

#### comparing models using anova -----------------------------------------

# gee_no_compound vs gee_no_rater
anova(gee_no_compound, gee_no_rater) # too heavy

# compare simpler model
test_rater <- geeglm(fresh ~ rater, 
                                     data = binary_outcome, 
                                     id = flowerID, 
                                     family = binomial, 
                                     corstr = "exchangeable")
test_no_rater <- geeglm(fresh ~ 1, 
                        data = binary_outcome, 
                        id = flowerID, 
                        family = binomial, 
                        corstr = "exchangeable")
anova(test_rater, test_no_rater) # 0.0001

#### conclusion -----------------------------------------------------------

# compound 6 & 14 are significantly better than water.

### Choose the working correlation --------------------------------------------

#### Independence ------------------------------------------------

gee_no_compound_ind <- geeglm(fresh ~ compound:day + day + species + subplotID + rater, 
                              data = binary_outcome, 
                              id = flowerID, 
                              family = binomial, 
                              corstr = "independence")

geepack::QIC(gee_no_compound_ind) # 36094 #significantly better

summary(gee_no_compound_ind)
anova(gee_no_compound_ind)

#### autoregressive ------------------------------------------------

gee_no_compound_ar <- geeglm(fresh ~ compound:day + day + species + subplotID + rater, 
                             data = binary_outcome, 
                             id = flowerID, 
                             family = binomial, 
                             corstr = "ar1")

geepack::QIC(gee_no_compound_ar) # 36066 #significantly better

summary(gee_no_compound_ar)
anova(gee_no_compound_ar)

#### unstructured ------------------------------------------------

gee_no_compound_un <- geeglm(fresh ~ compound:day + day + species + subplotID + rater, 
                             data = binary_outcome, 
                             id = flowerID, 
                             family = binomial, 
                             corstr = "unstructured")

geepack::QIC(gee_no_compound_un) # 36034 #significantly better

summary(gee_no_compound_un)
anova(gee_no_compound_un)

#### model selection with QIC --------------------------------------------

QIC <- MuMIn::QIC
## Not run: 
QIC <- function(x) geepack::QIC(x)[1]

model.sel(gee_no_compound_ind, gee_no_compound, gee_no_compound_ar, rank = QIC)

#### conclusion -------------------------

# autoregressive is the best working correlation

### Select the best compound -------------------------

#### Contrasting with glht function --------------------------------------------

## Compound 2, 6, 14 are significantly better than water. So we need to contrast
# them

# Extract all the coefficients from gee_no_compound_ar to a matrix
coef_matrix <- diag(length(coef(gee_no_compound_ar)))[-1,]
rownames(coef_matrix) <- names(coef(gee_no_compound_ar))[-1]
coef_matrix

## Finally, compare the three compounds
contrast_matrix <- rbind(
  "compound6:day - compound14:day = 0" = c(rep(0,33), 1, rep(0, 7), -1, 0),
  "compound14:day - compound2:day = 0" = c(rep(0,29), -1, rep(0, 11), 1, 0),
  "compound6:day - compound2:day = 0" = c(rep(0,29), -1, rep(0, 3), 1, rep(0, 9))
)
contrast_test <- glht(gee_no_compound_ar, linfct = contrast_matrix)
summary(contrast_test)

#### conclusion -------------------------

# compound 6 and 14 are not significantly different from each other.

### Output table for model selection --------------------------------

#### QIC table ------------------------------------------------

# first: dataframe
gee_dat <- data.frame(
  model_name = c(1, 2, 3, 4),
  qic = c(36034, 40663, 47884, 36147),
  description = c(
    "Compound*Day + Species + SubplotID + Rater",
    "Compound*Day + Species + Rater",
    "Compound*Day + Species + SubplotID",
    "Compound*Day + Rater + SubplotID"
  )
)

# create table
gee_tab_qic <- gee_dat %>%
  gt() %>%
  tab_header(
    title = "Model Selection Process",
    subtitle = "Comparison of Models based on Quasi-likelihood under the Independence Model Criterion (QIC)"
  ) %>%
  cols_label(
    model_name = "Model",
    qic = "QIC",
    description = "Formula"
  ) %>%
  cols_align(
    align = "center"
  )


setwd("./tables")
gtsave(gee_tab_qic, "tab_gee_qic.html")
webshot("tab_gee_qic.html", "tab_qic.pdf")

#### Tidy table for parameter estimates --------------------------------

tidy_gee <- tidy(gee_no_compound_ar, conf.int = TRUE, exponentiate = TRUE)

# Filter out rows corresponding to compound, subplotID, and rater effects
tidy_gee_filtered <- tidy_gee %>%
  filter(!grepl("^compound[0-9]+$|^subplotID[0-9]+$|^rater[0-9]+$", term))

# Round the estimates, confidence intervals, and standard errors to 2 decimals
tidy_gee_filtered$estimate <- round(tidy_gee_filtered$estimate, 2)
tidy_gee_filtered$conf.low <- round(tidy_gee_filtered$conf.low, 2)
tidy_gee_filtered$conf.high <- round(tidy_gee_filtered$conf.high, 2)
tidy_gee_filtered$std.error <- round(tidy_gee_filtered$std.error, 2)
tidy_gee_filtered$statistic <- round(tidy_gee_filtered$statistic, 2)

# Calculate 1-sided p-value
tidy_gee_filtered$p.value <- tidy_gee_filtered$p.value / 2

# adjust for multiple testing
tidy_gee_filtered$p.value <- ifelse(tidy_gee_filtered$p.value * 44 < 1, tidy_gee_filtered$p.value * 44, 1)

# Round the p-values to 3 decimals
tidy_gee_filtered$p.value <- round(tidy_gee_filtered$p.value, 3)

# Create the summary table excluding the filtered rows
tab_gee <- tidy_gee_filtered %>%
  gt() %>%
  tab_header(
    title = "Table 3: GEE Model Summary",
    subtitle = "Model: fresh ~ compound * day + day + species + subplotID + rater"
  ) %>%
  cols_label(
    term = "Parameters",
    estimate = "OR",
    conf.low = "2.5 % CI",
    conf.high = "97.5 % CI",
    std.error = "Std. Error",
    statistic = "Wald test statistic",
    p.value = html("P-value<sup>")
  ) %>%
  fmt_number(
    columns = vars(estimate, conf.low, conf.high, std.error),
    decimals = 2
  ) %>%
  fmt_number(
    columns = vars(p.value),
    decimals = 3
  ) %>%
  cols_align(
    align = "center",
    columns = vars(term, estimate, conf.low, conf.high, std.error, statistic, p.value)
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "lightblue")
    ),
    locations = cells_body(
      columns = vars(term),
      rows = estimate > 1 & p.value < 0.05
    )
  ) %>%
tab_footnote(
    footnote = "1-sided p-value bonferroni corrected for 44 parameters",
    locations = cells_column_labels(vars(p.value))
  )
tab_gee

setwd("./tables")
gtsave(tab_gee, "tab_gee.html")
gtsave(tab_gee, "tab_gee.pdf")
gtsave(tab_gee, "tab_gee.rtf")
gtsave(tab_gee, "tab_gee.docx")
gtsave(tab_gee, "tab_gee.tex")
webshot("tab_gee.html", "tab_gee.pdf")

## GLMM for binary outcome --------------------------------------------------

### Choose model structure --------------------------------------------

#### full model ------------------------------------------------
glmm_full <- glmer(fresh ~ day + compound:day + species + garden + (1|flowerID) + (1|subplotID) + (1|rater),
                          family = binomial(link = "logit"),
                          data = binary_outcome, nAGQ = 0)
summary(glmm_full)

AICc(glmm_full) # 15447
BIC(glmm_full) # 15644

#### remove garden ------------------------------------------------
glmm_no_garden <- glmer(fresh ~ day + compound:day + species + (1|flowerID) + (1|subplotID) + (1|rater),
                        family = binomial(link = "logit"),
                        data = binary_outcome, nAGQ = 0)

AICc(glmm_no_garden) # 15445 # can be removed
BIC(glmm_no_garden) # 15632 # can be removed

summary(glmm_no_garden)

vcov(glmm_no_garden) # variance-covariance matrix
print(glmm_no_garden, correlation = TRUE) # correlation matrix

#### remove interaction ------------------------------------------------
glmm_no_interaction <- glmer(fresh ~ compound + day + species + (1|flowerID) + (1|subplotID) + (1|rater),
                             family = binomial(link = "logit"),
                             data = binary_outcome, nAGQ = 0)

AICc(glmm_no_interaction) # 15435 # interaction worsens the model
BIC(glmm_no_interaction) # 15617 # interaction worsens the model

summary(glmm_no_interaction)

#### remove rater ------------------------------------------------
glmm_no_rater <- glmer(fresh ~ day + compound:day + species + (1|flowerID) + (1|subplotID),
      family = binomial(link = "logit"),
      data = binary_outcome, nAGQ = 0)

AICc(glmm_no_rater) # 17019 # can't be removed
BIC(glmm_no_rater) # 17197 # can't be removed

summary(glmm_no_rater)

#### remove species ------------------------------------------------
glmm_no_species <- glmer(fresh ~ day + compound:day + (1|flowerID) + (1|subplotID) + (1|rater),
                         family = binomial(link = "logit"),
                         data = binary_outcome, nAGQ = 0)

AICc(glmm_no_species) # 15469 # can't be removed
BIC(glmm_no_species) # 15647 # can't be removed

summary(glmm_no_species)

#### glmm_no_garden with random slope for compound:day -----------------------

control <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

glmm_no_garden_slope <- glmer(fresh ~ day + compound:day + (1 + compound:day|flowerID) + (1|subplotID) + (1|rater),
                              family = binomial(link = "logit"),
                              data = binary_outcome, 
                              nAGQ = 0,
                              control = control)

AICc(glmm_no_garden_slope) # 15524 # shouldn't be added
BIC(glmm_no_garden_slope) # 15712 # shouldn't be added

summary(glmm_no_garden_slope) # this model is way computationally heavy. doesn't work.

### Variance component test & likelihood ratio test -----------------------

#### compare glmm_no_garden with glmm_no_rater -----------------------

vt <- varCompTest(glmm_no_interaction, glmm_no_rater)
# doesn't work for more than one random effect

# compare the model with and without rater
test_glmm_rater <- glmer(fresh ~ day + compound:day + species + (1|rater),
                        family = binomial(link = "logit"),
                        data = binary_outcome, nAGQ = 0)
test_glmm_no_rater <- glm(fresh ~ day + compound:day + species,
                        family = binomial(link = "logit"),
                        data = binary_outcome)
vt <- varCompTest(test_glmm_rater, test_glmm_no_rater)
print(vt)
summary(vt)
# there's a very strong evidence to support the random effect of rater

#### compare the model with and without garden -----------------------

anova(glmm_full, glmm_no_garden) #garden is not significant # p: 0.63

### Select the best compound -------------------------

#### Contrasting with glht function --------------------------------------------

## Compound 6, 14 are significantly better than water. So we need to contrast
# them

# Extract the fixed effects coefficients
fixed_effects <- fixef(glmm_no_garden)

# Create a diagonal matrix for the fixed effects, excluding the intercept
coef_matrix <- diag(length(fixed_effects))[-1,]
rownames(coef_matrix) <- names(fixed_effects)[-1]
coef_matrix

## Finally, compare the three compounds
contrast_matrix <- matrix(c(rep(0,7), 1, rep(0, 7), -1, 0), 1)
contrast_test <- glht(glmm_no_garden, linfct = contrast_matrix)
summary(contrast_test)

#### conclusion -------------------------

# compound 6 and 14 are not significantly different from each other.

### Model output --------------------------------------------------------

#### Model selection ---------------------------------------------------

# Define the model structures
model_structures <- c(
  "day + compound:day + species + garden + (1|flowerID) + (1|subplotID) + (1|rater)",
  "day + compound:day + species + (1|flowerID) + (1|subplotID) + (1|rater)",
  "compound + day + species + (1|flowerID) + (1|subplotID) + (1|rater)",
  "day + compound:day + species + (1|flowerID) + (1|subplotID)",
  "day + compound:day + (1|flowerID) + (1|subplotID) + (1|rater)"
)

# models list
glmm_models <- list(
  glmm_full,
  glmm_no_garden,
  glmm_no_interaction,
  glmm_no_rater,
  glmm_no_species
)

# Calculate AIC and BIC
aic_values <- sapply(glmm_models, AIC)
bic_values <- sapply(glmm_models, BIC)

# Create a data frame
model_df <- data.frame(
  Model = 1:length(glmm_models),
  `Mixed-effects` = model_structures,
  AIC = aic_values,
  BIC = bic_values
)

# Sort by AIC
model_df <- model_df[order(model_df$AIC), ]
model_df$Model <- 1:nrow(model_df)

# Remove the first row, and re-number the models
model_df <- model_df[-1, ]
model_df$Model <- 1:nrow(model_df)

# Create the LaTeX table
print(xtable(model_df, caption = "Model comparison based on AIC and BIC",
             label = "tab:model_comparison"),
      include.rownames = FALSE,
      floating = TRUE,
      tabular.environment = "tabular",
      hline.after = c(-1, 0, nrow(model_df)))

#### Best model parameter estimate --------------------------------------------------

tidy_glmm <- tidy(glmm_no_garden, conf.int = TRUE, exponentiate = TRUE)
View(tidy_glmm)

# start from here
tidy_glmm <- broom.mixed::tidy(glmm_no_garden, conf.int = FALSE, exponentiate = TRUE)

# Remove the first (intercept) row
tidy_glmm <- tidy_glmm[-1, ]

# Round the estimates, confidence intervals, and standard errors to 2 decimals
tidy_glmm$estimate <- round(tidy_glmm$estimate, 2)
tidy_glmm$std.error <- round(tidy_glmm$std.error, 2)
tidy_glmm$statistic <- round(tidy_glmm$statistic, 2)

# # Calculate 1-sided p-value, adjusted for 20 parameters
tidy_glmm$p.value <- round(ifelse(tidy_glmm$p.value / 2 * 20 < 1, tidy_glmm$p.value / 2 * 20, 1), 3)

# change titles of columns 4 to 7 into: OR, Std. Error, Wald test statistic, P-value
colnames(tidy_glmm)[3:7] <- c("Parameters", "Estimate", "Std. Error", "Wald test statistic", "P-value")

# change the name of rows 17 to 19 of column 3 into: Variance(Flower), Variance(Subplot), Variance(Rater)
tidy_glmm[17:19, 3] <- c("Variance(Flower)", "Variance(Subplot)", "Variance(Rater)")

# similarly, change the name of rows 3 to 16 of column 3 into: Apathic Acid : day, 
# Beerse Brew : day, Concentrate of Caducues : day, Distilled of Discovery, 
# Essence of Epiphanea : day, Four in December : day, Granule of Geheref : day, 
# Kar-Hamel Mooh : day, Lucifer’s Liquid : day, Noospherol : day, 
# Oil of John’s Son : day, Powder of Perlimpinpin : day, Spirit of Scienza : day, 
# Zest of Zen : day

tidy_glmm[3:16, 3] <- c("Apathic Acid : day", "Beerse Brew : day", "Concentrate of Caducues : day", 
                         "Distilled of Discovery", "Essence of Epiphanea : day", "Four in December : day", 
                         "Granule of Geheref : day", "Kar-Hamel Mooh : day", "Lucifer’s Liquid : day", 
                         "Noospherol : day", "Oil of John’s Son : day", "Powder of Perlimpinpin : day", 
                         "Spirit of Scienza : day", "Zest of Zen : day")

# Remove the first two columns

tidy_glmm <- tidy_glmm[, -c(1,2)]

# square the rows 17-19 of the 2nd column
tidy_glmm[17:19, 2] <- tidy_glmm[17:19, 2]^2

# export into latex

print(xtable(tidy_glmm, caption = "GLMM model summary",
             label = "tab:GLMM_parameter_estimates"),
      include.rownames = FALSE,
      floating = TRUE,
      tabular.environment = "tabular",
      hline.after = c(-1, 0, nrow(tidy_glmm)))

## Model diagnostics --------------------------------------------------

# Question c: Explore/Visualize T0 up to T20 with a multivariate method -------

# Import the dataset

setwd("./final_report/data_analyses")
gaussian_data <- read.csv("gaussian_data_G17.csv")
# go back 2 levels to the original dr
Path = getwd()
setwd(dirname(dirname(Path)))

# Convert columns to factor
for (i in 1:6) {
  gaussian_data[[i]] <- as.factor(gaussian_data[[i]])
}

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

# Check gaussian data missingness
gaussian_long_NA <- gaussian_data |> 
  pivot_longer(
    cols = !(Flower_index:Subplot), 
    names_to = c(".value", "Day"), 
    names_sep = "_", 
    values_drop_na = FALSE
  )

sum(is.na(gaussian_long_NA$T))/nrow(gaussian_long_NA)*100 #3.17% missing data

# Group data by subplot and species, and count the number of flowers
strata <- gaussian_long_NA %>%
  group_by(Subplot, Type, Compound) %>%
  summarise(n = n()) 
View(strata)


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

# Use table 1 to get an idea about the balance of the data
table1(~ Width | Subplot * Type, data = gaussian_long)

table1(~ Compound | Subplot * Type, data = gaussian_data) # very balance design

# Reorder gaussian_long by Flower_index and Day
gaussian_long_arr <- gaussian_long |> 
  arrange(Flower_index, Day) |>
  mutate(Compound = factor(Compound,
             labels = c(
               "Distilled Water",
               "Apathic Acid",
               "Beerse Brew",
               "Concentrate of Caducues",
               "Distilled of Discovery",
               "Essence of Epiphanea",
               "Four in December",
               "Granule of Geheref",
               "Kar-Hamel Mooh",
               "Lucifer’s Liquid",
               "Noospherol",
               "Oil of John’s Son",
               "Powder of Perlimpinpin",
               "Spirit of Scienza",
               "Zest of Zen"
             )
           ))

# Plot the Width of the flowers over Day, color by Flower_index, facet by Compound
EDA_c <- ggplot(gaussian_long_arr, aes(x = Day, y = Width, group = Flower_index)) +
  geom_line(size = 0.3) +
  facet_wrap(~ Compound, ncol = 4, labeller = label_value) +
  theme_minimal() +
  ylab("Flower width (cm)") +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = "bold", family = "sans"), 
        strip.text.x = element_text(hjust = 0.5, size = 5, family = "sans"),
        axis.title = element_text(size = 7, family = "sans"),
        axis.text = element_text(size = 7, family = "sans"),
        legend.position = "none",
        panel.grid.minor = element_blank())

# export plots
setwd("./final_report/data_analyses/eda_plots")
ggsave("EDA_c.png", plot = EDA_c, width = 10, height = 8, units = "cm")

### Ermi's plot ---------------------------------------------------------------

### Plots ###

## Aggregated Datasets by Compound, Type, Garden, Subplot ## 
gaus_Comp <- gaussian_long %>%
  group_by(Day, Compound) %>%
  summarise(
    MeanWidth = mean(Width, na.rm = T),
    SDWidth   = sd(Width, na.rm = T),
    N         = n()
  )

gaus_Comp$Compound <-
  factor(
    gaus_Comp$Compound,
    labels = c(
      "Distilled Water",
      "Apathic Acid",
      "Beerse Brew",
      "Concentrate of Caducues",
      "Distilled of Discovery",
      "Essence of Epiphanea",
      "Four in December",
      "Granule of Geheref",
      "Kar-Hamel Mooh",
      "Lucifer’s Liquid",
      "Noospherol",
      "Oil of John’s Son",
      "Powder of Perlimpinpin",
      "Spirit of Scienza",
      "Zest of Zen"
    )
  )

gaus_Comp$MeanWidth <- round(gaus_Comp$MeanWidth, 2)
gaus_Comp$SDWidth <- round(gaus_Comp$SDWidth, 2)

gaus_Type <- gaussian_long %>%
  group_by(Day, Type) %>%
  summarise(
    MeanWidth = mean(Width, na.rm = T),
    SDWidth   = sd(Width, na.rm = T),
    N         = n()
  )

gaus_Garden <- gaussian_long %>%
  group_by(Day, Garden) %>%
  summarise(
    MeanWidth = mean(Width, na.rm = T),
    SDWidth   = sd(Width, na.rm = T),
    N         = n()
  )

gaus_Subplot <- gaussian_long %>%
  group_by(Day, Subplot) %>%
  summarise(
    MeanWidth = mean(Width, na.rm = T),
    SDWidth   = sd(Width, na.rm = T),
    N         = n()
  )

# Flower Width by Compound for T0-T20
ggplot(data = gaus_Comp) + 
  geom_line(mapping = aes(x = time, group = Compound, y = MeanWidth))
## with Colours
flowerwidth_compound <- ggplot(data = gaus_Comp, aes(x = Day, y = MeanWidth, group = Compound)) + 
  geom_line(aes(col = Compound), size = 0.3) +
  geom_point(aes(col = Compound), size = 1.9) +
  theme_minimal() +
  ylab("Mean flower width (cm)") +
  theme(axis.title = element_text(size = 9, family = "sans"),
        axis.text = element_text(size = 9, family = "sans"),
        legend.title = element_text(size = 8, family = "sans"),
        legend.text = element_text(size = 8, family = "sans"),
        panel.grid.minor = element_blank())
ggsave("flowerwidth_compound.png", plot = flowerwidth_compound, width = 12, height = 11, units = "cm")

# --> Compound 6 hast the smallest Width (considered the freshest at T20)


## Flower Width by Type for T0-T20
flowerwidth_type <- ggplot(data = gaus_Type, aes(x = Day, y = MeanWidth, group = Type)) + 
  geom_line(aes(col = Type), size = 0.3) +
  geom_point(aes(col = Type), size = 1.9) +
  theme_minimal() +
  ylab("Mean flower width (cm)") +
  theme(axis.title = element_text(size = 9, family = "sans"),
        axis.text = element_text(size = 9, family = "sans"),
        legend.title = element_text(size = 8, family = "sans"),
        legend.text = element_text(size = 8, family = "sans"),
        panel.grid.minor = element_blank())
ggsave("flowerwidth_type.png", plot = flowerwidth_type, width = 10, height = 9, units = "cm")
# --> species 2 withers faster than species 1

## Flower Width by Garden for T0-T20
flowerwidth_garden <- ggplot(data = gaus_Garden, aes(x = Day, y = MeanWidth, group = Garden)) + 
  geom_line(aes(col = Garden), size = 0.3) +
  geom_point(aes(col = Garden), size = 1.9) +
  theme_minimal() +
  ylab("Mean flower width (cm)") +
  theme(axis.title = element_text(size = 9, family = "sans"),
        axis.text = element_text(size = 9, family = "sans"),
        legend.title = element_text(size = 8, family = "sans"),
        legend.text = element_text(size = 8, family = "sans"),
        panel.grid.minor = element_blank())
ggsave("flowerwidth_garden.png", plot = flowerwidth_garden, width = 10, height = 9, units = "cm")

# No descriptive difference in Flower Width at any time point, similar pattern

## Flower Width by Subplot for T0-T20
flowerwidth_subplot <- ggplot(data = gaus_Subplot, aes(x = Day, y = MeanWidth, group = Subplot)) + 
  geom_line(aes(col = Subplot), size = 0.3) +
  geom_point(aes(col = Subplot), size = 1.9) +
  theme_minimal() +
  ylab("Mean flower width (cm)") +
  theme(axis.title = element_text(size = 12, family = "sans"),
        axis.text = element_text(size = 12, family = "sans"),
        legend.title = element_text(size = 11, family = "sans"),
        legend.text = element_text(size = 11, family = "sans"),
        panel.grid.minor = element_blank())
ggsave("flowerwidth_subplot.png", plot = flowerwidth_subplot, width = 13, height = 13, units = "cm")

# go back 3 levels to the original dr
Path = getwd()
setwd(dirname(dirname(dirname(Path))))

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

# Henry's code - question b ---------------------------------------------------

library(lme4)
library(reshape2)
library(ggplot2)
library(dplyr)
library(geepack)
library(glmmTMB)
library(MuMIn)
library(lattice)


#Question b
gaussianData <- read.csv("gaussian_data_G17.csv")
gaussianLong <- melt(gaussianData, id.vars = c("Flower_index", "Compound", "Rater", "Type", "Garden", "Subplot"), variable.name = "Day", value.name = "Width")
gaussianLong$Day <- as.numeric(gsub("T_", "", gaussianLong$Day))
gaussianLong$Compound <- as.factor(gaussianLong$Compound)
gaussianLong$Type <- as.factor(gaussianLong$Type)
gaussianLong$Subplot <- as.factor(gaussianLong$Subplot)
gaussianLong$Rater <- as.factor(gaussianLong$Rater)
gaussianLong$Flower_index <- as.factor(gaussianLong$Flower_index)
# Handling missing data using complete case analysis
gaussianLong <- na.omit(gaussianLong)
write.csv(gaussianLong, "gaussian_long.csv", row.names = FALSE)
lapply(gaussianLong[,sapply(gaussianLong, is.factor)], levels)
# Rater has only one level and is thus not considered

# m0: Additive model - analyzes the main effects of Day, Compound, Type, and Garden with random intercepts for Subplot
m0 <- lmer(Width ~ Day + Compound + Type + Garden + (1 | Subplot), data = gaussianLong)
summary(m0)
AIC(m0)
BIC(m0)
r2_m0 <- r.squaredGLMM(m0)
print(r2_m0)

# m1: Interaction model - includes the interaction between Day and Compound, plus the main effects of Type and Garden, with random intercepts for Subplot
m1 <- lmer(Width ~ Day * Compound + Type + Garden + (1 | Subplot), data = gaussianLong)
summary(m1)
AIC(m1)
BIC(m1)
r2_m1 <- r.squaredGLMM(m1)
print(r2_m1)

# m2: Additive model - Linear mixed model analyzing the effect of Day, Compound, and Type on Width, with random intercepts for Subplot
m2 <- lmer(Width ~ Day + Compound + Type + (1 | Subplot), data = gaussianLong)
summary(m2)
AIC(m2)
BIC(m2)
r2_m2 <- r.squaredGLMM(m2)
print(r2_m2)
# m3: Interaction model - includes the interaction between Day and Compound along with the main effect of Type, featuring random intercepts for Subplot
m3 <- lmer(Width ~ Day * Compound + Type + (1 | Subplot), data = gaussianLong)
summary(m3)
AIC(m3)
BIC(m3)
r2_m3 <- r.squaredGLMM(m3)
print(r2_m3)
# m4: Additive model - Adds random intercepts for Flower_index in addition to Subplot
m4 <- lmer(Width ~ Day + Compound + (1|Flower_index) + (1|Subplot), data = gaussianLong)
summary(m4)
AIC(m4)
BIC(m4)
r2_m4 <- r.squaredGLMM(m4)
print(r2_m4)
# m5: Interaction model - Includes interaction between Day and Compound and random intercepts for both Flower_index and Subplot
m5 <- lmer(Width ~ Day * Compound + (1|Flower_index) + (1|Subplot), data = gaussianLong)
summary(m5)
AIC(m5)
BIC(m5)
r2_m5 <- r.squaredGLMM(m5)
print(r2_m5)

# m6: Additive model - includes random slopes for Day:Compound within Subplot 
m6 <- lmer(Width ~ Day + Compound + Type + (1 + Compound:Day || Subplot), data = gaussianLong)
summary(m6)
AIC(m6)
BIC(m6)
r2_m6 <- r.squaredGLMM(m6)
print(r2_m6)
# m7: Interaction model - features the interaction between Day and Compound, main effect of Type, and random slopes for the Day:Compound interaction within Subplot
m7 <- lmer(Width ~ Day * Compound + Type + (1 + Compound:Day || Subplot), data = gaussianLong)
summary(m7)
AIC(m7)
BIC(m7)
r2_m7 <- r.squaredGLMM(m7)
print(r2_m7)

# m8: Interaction model - includes Day and Compound interaction with random slopes for this interaction on Flower_index and random intercepts for Subplot
m8 <- lmer(Width ~ Day * Compound + (1 + Compound:Day | Flower_index) + (1 | Subplot), data = gaussianLong)
summary(m8)
AIC(m8)
BIC(m8)
r2_m8 <- r.squaredGLMM(m8)
print(r2_m8)

# m9: Interaction model - examines main and interaction effects of Day and Compound, with random slopes for these interactions on Flower_index, and random intercepts for Subplot
m9 <- lmer(Width ~ Day * Compound + (1 + Compound:Day | Flower_index) + (1 | Subplot), data = gaussianLong)
summary(m9)
AIC(m9)
BIC(m9)
r2_m9 <- r.squaredGLMM(m9)
print(r2_m9)


# Residuals vs fitted values for m5
residuals_data <- data.frame(
  Fitted = fitted(m5),
  Residuals = resid(m5)
)

residualPlot <- ggplot(residuals_data, aes(x = Fitted, y = Residuals)) +
  geom_point() +
  geom_smooth(method = "loess", color = "red") +
  labs(title = "Residuals vs Fitted Values", x = "Fitted Values", y = "Residuals")

print(residualPlot)


# Cluster-specific effects check
randomEffectsPlot <- dotplot(ranef(m5, condVar=TRUE))
print(randomEffectsPlot)

### re-estimate models for comparison:
m0 <- lmer(Width ~ Day + Compound (1 | Flower_index) + (1 | Subplot), data = gaussianLong)
m1 <- lmer(Width ~ Day*Compound + (1 | Flower_index) + (1 | Subplot), data = gaussianLong)
m2 <- lmer(Width ~ Day:Compound + Day + (1 | Flower_index) + (1 | Subplot), data = gaussianLong)
m3 <- lmer(Width ~ Day:Compound + Day + Type + (1 | Flower_index) + (1 | Subplot), data = gaussianLong)
m4 <- lmer(Width ~ Day:Compound + Day + Type + Garden + (1 | Flower_index) + (1 | Subplot), data = gaussianLong)

m5 <- lmer(Width ~ Day * Compound + Type + (1|Flower_index) + (1|Subplot), data = gaussianLong)

anova(m0, m1)
anova(m1, m2) # not significantly better, but subplot is left in model anyway
anova(m2, m3)
anova(m3, m4)

aicm0 <- AIC(m0)
bicm0 <- BIC(m0)
r2m0 <- r.squaredGLMM(m0)

aicm1 <- AIC(m1)
bicm1 <- BIC(m1)
r2m1 <- r.squaredGLMM(m1)

aicm2 <- AIC(m2)
bicm2 <- BIC(m2)
r2m2 <- r.squaredGLMM(m2)

aicm3 <- AIC(m3)
bicm3 <- BIC(m3)
r2m3 <- r.squaredGLMM(m3)

aicm4 <- AIC(m4)
bicm4 <- BIC(m4)
r2m4 <- r.squaredGLMM(m4)

lmer_dat <- data.frame(
  model_name = c(0, 1, 2, 3, 4),
  description = c(
    "Compound + Day + FlowerID  + SubplotID",
    "Compound*Day + FlowerID + SubplotID",
    "Compound:Day + Day + FlowerID + SubplotID",
    "Compound:Day + Day + Species + Type + FlowerID + SubplotID",
    "Compound:Day + Day + Species + Garden + FlowerID + SubplotID"
  ),
  aic = c(aicm0, aicm1, aicm2, aicm3, aicm4),
  bic = c(bicm0, bicm1, bicm2, bicm3, bicm4),
  r2m = c(r2m0[1], r2m1[1], r2m2[1], r2m3[1], r2m4[1])
)



lmer_tab_sel <- lmer_dat %>%
  gt() %>%
  tab_header(
    title = "Model Selection Process",
    subtitle = "Comparison of Models based on AIC, BIC, R"
  ) %>%
  cols_label(
    model_name = "Model",
    aic = "AIC",
    bic = "BIC",
    r2m = "R^2",
    description = "Formula"
  ) %>%
  cols_align(
    align = "center"
  )


#setwd("./tables")
gtsave(lmer_tab_sel, "tab_lmer_sel.html")
gtsave(lmer_tab_sel, "tab_lmer_sel.tex")
webshot("tab_lmer_sel.html", "tab_lmer_sel.pdf")

##########################################################

# same analysis with renamed compounds
gaussianLongArr <- gaussianLong |> 
  arrange(Flower_index, Day) |>
  mutate(Compound = factor(Compound,
                           labels = c(
                             "Distilled Water",
                             "Apathic Acid",
                             "Beerse Brew",
                             "Concentrate of Caducues",
                             "Distilled of Discovery",
                             "Essence of Epiphanea",
                             "Four in December",
                             "Granule of Geheref",
                             "Kar-Hamel Mooh",
                             "Lucifer’s Liquid",
                             "Noospherol",
                             "Oil of John’s Son",
                             "Powder of Perlimpinpin",
                             "Spirit of Scienza",
                             "Zest of Zen"
                           )
  ))
# 
m3 <- lmer(Width ~ Day:Compound + Day + Type
           + (1 | Flower_index) + (1 | Subplot), data = gaussianLongArr)
lmer_gaussian <- m3


tidy_lmer <- broom.mixed::tidy(lmer_gaussian) %>%
  rename(
    Term = term,
    Estimate = estimate,
    `Standard Error` = std.error,
    `Degrees of Freedom` = df,
    `Statistic` = statistic,
    `p-value` = p.value
  ) %>%
  filter(!is.na(`Degrees of Freedom`)) %>%  # Exclude random effects
  dplyr::select(-group, -effect) %>%  # Drop unnecessary columns
  mutate(
    #Term = sub("^[^:]*:", "", Term),  # Remove variable name and keep only the label
    `p-value` = ifelse(`p-value` < .001, "<.001", as.character(`p-value`)),
    `Degrees of Freedom` = round(`Degrees of Freedom`, 0)
  )
# for table: remove "Compound" from each compound name
#tidy_lmer$Term[4:17] <- sub("Compound", "", tidy_lmer$Term[4:17]) %>% trimws()
tidy_lmer$Term[4:17] <- c("Apathic Acid : day", "Beerse Brew : day", "Concentrate of Caducues : day", 
                        "Distilled of Discovery", "Essence of Epiphanea : day", "Four in December : day", 
                        "Granule of Geheref : day", "Kar-Hamel Mooh : day", "Lucifer’s Liquid : day", 
                        "Noospherol : day", "Oil of John’s Son : day", "Powder of Perlimpinpin : day", 
                        "Spirit of Scienza : day", "Zest of Zen : day")


# Create and format output table
lmer_tab <- tidy_lmer %>%
  gt() %>%
  tab_header(
    title = "Linear Mixed-Effects Model Results",
    subtitle = "Summary of Fixed Effects"
  ) %>%
  fmt_number(
    columns = vars(Estimate, `Standard Error`, Statistic),
    decimals = 2
  ) %>%
  #fmt_scientific(
  #  columns = vars(`p-value`),
  #  decimals = 3
  #) %>%
  cols_label(
    Term = "Term",
    Estimate = "Estimate",
    `Standard Error` = "Std. Error",
    `Degrees of Freedom` = "DF",
    Statistic = "t-value",
    `p-value` = "p-value"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_options(
    table.font.size = "small",
    table.align = "left"
  )

setwd("./tables")
gtsave(lmer_tab, "tab_lmer.html")
gtsave(lmer_tab, "tab_lmer.tex")
webshot("tab_lmer.html", "tab_lmer.pdf")
######

lmer_gaussian <- m3
summary(lmer_gaussian)

# compounds 2, 3, 5, 6, 9,11, 13, 14, 15 show smaller width over time compared to water
# --> contrast for the 4 largest coefficients


negative_compounds <- c("Compound2", "Compound3", "Compound5", "Compound6")

margmeans <- emmeans(lmer_gaussian, ~ Compound * Day)
# margmeans <- emmeans(lmer_gaussian, ~ Compound | Day, pbkrtest.limit = 3660)

contrasts <- emmeans::contrast(margmeans, method = "pairwise", simple = "each", combine = TRUE, 
                               adjust = "bonferroni")

negative_compounds <- c("Apathic Acid - Beerse Brew", 
                        "Apathic Acid - Distilled of Discovery", 
                        "Apathic Acid - Essence of Epiphanea", 
                        "Beerse Brew - Distilled of Discovery", 
                        "Beerse Brew - Essence of Epiphanea", 
                        "Distilled of Discovery - Essence of Epiphanea")

filtered_contrasts <- contrasts %>%
  as.data.frame() %>%
  mutate(p.value = ifelse(p.value < .001, "<.001", p.value)) %>%
  filter(contrast %in% negative_compounds)

filtered_contrasts


#filtered_contrasts$p.value <- ifelse(filtered_contrasts$p.value < .001, "<.001", filtered_contrasts$p.value)

# filtered_contrasts <- filtered_contrasts %>%
#   mutate(p.value = paste0("$", p.value, "$"))


# Contrast table
contrast_results <- data.frame(
  contrast = filtered_contrasts$contrast,
  estimate = filtered_contrasts$estimate,
  SE = filtered_contrasts$SE,
  df = filtered_contrasts$df,
  t.ratio = filtered_contrasts$t.ratio,
  p.value = filtered_contrasts$p.value
)

tab_lmer_contr <- contrast_results %>%
  gt() %>%
  tab_header(
    title = "Contrast Results",
    subtitle = "Estimated Marginal Means Contrasts"
  ) %>%
  cols_label(
    contrast = "Contrast",
    estimate = "Estimate",
    SE = "SE",
    df = "DF",
    t.ratio = "t-value",
    p.value = "p-value"
  ) %>%
  fmt_number(
    columns = vars(estimate, SE, t.ratio),
    decimals = 2,
    #drop_trailing_zeros = TRUE
  ) %>%
  fmt_number(
    columns = vars(p.value),
    decimals = 3
  ) %>%
  cols_width(contrast ~ px(180))

# save
setwd("./tables")
gtsave(tab_lmer_contr, "tab_lmer_contr.html")
gtsave(tab_lmer_contr, "tab_lmer_contr.tex")
webshot("tab_lmer_contr.html", "tab_lmer_contr.pdf")


####### Question d ######

# Transforming data for PCA
pcaData <- dcast(gaussianLong, Flower_index + Compound ~ Day, value.var = "Width")

if (sum(is.na(pcaData)) > 0) {
  pcaData <- pcaData %>%
    mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .)))
}

pcaData[, -c(1, 2)] <- lapply(pcaData[, -c(1, 2)], as.numeric)

pcaResults <- prcomp(pcaData[, -c(1, 2)], center = TRUE, scale. = TRUE)

# Check results
summary(pcaResults)

# Screeplot and Biplot for PCA
# We notice that most of the variance is explained by the first Principal Component
screeplot(pcaResults, type = "lines")
# PC1 captures a significant variance component, given the clustering of points along PC1 axis
biplot(pcaResults)

# PCA scores
scores <- data.frame(pcaResults$x)
names(scores) <- paste("PC", 1:ncol(scores), sep = "")
scores <- cbind(pcaData[, 1:2], scores)  # Flower_index and Compound as first two columns in pcaData

# Visualizing patterns
# Further detailed analysis is needed to understand the role of 'Compound' in the dataset, which cannot be explained by PCA Scores only
ggplot(scores, aes(x = PC1, y = PC2, color = Compound)) +
  geom_point() +
  labs(title = "PCA Scores by Compound", x = "PC1", y = "PC2") +
  theme_minimal()

pcaData$meanWidth <- aggregate(Width ~ Flower_index + Compound, data = gaussianLong, FUN = mean)$Width

combined_data <- merge(scores, pcaData, by = c("Flower_index", "Compound"))

# Regression analysis using PCA components
lm_pca <- lm(meanWidth ~ PC1 + PC2, data = combined_data)
summary(lm_pca)

# Question d: Analyze T0 up to T20 with a multivariate method -------

## Multivariate regression analysis
# follow the post at 
# https://library.virginia.edu/data/articles/getting-started-with-multivariate-multiple-regression

# Regress T0 to T20 on Compound, Type, and Subplot (the saturated model)
# , with gaussian_data
mlm1 <- lm(cbind(T_0, T_1, T_2, T_3, T_4, T_5, T_6, T_7, T_8, T_9, T_10, T_11, T_12, T_13, T_14, T_15, T_16, T_17, T_18, T_19, T_20) ~ Compound + Type + Subplot, data = gaussian_data)
summary(mlm1)

# variance - covariance matrix
vcov(mlm1) # too complicated

# check what predictor is multivariately significant
# Either Anova() or Manova() can be used
Anova <- Anova(mlm1) # all 3 factors are significant

# check if we can remove Subplot
mlm2 <- update(mlm1, . ~ . - Subplot)
anova(mlm1, mlm2) #different from Anova() # Can't remove subplot

# mlm1 is the best model. 

# we are interested in Compound with negative coefficients in mlm1
# extract these compounds


# Questions e and f ---------------------------------------------------------

# model 01
fixed_slope_lmer <- lmer(
  Width ~ Day +
    Day:Compound + 
    Garden +
    Type + 
    (1 | Subplot) +
    (1 | Flower_index),
  data = gaussian_long 
)
# why there's no compound here?

anova(fixed_slope_lmer)

print(summary(fixed_slope_lmer), correlation = FALSE)

AIC(fixed_slope_lmer) #4717.925
BIC(fixed_slope_lmer) #4848.235

# model 02
random_slope_lmer <- lmer(
  Width ~ Day +
    Compound:Day + 
    Type + 
    Garden +
    (1 | Subplot) +
    (Day | Flower_index),
  data = gaussian_long,
  control = lmerControl(optimizer = "optimx",
                      optCtrl = list(method = "nlminb")) 
)

anova(random_slope_lmer)

print(summary(random_slope_lmer), correlation = FALSE)

AIC(random_slope_lmer) #4175.009
BIC(random_slope_lmer) #4317.729
