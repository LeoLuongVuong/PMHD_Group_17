library(lme4)
library(reshape2)
library(ggplot2)
library(dplyr)
library(geepack)
library(glmmTMB)
library(MuMIn)
library(lattice)


#Question b
gaussianData <- read.csv("C:/Users/Moi/Downloads/gaussian_data_G17.csv")
gaussianLong <- melt(gaussianData, id.vars = c("Flower_index", "Compound", "Rater", "Type", "Garden", "Subplot"), variable.name = "Day", value.name = "Width")
gaussianLong$Day <- as.numeric(gsub("T_", "", gaussianLong$Day))
gaussianLong$Compound <- as.factor(gaussianLong$Compound)
gaussianLong$Type <- as.factor(gaussianLong$Type)
gaussianLong$Subplot <- as.factor(gaussianLong$Subplot)
gaussianLong$Rater <- as.factor(gaussianLong$Rater)

for(col in names(gaussianLong)) {
  if(is.numeric(gaussianLong[[col]])) {
    gaussianLong[[col]][is.na(gaussianLong[[col]])] <- median(gaussianLong[[col]], na.rm = TRUE)
  }
}
# Simple Random Intercept Model 
m1 <- lmer(Width ~ Day + Compound + Type + (1 | Subplot), data = gaussianLong, na.action = na.fail)
summary(m1)
AIC(m1)  
BIC(m1)  

# Including Random Slopes
m2 <- lmer(Width ~ Day + Compound + Type + (1 + Day || Subplot), data = gaussianLong, na.action = na.fail)
summary(m2)
AIC(m2)  
BIC(m2)  

# Crossed Random Effects 
m3 <- lmer(Width ~ Day + Compound + Type + (1 | Subplot), data = gaussianLong, na.action = na.fail)
summary(m3)
AIC(m3)  
BIC(m3)  

# Model with Interaction Effects
m4 <- lmer(Width ~ Day * Compound + Type + (1 | Subplot), data = gaussianLong, na.action = na.fail)
summary(m4)
AIC(m4)  
BIC(m4)  

# AIC and BIC values
aicBicValues <- data.frame(
  Model = c("Model 1", "Model 2", "Model 3", "Model 4"),
  AIC = c(AIC(m1), AIC(m2), AIC(m3), AIC(m4)),
  BIC = c(BIC(m1), BIC(m2), BIC(m3), BIC(m4))
)
print(aicBicValues)

# Stepwise model selection for GLMM, showing m4 as the best model
candidateModels <- dredge(m4)

# Checking the best models based on AIC
topModels <- get.models(candidateModels, subset = delta < 2)  # Selecting models with ΔAIC < 2


lapply(topModels, summary)

# AIC scores of all candidate models
print(candidateModels)

# Generalized Estimating Equations
geeModel <- geeglm(Width ~ Day + Compound + Type, id=Subplot, data=gaussianLong, family=gaussian, corstr="exchangeable")
summary(geeModel)

if (length(topModels) > 0) {
  stepwiseModel <- topModels[[1]]  
} else {
  stepwiseModel <- m4  
}
# Marginal R-squared from glmmTMB
marginalR2 <- r.squaredGLMM(stepwiseModel)
print(marginalR2)

# Residuals vs fitted values for m4
residuals_data <- data.frame(
  Fitted = fitted(stepwiseModel),
  Residuals = resid(stepwiseModel)
)

residualPlot <- ggplot(residuals_data, aes(x = Fitted, y = Residuals)) +
  geom_point() +
  geom_smooth(method = "loess", color = "red") +
  labs(title = "Residuals vs Fitted Values", x = "Fitted Values", y = "Residuals")

print(residualPlot)


# Cluster-specific effects check
randomEffectsPlot <- dotplot(ranef(m4, condVar=TRUE))
print(randomEffectsPlot)

#Question d

# Transforming data for PCA
pcaData <- dcast(gaussianLong, Flower_index + Compound ~ Day, value.var = "Width")
pcaResults <- prcomp(pcaData[, -c(1, 2)], center = TRUE, scale. = TRUE)

# Screeplot and Biplot for PCA
screeplot(pcaResults, type = "lines")
biplot(pcaResults)

# PCA scores
scores <- data.frame(pcaResults$x)
names(scores) <- paste("PC", 1:ncol(scores), sep = "")
scores <- cbind(pcaData[, 1:2], scores)  # Flower_index and Compound as first two columns in pcaData

# Visualizing patterns
ggplot(scores, aes(x = PC1, y = PC2, color = Compound)) +
  geom_point() +
  labs(title = "PCA Scores by Compound", x = "PC1", y = "PC2") +
  theme_minimal()

pcaData$meanWidth <- aggregate(Width ~ Flower_index + Compound, data = gaussianLong, FUN = mean)$Width

combined_data <- merge(scores, pcaData, by = c("Flower_index", "Compound"))

# Regression analysis using PCA components
lm_pca <- lm(meanWidth ~ PC1 + PC2, data = combined_data)
summary(lm_pca)

# Residuals vs fitted values to check for patterns of non-linearity or heteroscedasticity
residuals_data <- data.frame(Fitted = fitted(lm_pca), Residuals = resid(lm_pca))
ggplot(residuals_data, aes(x = Fitted, y = Residuals)) +
  geom_point() +
  geom_smooth(method = "loess", color = "red") +
  labs(title = "Residuals vs Fitted Values", x = "Fitted Values", y = "Residuals")