library(gee)
library(tidyr)
library(tidyverse)
library(lme4)
library(ggplot2)

dat_gaus <- read.csv("final_report/data_analyses/gaussian_data_G17.csv", sep = ",", header = T)

head(dat_gaus)
summary(dat_gaus)

dat_gaus$Flower_index <- as.factor(dat_gaus$Flower_index)
dat_gaus$Compound <- as.factor(dat_gaus$Compound)
dat_gaus$Rater <- as.factor(dat_gaus$Rater)
dat_gaus$Type <- as.factor(dat_gaus$Type)
dat_gaus$Garden <- as.factor(dat_gaus$Garden)
dat_gaus$Subplot <- as.factor(dat_gaus$Subplot)

gaus_long <- dat_gaus %>% pivot_longer(cols = starts_with("T_"), names_to = "time", values_to = "Width")

gaus_long$time <- factor(gaus_long$time, levels = c(
  "T_0", "T_1", "T_2", "T_3", "T_4", "T_5", "T_6", "T_7", "T_8", "T_9", "T_10", 
  "T_11","T_12","T_13","T_14","T_15","T_16","T_17","T_18","T_19","T_20"))


gaus_long$time <- as.factor(gaus_long$time)

## check balanced design
table(gaus_long$time)  # N = 180, complete observations for all 20 time points

gaus_long %>% count(Flower_index)  # n = 21, equal number of time points per flower 

ggplot(data = gaus_long) +
  geom_bar(mapping = aes(x = Compound))  # n = 250, equal number of flowers per compoung



## check out distribution of flower Width
ggplot(data = gaus_long) + 
  geom_histogram(mapping = aes(x = Width), binwidth = 0.5)


## Missings
missings_gauss <- gaus_long %>%
  group_by(time) %>%
  summarise(MissingValues = sum(is.na(Width)))
# --> N = 10 missings for T13-T20

# check distribution grouped by flower (summarized over time)
gaus_agg <- gaus_long %>%
  group_by(Flower_index) %>%
  summarise(
    MeanWidth = mean(Width, na.rm = T),
    SDWidth   = sd(Width, na.rm = T),
    N         = n()
  ) %>%
  round(2)


ggplot(data = gaus_agg) + 
  geom_histogram(mapping = aes(x = MeanWidth), binwidth = 0.5)


## Aggregated Datasets by Compound, Type, Garden, Subplot ## 
gaus_Comp <- gaus_long %>%
  group_by(time, Compound) %>%
  summarise(
    MeanWidth = mean(Width, na.rm = T),
    SDWidth   = sd(Width, na.rm = T),
    N         = n()
    )
gaus_Comp$MeanWidth <- round(gaus_Comp$MeanWidth, 2)
gaus_Comp$SDWidth <- round(gaus_Comp$SDWidth, 2)

gaus_Type <- gaus_long %>%
  group_by(time, Type) %>%
  summarise(
    MeanWidth = mean(Width, na.rm = T),
    SDWidth   = sd(Width, na.rm = T),
    N         = n()
  )

gaus_Garden <- gaus_long %>%
  group_by(time, Garden) %>%
  summarise(
    MeanWidth = mean(Width, na.rm = T),
    SDWidth   = sd(Width, na.rm = T),
    N         = n()
  )

gaus_Subplot <- gaus_long %>%
  group_by(time, Subplot) %>%
  summarise(
    MeanWidth = mean(Width, na.rm = T),
    SDWidth   = sd(Width, na.rm = T),
    N         = n()
  )


### Plots ###

# Flower Width by Compound for T0-T20
ggplot(data = gaus_Comp) + 
  geom_line(mapping = aes(x = time, group = Compound, y = MeanWidth))
## with Colours
ggplot(data = gaus_Comp, aes(x = time, y = MeanWidth, group = Compound)) + 
  geom_line(aes(col = Compound), size = 0.3) +
  geom_point(aes(col = Compound), size = 1.9)
# --> Compound 6 hast the smallest Width (considered the freshest at T20)

## Flower Width by Type for T0-T20
ggplot(data = gaus_Type, aes(x = time, y = MeanWidth, group = Type)) + 
  geom_line(aes(col = Type), size = 0.3) +
  geom_point(aes(col = Type), size = 1.9)
# --> species 2 withers faster than species 1

## Flower Width by Garden for T0-T20
ggplot(data = gaus_Garden, aes(x = time, y = MeanWidth, group = Garden)) + 
  geom_line(aes(col = Garden), size = 0.3) +
  geom_point(aes(col = Garden), size = 1.9)
# No descriptive difference in Flower Width at any time point, similar pattern


## Flower Width by Subplot for T0-T20
ggplot(data = gaus_Subplot, aes(x = time, y = MeanWidth, group = Subplot)) + 
  geom_line(aes(col = Subplot), size = 0.3) +
  geom_point(aes(col = Subplot), size = 1.9)


# width_gee <- gee(Width ~ Compound*time, data = dat_long, id = Flower_index, corstr = "independence")

################################################################################
# Count Data


dat <- read.csv("count_data_G17.csv", sep = ",", header = T)

dat <- dat[!is.na(dat$tot.vase.days),]
dummy_matrix <- matrix(0, nrow = nrow(dat), ncol = 30)

for (i in 1:nrow(dat)) {
  dummy_matrix[i, 1:(dat$tot.vase.days[i])] <- 1
}

dat_dummy <- cbind(dat, dummy_matrix)

#names <- c(colnames(dat), paste0("T", 0:29))
names <- c(colnames(dat), 0:29)
colnames(dat_dummy) <- names
# dat_dummy <- dat_dummy %>% rename(Xtot.vase.days = tot.vase.days) # begins with "t"

dat_dummy_long <- dat_dummy %>% pivot_longer(cols = 8:37, names_to = "Day", values_to = "Freshness")
#dat_dummy_long <- dat_dummy %>% pivot_longer(cols = starts_with("T"), names_to = "Day", values_to = "Freshness")


fresh_gee <- gee(Freshness ~ Day * compound, 
                 family = "binomial", 
                 data = dat_dummy_long, 
                 id = flowerID, 
                 corstr = "exchangeable")

fresh_glmm <- glmer(Freshness ~ Day * compound + (1 | flowerID), dat_dummy_long,
                   family = binomial)
