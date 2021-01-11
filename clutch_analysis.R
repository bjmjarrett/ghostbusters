library(tidyverse)
#library(rethinking) # not required for this analysis
#library(gridExtra) # not required for this analysis
library(brms)
#library(gameofthrones) # not required for this analysis
#library(emmeans) # not required for this analysis
#library(modelr) # not required for this analysis
#library(tidybayes) # not required for this analysis

#### ANALYSIS for the clutch size data (number of eggs visible on the bottom of the breeding box) for the Ghosbusters experiment.

# import data
data <- read.csv("epistasis_raw_data.csv",header=T)

# convert variables
data$MALE_SIZE <- as.numeric(levels(data$MALE_SIZE))[data$MALE_SIZE]
data$FEMALE_SIZE <- as.numeric(levels(data$FEMALE_SIZE))[data$FEMALE_SIZE]
data$MALE_AGE <- as.numeric(levels(data$MALE_AGE))[data$MALE_AGE]
data$FEMALE_AGE <- as.numeric(levels(data$FEMALE_AGE))[data$FEMALE_AGE]
data$hole <- as.factor(data$hole)
data$no_larvae <- as.numeric(levels(data$no_larvae))[data$no_larvae]
data$dev_hours <- as.numeric(levels(data$dev_hours))[data$dev_hours]

# check data dimensions
str(data)
table(data$treat_cross, data$treatment)
table(data$cross, data$treatment)

# clean up data by removing boxes without eggs
data <- data[!is.na(data$no_larvae),]
data <- filter(data, no_eggs > 0)

# create hybrid variable, where FN and NF are hybrids
data <- data %>% 
  mutate(hybrid_fac = as.factor(ifelse(treat_cross == "FN", "hybrid",
                                       ifelse(treat_cross == "NF", "hybrid",
                                              ifelse(treat_cross == "FF", "full_full", "no_no")))))

# dummy variables for male and female backgrounds
data$male_F <- ifelse(data$male_treatment == "F", 1, 0)
data$female_F <- ifelse(data$female_treatment == "F", 1, 0)

#### MODELS

#remove NAs from no_larvae
d <- data %>% 
  drop_na(no_larvae, FEMALE_SIZE, MALE_SIZE, carcass)

# check data dimensions
table(d$treat_cross, d$treatment)
table(d$cross, d$treatment)


# standardise variables
d$carcass.s <- (d$carcass - mean(d$carcass)) / sd(d$carcass)
d$female_size.s <- (d$FEMALE_SIZE - mean(d$FEMALE_SIZE)) / sd(d$FEMALE_SIZE)
d$male_size.s <- (d$MALE_SIZE - mean(d$MALE_SIZE)) / sd(d$MALE_SIZE)


### 1. clutch size hybrid analysis

# dataframe for clutch size - all broods are included because the care treatment has not been imposed
d_eggs <- d %>%
  select(male_F, female_F, carcass.s, female_size.s, male_size.s, no_eggs, hybrid_fac)

# initial data plotting
ggplot(d_eggs, aes(y = no_eggs, x = hybrid_fac)) +
  geom_boxplot()

# first model with all variables
m_hyb_eggs_1 <- brm(no_eggs ~ hybrid_fac + carcass.s + female_size.s + male_size.s,
                    data = d_eggs, 
                    family = "negbinomial",
                    prior = set_prior('normal(0,5)'),
                    chains = 4, cores = 4,
                    control = list(adapt_delta = 0.99, max_treedepth = 15),
                    save_all_pars = TRUE)

summary(m_hyb_eggs_1)
pp_check(m_hyb_eggs_1)

# hypothesis testing for contrasts between hybrid factor levels
hyp <- hypothesis(m_hyb_eggs_1, c("Intercept = Intercept + hybrid_fachybrid", "Intercept = Intercept + hybrid_facno_no", "Intercept + hybrid_facno_no = Intercept + hybrid_fachybrid"), alpha = 0.05)

plot(hyp)

# updated model removing the hybrid factor
m_hyb_eggs_2 <- update(m_hyb_eggs_1,
                       formula. = no_eggs ~ carcass.s + female_size.s + male_size.s)
summary(m_hyb_eggs_2)

# adding loo and waic criteria
m_hyb_eggs_1 <- add_criterion(m_hyb_eggs_1, c("waic", "loo"))
m_hyb_eggs_2 <- add_criterion(m_hyb_eggs_2, c("waic", "loo"))

# compares models
w <- loo_compare(m_hyb_eggs_1, m_hyb_eggs_2, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

# model weights about how likely the model is to best predict future data
model_weights(m_hyb_eggs_1, m_hyb_eggs_2, weights = "loo")

### SECOND analysis focusing on the background of the males and females

### CLUTCH 2

# check data
ggplot(aes(y = no_eggs, x = male_F, col = as.factor(female_F)),
       data = d_eggs) +
  geom_point(position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2),
             alpha = 0.5) 

# first model with interaction
m_eggs_1 <- brm(no_eggs ~ female_F * male_F + carcass.s + female_size.s + male_size.s,
                data = d_eggs, 
                family = "negbinomial",
                prior = set_prior('normal(0,5)'),
                chains = 4, cores = 4,
                control = list(adapt_delta = 0.99, max_treedepth = 15),
                save_all_pars = TRUE)

summary(m_eggs_1)
pp_check(m_eggs_1)

m_eggs_2 <- update(m_eggs_1,
                   formula. = no_eggs ~ female_F + male_F + 
                     carcass.s + female_size.s + male_size.s)

m_eggs_3 <- update(m_eggs_1,
                   formula. = no_eggs ~ female_F + 
                     carcass.s + female_size.s + male_size.s)

m_eggs_4 <- update(m_eggs_1,
                   formula. = no_eggs ~ male_F + 
                     carcass.s + female_size.s + male_size.s)

m_eggs_5 <- update(m_eggs_1,
                   formula. = no_eggs ~ carcass.s + female_size.s + male_size.s)


m_eggs_1 <- add_criterion(m_eggs_1, c("waic", "loo"))
m_eggs_2 <- add_criterion(m_eggs_2, c("waic", "loo"))
m_eggs_3 <- add_criterion(m_eggs_3, c("waic", "loo"))
m_eggs_4 <- add_criterion(m_eggs_4, c("waic", "loo"))
m_eggs_5 <- add_criterion(m_eggs_5, c("waic", "loo"))

w <- loo_compare(m_eggs_1, m_eggs_2, m_eggs_3, m_eggs_4, m_eggs_5, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_eggs_1, m_eggs_2, m_eggs_3, m_eggs_4, m_eggs_5, weights = "loo")
