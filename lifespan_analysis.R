library(tidyverse)
library(rethinking)
library(gridExtra)
library(brms)
library(gameofthrones)
library(emmeans)
library(modelr)
library(tidybayes)


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

data <- data %>% 
  drop_na(no_larvae, FEMALE_SIZE, MALE_SIZE, carcass)

# create hybrid variable, where FN and NF are hybrids
data <- data %>% 
  mutate(hybrid_fac = as.factor(ifelse(treat_cross == "FN", "hybrid",
                                       ifelse(treat_cross == "NF", "hybrid",
                                              ifelse(treat_cross == "FF", "full_full", "no_no")))))

# dummy variables for male and female backgrounds
data$male_F <- ifelse(data$male_treatment == "F", 1, 0)
data$female_F <- ifelse(data$female_treatment == "F", 1, 0)

# dummy variable for successful broods
data$success <- ifelse(data$no_larvae > 0, 1, 0)

data$carcass.s <- (data$carcass - mean(data$carcass)) / sd(data$carcass)
data$female_size.s <- (data$FEMALE_SIZE - mean(data$FEMALE_SIZE)) / sd(data$FEMALE_SIZE)
data$male_size.s <- (data$MALE_SIZE - mean(data$MALE_SIZE)) / sd(data$MALE_SIZE)


#####

### FULL CARE

d_age_FC <- data %>%
  select(male_F, female_F, carcass.s, female_size.s, male_size.s, 
         no_larvae, hybrid_fac, treatment, success, brood_mass,
         FEMALE_AGE, MALE_AGE) %>% 
  filter(treatment == "F",
         success == 1,
         MALE_AGE > 0,
         FEMALE_AGE > 0) 

#### MALE AGE

ggplot(aes(y = MALE_AGE, x = male_F, col = as.factor(female_F)),
       data = d_age_FC) +
  geom_point(position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2),
             alpha = 0.5) 

m_mage_FC_1 <- brm(MALE_AGE ~ female_F * male_F + 
                     carcass.s + female_size.s + male_size.s +
                     brood_mass,
                   data = d_age_FC, 
                   family = "weibull",
                   prior = set_prior('normal(0,5)'),
                   chains = 4, cores = 4,
                   control = list(adapt_delta = 0.99, max_treedepth = 15),
                   save_all_pars = TRUE)

summary(m_mage_FC_1)
pp_check(m_mage_FC_1)

m_mage_FC_2 <- update(m_mage_FC_1,
                      formula. = MALE_AGE ~ male_F + female_F + 
                        carcass.s + female_size.s + male_size.s + brood_mass)

m_mage_FC_3 <- update(m_mage_FC_1,
                      formula. = MALE_AGE ~ female_F + 
                        carcass.s + female_size.s + male_size.s + brood_mass)

m_mage_FC_4 <- update(m_mage_FC_1,
                      formula. = MALE_AGE ~ male_F + 
                        carcass.s + female_size.s + male_size.s + brood_mass)

m_mage_FC_5 <- update(m_mage_FC_1,
                      formula. = MALE_AGE ~ carcass.s + female_size.s + male_size.s + 
                        brood_mass)

m_mage_FC_1 <- add_criterion(m_mage_FC_1, c("waic", "loo"))
m_mage_FC_2 <- add_criterion(m_mage_FC_2, c("waic", "loo"))
m_mage_FC_3 <- add_criterion(m_mage_FC_3, c("waic", "loo"))
m_mage_FC_4 <- add_criterion(m_mage_FC_4, c("waic", "loo"))
m_mage_FC_5 <- add_criterion(m_mage_FC_5, c("waic", "loo"))

w <- loo_compare(m_mage_FC_1, m_mage_FC_2, m_mage_FC_3, 
                 m_mage_FC_4, m_mage_FC_5, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_mage_FC_1, m_mage_FC_2, m_mage_FC_3, 
              m_mage_FC_4, m_mage_FC_5, weights = "loo")

conditional_effects(m_mage_FC_2)
summary(m_mage_FC_2)

### FEMALE LIFESPAN

ggplot(aes(y = FEMALE_AGE, x = male_F, col = as.factor(female_F)),
       data = d_age_FC) +
  geom_point(position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2),
             alpha = 0.5) 

m_fage_FC_1 <- brm(FEMALE_AGE ~ female_F * male_F + 
                     carcass.s + female_size.s + male_size.s +
                     brood_mass,
                   data = d_age_FC, 
                   family = "weibull",
                   prior = set_prior('normal(0,5)'),
                   chains = 4, cores = 4,
                   control = list(adapt_delta = 0.99, max_treedepth = 15),
                   save_all_pars = TRUE)

summary(m_fage_FC_1)
pp_check(m_fage_FC_1)

m_fage_FC_2 <- update(m_fage_FC_1,
                      formula. = FEMALE_AGE ~ male_F + female_F + 
                        carcass.s + female_size.s + male_size.s + brood_mass)

m_fage_FC_3 <- update(m_fage_FC_1,
                      formula. = FEMALE_AGE ~ female_F + 
                        carcass.s + female_size.s + male_size.s + brood_mass)

m_fage_FC_4 <- update(m_fage_FC_1,
                      formula. = FEMALE_AGE ~ male_F + 
                        carcass.s + female_size.s + male_size.s + brood_mass)

m_fage_FC_5 <- update(m_fage_FC_1,
                      formula. = FEMALE_AGE ~ carcass.s + female_size.s + male_size.s + 
                        brood_mass)

m_fage_FC_1 <- add_criterion(m_fage_FC_1, c("waic", "loo"), reloo = T)
m_fage_FC_2 <- add_criterion(m_fage_FC_2, c("waic", "loo"), reloo = T)
m_fage_FC_3 <- add_criterion(m_fage_FC_3, c("waic", "loo"), reloo = T)
m_fage_FC_4 <- add_criterion(m_fage_FC_4, c("waic", "loo"), reloo = T)
m_fage_FC_5 <- add_criterion(m_fage_FC_5, c("waic", "loo"), reloo = T)

w <- loo_compare(m_fage_FC_1, m_fage_FC_2, m_fage_FC_3, 
                 m_fage_FC_4, m_fage_FC_5, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_fage_FC_1, m_fage_FC_2, m_fage_FC_3, 
              m_fage_FC_4, m_fage_FC_5, weights = "loo")

conditional_effects(m_fage_FC_3)



#### NO CARE 


d_age_NC <- data %>%
  select(male_F, female_F, carcass.s, female_size.s, male_size.s, 
         no_larvae, hybrid_fac, treatment, success, brood_mass,
         FEMALE_AGE, MALE_AGE) %>% 
  filter(treatment == "N",
         success == 1,
         MALE_AGE > 0,
         FEMALE_AGE > 0) 

#### MALE AGE

ggplot(aes(y = MALE_AGE, x = male_F, col = as.factor(female_F)),
       data = d_age_NC) +
  geom_point(position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2),
             alpha = 0.5) 

m_mage_NC_1 <- brm(MALE_AGE ~ female_F * male_F + 
                     carcass.s + female_size.s + male_size.s +
                     brood_mass,
                   data = d_age_NC, 
                   family = "weibull",
                   prior = set_prior('normal(0,5)'),
                   chains = 4, cores = 4,
                   control = list(adapt_delta = 0.99, max_treedepth = 15),
                   save_all_pars = TRUE)

summary(m_mage_NC_1)
pp_check(m_mage_NC_1)

m_mage_NC_2 <- update(m_mage_NC_1,
                      formula. = MALE_AGE ~ male_F + female_F + 
                        carcass.s + female_size.s + male_size.s + brood_mass)

m_mage_NC_3 <- update(m_mage_NC_1,
                      formula. = MALE_AGE ~ female_F + 
                        carcass.s + female_size.s + male_size.s + brood_mass)

m_mage_NC_4 <- update(m_mage_NC_1,
                      formula. = MALE_AGE ~ male_F + 
                        carcass.s + female_size.s + male_size.s + brood_mass)

m_mage_NC_5 <- update(m_mage_NC_1,
                      formula. = MALE_AGE ~ carcass.s + female_size.s + male_size.s + 
                        brood_mass)

m_mage_NC_1 <- add_criterion(m_mage_NC_1, c("waic", "loo"))
m_mage_NC_2 <- add_criterion(m_mage_NC_2, c("waic", "loo"))
m_mage_NC_3 <- add_criterion(m_mage_NC_3, c("waic", "loo"))
m_mage_NC_4 <- add_criterion(m_mage_NC_4, c("waic", "loo"))
m_mage_NC_5 <- add_criterion(m_mage_NC_5, c("waic", "loo"))

w <- loo_compare(m_mage_NC_1, m_mage_NC_2, m_mage_NC_3, 
                 m_mage_NC_4, m_mage_NC_5, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_mage_NC_1, m_mage_NC_2, m_mage_NC_3, 
              m_mage_NC_4, m_mage_NC_5, weights = "loo")

conditional_effects(m_mage_NC_2)
summary(m_mage_NC_2)


### FEMALE LIFESPAN


m_fage_NC_1 <- brm(FEMALE_AGE ~ female_F * male_F + 
                     carcass.s + female_size.s + male_size.s +
                     brood_mass,
                   data = d_age_NC, 
                   family = "weibull",
                   prior = set_prior('normal(0,5)'),
                   chains = 4, cores = 4,
                   control = list(adapt_delta = 0.99, max_treedepth = 15),
                   save_all_pars = TRUE)

summary(m_fage_NC_1)
pp_check(m_fage_NC_1)

m_fage_NC_2 <- update(m_fage_NC_1,
                      formula. = FEMALE_AGE ~ male_F + female_F + 
                        carcass.s + female_size.s + male_size.s + brood_mass)

m_fage_NC_3 <- update(m_fage_NC_1,
                      formula. = FEMALE_AGE ~ female_F + 
                        carcass.s + female_size.s + male_size.s + brood_mass)

m_fage_NC_4 <- update(m_fage_NC_1,
                      formula. = FEMALE_AGE ~ male_F + 
                        carcass.s + female_size.s + male_size.s + brood_mass)

m_fage_NC_5 <- update(m_fage_NC_1,
                      formula. = FEMALE_AGE ~ carcass.s + female_size.s + male_size.s + 
                        brood_mass)

m_fage_NC_1 <- add_criterion(m_fage_NC_1, c("waic", "loo"))
m_fage_NC_2 <- add_criterion(m_fage_NC_2, c("waic", "loo"))
m_fage_NC_3 <- add_criterion(m_fage_NC_3, c("waic", "loo"))
m_fage_NC_4 <- add_criterion(m_fage_NC_4, c("waic", "loo"))
m_fage_NC_5 <- add_criterion(m_fage_NC_5, c("waic", "loo"))

w <- loo_compare(m_fage_NC_1, m_fage_NC_2, m_fage_NC_3, 
                 m_fage_NC_4, m_fage_NC_5, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_fage_NC_1, m_fage_NC_2, m_fage_NC_3, 
              m_fage_NC_4, m_fage_NC_5, weights = "loo")

summary(m_fage_NC_3)
pp_check(m_fage_NC_3)
conditional_effects(m_fage_NC_3)


ggplot(aes(x = MALE_AGE, y = FEMALE_AGE, 
           col = as.factor(female_F), shape = as.factor(male_F)),
       data = d_age_NC) +
  geom_point() +
  geom_smooth(method = "lm")
