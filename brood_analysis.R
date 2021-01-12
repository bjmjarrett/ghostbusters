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

#### MODELS

# hybrid factor in FULL CARE

d_brood_FC <- data %>%
  select(male_F, female_F, carcass.s, female_size.s, male_size.s, 
         no_larvae, hybrid_fac, treatment, hole, success) %>% 
  filter(treatment == "F")


ggplot(aes(y = no_larvae, x = male_F, col = as.factor(female_F)),
       data = d_brood_FC) +
  geom_point(position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2),
             alpha = 0.5) 

m_hyb_brood_fc_1 <- brm(no_larvae ~ hybrid_fac + hole + 
                          carcass.s + female_size.s + male_size.s,
                        data = d_brood_FC, 
                        family = "gaussian",
                        prior = set_prior('normal(0,5)'),
                        chains = 4, cores = 4,
                        control = list(adapt_delta = 0.99, max_treedepth = 15),
                        save_all_pars = TRUE)

summary(m_hyb_brood_fc_1)
pp_check(m_hyb_brood_fc_1)

hyp <- hypothesis(m_hyb_brood_fc_1, c("Intercept = Intercept + hybrid_fachybrid", 
                                      "Intercept = Intercept + hybrid_facno_no", 
                                      "Intercept + hybrid_facno_no = Intercept + hybrid_fachybrid"), 
                  alpha = 0.05)

plot(hyp)
hyp


m_hyb_brood_fc_2 <- update(m_hyb_brood_fc_1,
                           formula. = no_larvae ~ hole + 
                             carcass.s + female_size.s + male_size.s)
summary(m_hyb_brood_fc_2)

m_hyb_brood_fc_1 <- add_criterion(m_hyb_brood_fc_1, c("waic", "loo"))
m_hyb_brood_fc_2 <- add_criterion(m_hyb_brood_fc_2, c("waic", "loo"))

w <- loo_compare(m_hyb_brood_fc_1, m_hyb_brood_fc_2, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_hyb_brood_fc_1, m_hyb_brood_fc_2, weights = "loo")


#####

m_brood_fc_1 <- brm(no_larvae ~ male_F * female_F + hole +
                      carcass.s + female_size.s + male_size.s,
                    data = d_brood_FC, 
                    family = "gaussian",
                    prior = set_prior('normal(0,5)'),
                    chains = 4, cores = 4,
                    control = list(adapt_delta = 0.99, max_treedepth = 15),
                    save_all_pars = TRUE)

summary(m_brood_fc_1)
pp_check(m_brood_fc_1)

m_brood_fc_2 <- update(m_brood_fc_1,
                       formula. = no_larvae ~ male_F + female_F + hole + carcass.s + 
                         female_size.s + male_size.s)

m_brood_fc_3 <- update(m_brood_fc_1,
                       formula. = no_larvae ~ female_F + hole + carcass.s + 
                         female_size.s + male_size.s)

m_brood_fc_4 <- update(m_brood_fc_1,
                       formula. = no_larvae ~ male_F + hole + carcass.s + 
                         female_size.s + male_size.s)

m_brood_fc_5 <- update(m_brood_fc_1,
                       formula. = no_larvae ~ carcass.s + hole +
                         female_size.s + male_size.s)

m_brood_fc_1 <- add_criterion(m_brood_fc_1, c("waic", "loo"))
m_brood_fc_2 <- add_criterion(m_brood_fc_2, c("waic", "loo"))
m_brood_fc_3 <- add_criterion(m_brood_fc_3, c("waic", "loo"))
m_brood_fc_4 <- add_criterion(m_brood_fc_4, c("waic", "loo"))
m_brood_fc_5 <- add_criterion(m_brood_fc_5, c("waic", "loo"))

w <- loo_compare(m_brood_fc_1, m_brood_fc_2, m_brood_fc_3, m_brood_fc_4, m_brood_fc_5, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_brood_fc_1, m_brood_fc_2, m_brood_fc_3, m_brood_fc_4, m_brood_fc_5, weights = "loo")


#### NO CARE ENVIRONMENT

d_brood_NC <- data %>%
  select(male_F, female_F, carcass.s, female_size.s, male_size.s, 
         no_larvae, hybrid_fac, treatment, hole, success) %>% 
  filter(treatment == "N")

ggplot(aes(y = no_larvae, x = male_F, col = as.factor(female_F)),
       data = d_brood_NC) +
  geom_point(position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2),
             alpha = 0.5) 

### hybrid factor on successes and failures
### including hole
m_hyb_brood_1 <- brm(success ~ hybrid_fac + hole + carcass.s + female_size.s + male_size.s,
                     data = d_brood_NC, 
                     family = "bernoulli",
                     prior = set_prior('normal(0,5)'),
                     chains = 4, cores = 4,
                     control = list(adapt_delta = 0.99, max_treedepth = 15),
                     save_all_pars = TRUE)

summary(m_hyb_brood_1)
pp_check(m_hyb_brood_1)

hyp <- hypothesis(m_hyb_brood_1, c("Intercept = Intercept + hybrid_fachybrid", "Intercept = Intercept + hybrid_facno_no", "Intercept + hybrid_facno_no = Intercept + hybrid_fachybrid"), alpha = 0.05)

plot(hyp)
hyp


m_hyb_brood_2 <- update(m_hyb_brood_1,
                        formula. = success ~ hole + carcass.s + female_size.s + male_size.s)
summary(m_hyb_brood_2)

m_hyb_brood_1 <- add_criterion(m_hyb_brood_1, c("waic", "loo"))
m_hyb_brood_2 <- add_criterion(m_hyb_brood_2, c("waic", "loo"))

w <- loo_compare(m_hyb_brood_1, m_hyb_brood_2, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_hyb_brood_1, m_hyb_brood_2, weights = "loo")


#### male and female background
m_brood_1 <- brm(success ~ male_F * female_F + hole +
                      carcass.s + female_size.s + male_size.s,
                 data = d_brood_NC, 
                 family = "bernoulli",
                 prior = set_prior('normal(0,5)'),
                 chains = 4, cores = 4,
                 control = list(adapt_delta = 0.99, max_treedepth = 15),
                 save_all_pars = TRUE)

summary(m_brood_1)
pp_check(m_brood_1)

m_brood_2 <- update(m_brood_1,
                    formula. = success ~ male_F + female_F + hole +
                      carcass.s + female_size.s + male_size.s)

m_brood_3 <- update(m_brood_1,
                    formula. = success ~ female_F + hole +
                      carcass.s + female_size.s + male_size.s)

m_brood_4 <- update(m_brood_1,
                    formula. = success ~ male_F + hole +
                      carcass.s + female_size.s + male_size.s)

m_brood_5 <- update(m_brood_1,
                    formula. = success ~ hole +
                      carcass.s + female_size.s + male_size.s)

m_brood_1 <- add_criterion(m_brood_1, c("waic", "loo"))
m_brood_2 <- add_criterion(m_brood_2, c("waic", "loo"))
m_brood_3 <- add_criterion(m_brood_3, c("waic", "loo"))
m_brood_4 <- add_criterion(m_brood_4, c("waic", "loo"))
m_brood_5 <- add_criterion(m_brood_5, c("waic", "loo"))

w <- loo_compare(m_brood_1, m_brood_2, m_brood_3, m_brood_4, m_brood_5, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_brood_1, m_brood_2, m_brood_3, m_brood_4, m_brood_5, weights = "loo")

### hybrid factor on successes and failures
### EXCLUDING hole
m_hyb_brood_1_nohole <- brm(success ~ hybrid_fac + carcass.s + female_size.s + male_size.s,
                     data = d_brood_NC, 
                     family = "bernoulli",
                     prior = set_prior('normal(0,5)'),
                     chains = 4, cores = 4,
                     control = list(adapt_delta = 0.99, max_treedepth = 15),
                     save_all_pars = TRUE)

summary(m_hyb_brood_1_nohole)
pp_check(m_hyb_brood_1_nohole)

hyp <- hypothesis(m_hyb_brood_1_nohole, c("Intercept = Intercept + hybrid_fachybrid", 
                                          "Intercept = Intercept + hybrid_facno_no", 
                                          "Intercept + hybrid_facno_no = Intercept + hybrid_fachybrid"),
                  alpha = 0.05)

plot(hyp)
hyp


m_hyb_brood_2_nohole <- update(m_hyb_brood_1_nohole,
                        formula. = success ~ carcass.s + female_size.s + male_size.s)
summary(m_hyb_brood_2_nohole)

m_hyb_brood_1_nohole <- add_criterion(m_hyb_brood_1_nohole, c("waic", "loo"))
m_hyb_brood_2_nohole <- add_criterion(m_hyb_brood_2_nohole, c("waic", "loo"))

w <- loo_compare(m_hyb_brood_1_nohole, m_hyb_brood_2_nohole, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_hyb_brood_1_nohole, m_hyb_brood_2_nohole, weights = "loo")


#### male and female background
m_brood_1_nohole <- brm(success ~ male_F * female_F +
                   carcass.s + female_size.s + male_size.s,
                 data = d_brood_NC, 
                 family = "bernoulli",
                 prior = set_prior('normal(0,5)'),
                 chains = 4, cores = 4,
                 control = list(adapt_delta = 0.99, max_treedepth = 15),
                 save_all_pars = TRUE)

summary(m_brood_1_nohole)
pp_check(m_brood_1_nohole)

m_brood_2_nohole <- update(m_brood_1_nohole,
                    formula. = success ~ male_F + female_F +
                      carcass.s + female_size.s + male_size.s)

m_brood_3_nohole <- update(m_brood_1_nohole,
                    formula. = success ~ female_F +
                      carcass.s + female_size.s + male_size.s)

m_brood_4_nohole <- update(m_brood_1_nohole,
                    formula. = success ~ male_F +
                      carcass.s + female_size.s + male_size.s)

m_brood_5_nohole <- update(m_brood_1_nohole,
                    formula. = success ~
                      carcass.s + female_size.s + male_size.s)

m_brood_1_nohole <- add_criterion(m_brood_1_nohole, c("waic", "loo"))
m_brood_2_nohole <- add_criterion(m_brood_2_nohole, c("waic", "loo"))
m_brood_3_nohole <- add_criterion(m_brood_3_nohole, c("waic", "loo"))
m_brood_4_nohole <- add_criterion(m_brood_4_nohole, c("waic", "loo"))
m_brood_5_nohole <- add_criterion(m_brood_5_nohole, c("waic", "loo"))

w <- loo_compare(m_brood_1_nohole, m_brood_2_nohole, m_brood_3_nohole, 
                 m_brood_4_nohole, m_brood_5_nohole, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_brood_1_nohole, m_brood_2_nohole, m_brood_3_nohole, 
              m_brood_4_nohole, m_brood_5_nohole, weights = "loo")


### BROOD SIZE ANALYSIS with only successful broods

d_brood_NCs <- d %>%
  select(male_F, female_F, carcass.s, female_size.s, male_size.s, 
         no_larvae, hybrid_fac, treatment, hole, success) %>% 
  filter(treatment == "N",
         success == 1)


ggplot(aes(y = no_larvae, x = male_F, col = as.factor(female_F)),
       data = d_brood_NCs) +
  geom_point(position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2),
             alpha = 0.5) 

m_hyb_brood_1ncs <- brm(no_larvae ~ hybrid_fac + carcass.s + female_size.s + male_size.s + hole,
                        data = d_brood_NCs, 
                        family = "negbinomial",
                        prior = set_prior('normal(0,5)'),
                        chains = 4, cores = 4,
                        control = list(adapt_delta = 0.99, max_treedepth = 15),
                        save_all_pars = TRUE)

summary(m_hyb_brood_1ncs)
pp_check(m_hyb_brood_1ncs)

hyp <- hypothesis(m_hyb_brood_1ncs, c("Intercept = Intercept + hybrid_fachybrid", 
                                      "Intercept = Intercept + hybrid_facno_no", 
                                      "Intercept + hybrid_facno_no = Intercept + hybrid_fachybrid"),
                  alpha = 0.05)

plot(hyp)
hyp


m_hyb_brood_2ncs <- update(m_hyb_brood_1ncs,
                           formula. = no_larvae ~ carcass.s + female_size.s + male_size.s + hole)
summary(m_hyb_brood_2ncs)

m_hyb_brood_1ncs <- add_criterion(m_hyb_brood_1ncs, c("waic", "loo"))
m_hyb_brood_2ncs <- add_criterion(m_hyb_brood_2ncs, c("waic", "loo"))

w <- loo_compare(m_hyb_brood_1ncs, m_hyb_brood_2ncs, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_hyb_brood_1ncs, m_hyb_brood_2ncs, weights = "loo")

#####

m_brood_1ncs <- brm(no_larvae ~ male_F * female_F + hole +
                      carcass.s + female_size.s + male_size.s,
                    data = d_brood_NCs, 
                    family = "negbinomial",
                    prior = set_prior('normal(0,5)'),
                    chains = 4, cores = 4,
                    control = list(adapt_delta = 0.99, max_treedepth = 15),
                    save_all_pars = TRUE)

summary(m_brood_1ncs)
pp_check(m_brood_1ncs)

m_brood_2ncs <- update(m_brood_1ncs,
                       formula. = no_larvae ~ male_F + female_F + hole + carcass.s + 
                         female_size.s + male_size.s)

m_brood_3ncs <- update(m_brood_1ncs,
                       formula. = no_larvae ~ female_F + hole + carcass.s + 
                         female_size.s + male_size.s)

m_brood_4ncs <- update(m_brood_1ncs,
                       formula. = no_larvae ~ male_F + hole + carcass.s + 
                         female_size.s + male_size.s)

m_brood_5ncs <- update(m_brood_1ncs,
                       formula. = no_larvae ~ carcass.s + hole +
                         female_size.s + male_size.s)

m_brood_1ncs <- add_criterion(m_brood_1ncs, c("waic", "loo"))
m_brood_2ncs <- add_criterion(m_brood_2ncs, c("waic", "loo"))
m_brood_3ncs <- add_criterion(m_brood_3ncs, c("waic", "loo"))
m_brood_4ncs <- add_criterion(m_brood_4ncs, c("waic", "loo"))
m_brood_5ncs <- add_criterion(m_brood_5ncs, c("waic", "loo"))

w <- loo_compare(m_brood_1ncs, m_brood_2ncs, m_brood_3ncs, 
                 m_brood_4ncs, m_brood_5ncs, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_brood_1ncs, m_brood_2ncs, m_brood_3ncs, 
              m_brood_4ncs, m_brood_5ncs, weights = "loo")

### plot
m_brood_3ncs_plot <- brm(no_larvae ~ as.factor(female_F) + hole +
                      carcass.s + female_size.s + male_size.s,
                    data = d_brood_NCs, 
                    family = "negbinomial",
                    prior = set_prior('normal(0,5)'),
                    chains = 4, cores = 4,
                    control = list(adapt_delta = 0.99, max_treedepth = 15),
                    save_all_pars = TRUE)

Avg.carcass <- mean(d_brood_NCs$carcass.s)
Avg.male <- mean(d_brood_NCs$male_size.s)
Avg.female <- mean(d_brood_NCs$female_size.s)

post2 <- d_brood_NCs %>%
  data_grid(carcass.s = Avg.carcass,
            male_size.s = Avg.male,
            female_size.s = Avg.female,
            hole = 1,
            female_F = c(0, 1)) %>%
  add_fitted_draws(model = m_brood_3ncs_plot,
                   re_formula = NA,
                   allow_new_levels = FALSE)

inters <- post2 %>%
  mean_hdci(.value, .width = 0.95)

# plot using predicted values from the model as boxplots
ggplot(post2, aes(y = .value, x = as.factor(female_F), col = as.factor(female_F))) +
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(data = d_brood_NCs, alpha = 0.3, 
              aes(y = no_larvae, x = as.factor(female_F), col = as.factor(female_F))) +
  scale_color_got(discrete = TRUE, option = "Tully", direction = 1) +
  theme_linedraw() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Number of larvae") +
  scale_x_discrete("Female background", breaks = c(0,1), 
                   labels = c("No Care", "Full Care"))

# plot without model predictions
ggplot(d_brood_NCs, aes(y = no_larvae, x = as.factor(female_F), col = as.factor(female_F))) +
  geom_boxplot() + 
  geom_jitter(alpha = 0.3, aes(col = as.factor(female_F))) + 
  scale_color_got(discrete = TRUE, option = "Tully", direction = 1) +
  theme_linedraw() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Number of larvae") +
  scale_x_discrete("Female background", breaks = c(0,1), 
                   labels = c("No Care", "Full Care"))
