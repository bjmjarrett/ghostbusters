library(tidyverse)
library(gridExtra) 
library(brms)
library(gameofthrones) 
library(modelr) 
library(tidybayes) 

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

# data for the hole analysis
d_hole <- d %>%
  select(male_F, female_F, carcass.s, female_size.s, male_size.s, hole, hybrid_fac)

# plot hole data
ggplot(aes(y = hole, x = male_F, col = as.factor(female_F)),
       data = d_hole) +
  geom_point(position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2),
             alpha = 0.5) 

# model with hybrid factor
m_hyb_hole_1 <- brm(hole ~ hybrid_fac + carcass.s + female_size.s + male_size.s,
                    data = d_hole, 
                    family = "bernoulli",
                    prior = set_prior('normal(0,5)'),
                    chains = 4, cores = 4,
                    control = list(adapt_delta = 0.99, max_treedepth = 15),
                    save_all_pars = TRUE)

summary(m_hyb_hole_1)
pp_check(m_hyb_hole_1)

# contrasts for the three levels of the hybrid factor
hyp <- hypothesis(m_hyb_hole_1, c("Intercept = Intercept + hybrid_fachybrid", "Intercept = Intercept + hybrid_facno_no", "Intercept + hybrid_facno_no = Intercept + hybrid_fachybrid"), alpha = 0.05)

plot(hyp)
hyp

# simpler model
m_hyb_hole_2 <- update(m_hyb_hole_1,
                       formula. = hole ~ carcass.s + female_size.s + male_size.s)
summary(m_hyb_hole_2)

# model comparisons
m_hyb_hole_1 <- add_criterion(m_hyb_hole_1, c("waic", "loo"))
m_hyb_hole_2 <- add_criterion(m_hyb_hole_2, c("waic", "loo"))

w <- loo_compare(m_hyb_hole_1, m_hyb_hole_2, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_hyb_hole_1, m_hyb_hole_2, weights = "loo")


#### MODEL with male and female background

# full model
m_hole_1 <- brm(hole ~ female_F * male_F + carcass.s + female_size.s + male_size.s,
                data = d_hole, 
                family = "bernoulli",
                prior = set_prior('normal(0,5)'),
                chains = 4, cores = 4,
                control = list(adapt_delta = 0.99, max_treedepth = 15),
                save_all_pars = TRUE)

summary(m_hole_1)
pp_check(m_hole_1)

m_hole_2 <- update(m_hole_1,
                   formula. = hole ~ female_F + male_F + 
                     carcass.s + female_size.s + male_size.s)

m_hole_3 <- update(m_hole_1,
                   formula. = hole ~ female_F + 
                     carcass.s + female_size.s + male_size.s)

m_hole_4 <- update(m_hole_1,
                   formula. = hole ~ male_F + 
                     carcass.s + female_size.s + male_size.s)

m_hole_5 <- update(m_hole_1,
                   formula. = hole ~ carcass.s + female_size.s + male_size.s)

# model comparisons
m_hole_1 <- add_criterion(m_hole_1, c("waic", "loo"))
m_hole_2 <- add_criterion(m_hole_2, c("waic", "loo"))
m_hole_3 <- add_criterion(m_hole_3, c("waic", "loo"))
m_hole_4 <- add_criterion(m_hole_4, c("waic", "loo"))
m_hole_5 <- add_criterion(m_hole_5, c("waic", "loo"))

w <- loo_compare(m_hole_1, m_hole_2, m_hole_3, m_hole_4, m_hole_5, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_hole_1, m_hole_2, m_hole_3, m_hole_4, m_hole_5, weights = "loo")


#### PLOT THE EFFECTS

# average values of covariates for model prediction to plot the effects we are interested in
Avg.carcass <- mean(d_hole$carcass.s)
Avg.male <- mean(d_hole$male_size.s)
Avg.female <- mean(d_hole$female_size.s)

# get predictions from the interaction model using an empty dataframe. Uses functions from todybayes and modelr
post2 <- d_hole %>%
  data_grid(carcass.s = Avg.carcass,
            male_size.s = Avg.male,
            female_size.s = Avg.female,
            male_F = c(0, 1),
            female_F = seq(0, 1, 0.1)) %>%
  add_fitted_draws(model = m_hole_1,
                   re_formula = NA,
                   allow_new_levels = FALSE)

inters <- post2 %>%
  mean_hdci(.value, .width = 0.95)

# plot the predicted values and 95% ribbon
ggplot(post2, aes(y = .value, x = female_F)) +
  geom_smooth(se = F, aes(col = as.factor(male_F))) +
  geom_ribbon(data = inters,
              aes(ymin = .lower, ymax = .upper,
                  fill = as.factor(male_F)), alpha = 0.1) +
  theme_linedraw() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = c(0.95, 0.95),
        legend.justification = c(1,1)
  ) +
  scale_y_continuous("Proportion of broods with incision", limits = c(0, 0.6)) +
  scale_x_continuous("Female background", breaks = c(0,1), 
                     labels = c("No Care", "Full Care"), limits = c(-0.1, 1.1)) +
  guides(col = guide_legend(nrow = 1)) +
  scale_color_got(name = "Male background",
                  labels = c("No Care", "Full Care"),
                  discrete = TRUE, option = "Tully", direction = 1) +
  scale_fill_got_d(name = "Male background",
                   labels = c("No Care", "Full Care"), option = "Tully", direction = 1) 

