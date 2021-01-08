library(tidyverse)
library(rethinking)
library(gridExtra)
library(brms)
library(gameofthrones)

###
### HETEROSIS
###

# check for heterosis
heterosis <- read.csv("heterosis_check.csv")

ggplot(heterosis,
       aes(y = no_larvae, x = experiment)) +
  geom_boxplot() +
  facet_grid(.~care)

levels(heterosis$care) <- c("Full Care", "No Care")

p1 <- heterosis %>% 
  ggplot(aes(y = no_larvae, x = experiment, col = care)) +
  geom_boxplot() +
  facet_grid(.~care) + 
  ylab("Number of larvae") +
  xlab("Cross") +
  scale_x_discrete(labels=c("Within-line", "Between-line")) +
  theme_linedraw() +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  ) +
  scale_color_got(discrete = TRUE, option = "Tully", direction = -1) +
  labs(tag = "A")


p2 <- heterosis %>% 
  ggplot(aes(y = brood_mass/no_larvae, x = experiment, col = care)) +
  geom_boxplot() +
  facet_grid(.~care) +
  ylab("Average larval mass (g)") +
  xlab("Cross") +
  scale_x_discrete(labels=c("Within-line", "Between-line")) +
  theme_linedraw() +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
  ) +
  scale_color_got(discrete = TRUE, option = "Tully", direction = -1) +
  labs(tag = "B")


### heterosis plot

grid.arrange(p1, p2, nrow = 2)

### analysis

heterosis$carcass.s <- (heterosis$carcass - mean(heterosis$carcass)) / sd(heterosis$carcass)

# heterosis model 1: FULL CARE

d_h_1 <- heterosis %>% filter(care == "Full Care")
hist(d_h_1$no_larvae)

m_h_1.1 <- brm(no_larvae ~ experiment + carcass.s,
               data = d_h_1, 
               family = "zero_inflated_negbinomial",
               prior = set_prior('normal(0,5)'),
               chains = 4, cores = 4,
               control = list(adapt_delta = 0.99, max_treedepth = 15),
               save_all_pars = TRUE)

summary(m_h_1.1)
pp_check(m_h_1.1)


m_h_1.2 <- update(m_h_1.1,
                  formula. = no_larvae ~ carcass.s)

m_h_1.1 <- add_criterion(m_h_1.1, c("waic", "loo"))
m_h_1.2 <- add_criterion(m_h_1.2, c("waic", "loo"))

w <- loo_compare(m_h_1.1, m_h_1.2, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_h_1.1, m_h_1.2, weights = "loo")

# heterosis model 2: FULL CARE

d_h_1 <- heterosis %>% filter(care == "Full Care")
hist(d_h_1$no_larvae)

m_h_2.1 <- brm(brood_mass/no_larvae ~ experiment + carcass.s,
               data = d_h_1, 
               family = "gaussian",
               prior = set_prior('normal(0,5)'),
               chains = 4, cores = 4,
               control = list(adapt_delta = 0.99, max_treedepth = 15),
               save_all_pars = TRUE)

summary(m_h_2.1)
pp_check(m_h_2.1)


m_h_2.2 <- update(m_h_2.1,
                  formula. = brood_mass/no_larvae ~ carcass.s)

m_h_2.1 <- add_criterion(m_h_2.1, c("waic", "loo"))
m_h_2.2 <- add_criterion(m_h_2.2, c("waic", "loo"))

w <- loo_compare(m_h_2.1, m_h_2.2, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_h_2.1, m_h_2.2, weights = "loo")

# heterosis model 3: NO CARE

d_h_2 <- heterosis %>% filter(care == "No Care")
hist(d_h_2$no_larvae)

m_h_3.1 <- brm(no_larvae ~ experiment + carcass.s,
               data = d_h_2, 
               family = "zero_inflated_negbinomial",
               prior = set_prior('normal(0,5)'),
               chains = 4, cores = 4,
               control = list(adapt_delta = 0.99, max_treedepth = 15),
               save_all_pars = TRUE)

summary(m_h_3.1)
pp_check(m_h_3.1)

m_h_3.2 <- update(m_h_3.1,
                  formula. = no_larvae ~ carcass.s)

m_h_3.1 <- add_criterion(m_h_3.1, c("waic", "loo"))
m_h_3.2 <- add_criterion(m_h_3.2, c("waic", "loo"))

w <- loo_compare(m_h_3.1, m_h_3.2, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_h_3.1, m_h_3.2, weights = "loo")

# heterosis model 4: NO CARE

d_h_2 <- heterosis %>% filter(care == "No Care")
hist(d_h_2$no_larvae)

m_h_4.1 <- brm(brood_mass/no_larvae ~ experiment + carcass.s,
               data = d_h_2, 
               family = "gaussian",
               prior = set_prior('normal(0,5)'),
               chains = 4, cores = 4,
               control = list(adapt_delta = 0.99, max_treedepth = 15),
               save_all_pars = TRUE)

summary(m_h_4.1)
pp_check(m_h_4.1)

m_h_4.2 <- update(m_h_4.1,
                  formula. = brood_mass/no_larvae ~ carcass.s)

m_h_4.1 <- add_criterion(m_h_4.1, c("waic", "loo"))
m_h_4.2 <- add_criterion(m_h_4.2, c("waic", "loo"))

w <- loo_compare(m_h_4.1, m_h_4.2, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_h_4.1, m_h_4.2, weights = "loo")