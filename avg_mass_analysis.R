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

##### AVERAGE LARVAL MASS
# FULL CARE

d_avg_FC <- data %>%
  select(male_F, female_F, carcass.s, female_size.s, male_size.s, 
         no_larvae, hybrid_fac, treatment, hole, success, brood_mass) %>% 
  filter(treatment == "F",
         success == 1) %>% 
  mutate(avg_mass = brood_mass / no_larvae)

ggplot(aes(y = avg_mass, x = no_larvae),
       data = d_avg_FC) +
  geom_point()

ggplot(aes(y = avg_mass, x = male_F, col = as.factor(female_F)),
       data = d_avg_FC) +
  geom_point(position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2),
             alpha = 0.5) 

m_hyb_avg_1 <- brm(avg_mass ~ hybrid_fac + no_larvae + carcass.s + 
                     female_size.s + male_size.s,
                   data = d_avg_FC, 
                   family = "gaussian",
                   prior = set_prior('normal(0,5)'),
                   chains = 1, cores = 4,
                   control = list(adapt_delta = 0.99, max_treedepth = 15),
                   save_all_pars = TRUE)

summary(m_hyb_avg_1)
pp_check(m_hyb_avg_1)

hyp <- hypothesis(m_hyb_brood_1, c("Intercept = Intercept + hybrid_fachybrid", 
                                   "Intercept = Intercept + hybrid_facno_no", 
                                   "Intercept + hybrid_facno_no = Intercept + hybrid_fachybrid"),
                  alpha = 0.05)

plot(hyp)
hyp


m_hyb_brood_2 <- update(m_hyb_brood_1,
                        formula. = avg_mass ~ no_larvae + carcass.s + female_size.s + male_size.s)
summary(m_hyb_brood_2)

m_hyb_brood_1 <- add_criterion(m_hyb_brood_1, c("waic", "loo"))
m_hyb_brood_2 <- add_criterion(m_hyb_brood_2, c("waic", "loo"))

w <- loo_compare(m_hyb_brood_1, m_hyb_brood_2, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_hyb_brood_1, m_hyb_brood_2, weights = "loo")

####

m_avg_fc_1 <- brm(avg_mass ~ female_F * male_F + no_larvae + carcass.s + 
                    female_size.s + male_size.s,
                  data = d_avg_FC, 
                  family = "gaussian",
                  prior = set_prior('normal(0,5)'),
                  chains = 4, cores = 4,
                  control = list(adapt_delta = 0.99, max_treedepth = 15),
                  save_all_pars = TRUE)

summary(m_avg_fc_1)
pp_check(m_avg_fc_1)


m_avg_fc_2 <- update(m_avg_fc_1,
                     formula. = avg_mass ~ male_F + female_F + no_larvae + carcass.s + 
                       female_size.s + male_size.s)

m_avg_fc_3 <- update(m_avg_fc_1,
                     formula. = avg_mass ~ female_F + no_larvae + carcass.s + 
                       female_size.s + male_size.s)

m_avg_fc_4 <- update(m_avg_fc_1,
                     formula. = avg_mass ~ male_F + no_larvae + carcass.s + 
                       female_size.s + male_size.s)

m_avg_fc_5 <- update(m_avg_fc_1,
                     formula. = avg_mass ~ no_larvae + carcass.s +
                       female_size.s + male_size.s)

m_avg_fc_1 <- add_criterion(m_avg_fc_1, c("waic", "loo"))
m_avg_fc_2 <- add_criterion(m_avg_fc_2, c("waic", "loo"))
m_avg_fc_3 <- add_criterion(m_avg_fc_3, c("waic", "loo"))
m_avg_fc_4 <- add_criterion(m_avg_fc_4, c("waic", "loo"))
m_avg_fc_5 <- add_criterion(m_avg_fc_5, c("waic", "loo"))

w <- loo_compare(m_avg_fc_1, m_avg_fc_2, m_avg_fc_3, 
                 m_avg_fc_4, m_avg_fc_5, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_avg_fc_1, m_avg_fc_2, m_avg_fc_3, 
              m_avg_fc_4, m_avg_fc_5, weights = "loo")


#### NO CARE ENVIRONMENT

# subset for the No Care data
d_avg_NC <- data %>%
  select(male_F, female_F, carcass.s, female_size.s, male_size.s, 
         no_larvae, hybrid_fac, treatment, hole, success, brood_mass) %>% 
  filter(treatment == "N",
         success == 1) %>% 
  mutate(avg_mass = brood_mass / no_larvae)

ggplot(aes(y = avg_mass, x = no_larvae),
       data = d_avg_NC) +
  geom_point()

ggplot(aes(y = avg_mass, x = male_F, col = as.factor(female_F)),
       data = d_avg_NC) +
  geom_point(position = position_jitterdodge(dodge.width = 0.4, jitter.width = 0.2),
             alpha = 0.5) 

#### MODEL TO EXTRACT RESIDUALS

# the relationship between average larval mass and brood size in the No Care environment is best explained by a cubic relationship (Schrader, Jarrett and Kilner 2015 Evolution). Including this in the models was causing the chains to not converge so I used the residuals from the below model. The residuals are pasted at the bottom on this file to save running this model

mod_avg_res <- brm(avg_mass ~ no_larvae + I(no_larvae^2) + I(no_larvae^3),
                   data = d_avg_NC,
                   family = "gaussian",
                   prior = set_prior('normal(0,1)'),
                   chains = 4, cores = 4,
                   control = list(adapt_delta = 0.99, max_treedepth = 20),
                   save_all_pars = TRUE,
                   iter = 20000, warmup = 2000)
summary(mod_avg_res)

# extract residuals
avg_mass_residuals <- residuals(mod_avg_res)

# append residuals to the NC dataset
d_avg_NC$avg_mass_res <- avg_mass_residuals[,1]

####

m_hyb_avg_nc_1 <- brm(avg_mass_res ~ hybrid_fac + 
                        carcass.s + female_size.s + male_size.s,
                      data = d_avg_NC, 
                      family = "gaussian",
                      prior = set_prior('normal(0,5)'),
                      chains = 4, cores = 4,
                      control = list(adapt_delta = 0.99, max_treedepth = 15),
                      save_all_pars = TRUE)

summary(m_hyb_avg_nc_1)
pp_check(m_hyb_avg_nc_1)

hyp <- hypothesis(m_hyb_avg_nc_1, c("Intercept = Intercept + hybrid_fachybrid", 
                                    "Intercept = Intercept + hybrid_facno_no", 
                                    "Intercept + hybrid_facno_no = Intercept + hybrid_fachybrid"),
                  alpha = 0.05)

plot(hyp)
hyp


m_hyb_avg_nc_2 <- update(m_hyb_avg_nc_1,
                         formula. = avg_mass_res ~ carcass.s + female_size.s + male_size.s)
summary(m_hyb_avg_nc_2)

m_hyb_avg_nc_1 <- add_criterion(m_hyb_avg_nc_1, c("waic", "loo"))
m_hyb_avg_nc_2 <- add_criterion(m_hyb_avg_nc_2, c("waic", "loo"))

w <- loo_compare(m_hyb_avg_nc_1, m_hyb_avg_nc_2, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_hyb_avg_nc_1, m_hyb_avg_nc_2, weights = "loo")

#### hybrid plot
Avg.carcass <- mean(d_avg_NC$carcass.s)
Avg.male <- mean(d_avg_NC$male_size.s)
Avg.female <- mean(d_avg_NC$female_size.s)

post2 <- d_avg_NC %>%
  data_grid(carcass.s = Avg.carcass,
            male_size.s = Avg.male,
            female_size.s = Avg.female,
            hybrid_fac = c("full_full", "hybrid", "no_no")) %>%
  add_fitted_draws(model = m_hyb_avg_nc_1,
                   re_formula = NA,
                   allow_new_levels = FALSE)

inters <- post2 %>%
  mean_hdci(.value, .width = 0.95)

# plot using predicted values from the model as boxplots
p1 <- ggplot(post2, aes(y = .value, x = hybrid_fac, col = hybrid_fac)) +
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(data = d_avg_NC, alpha = 0.3, 
              aes(y = avg_mass_res, x = hybrid_fac, col = hybrid_fac)) +
  scale_color_got(discrete = TRUE, option = "Tully", direction = -1, begin = 0) +
  theme_linedraw() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Residual average larval mass") +
  scale_x_discrete("Cross", 
                   labels = c("Full Care x Full Care", "Hybrid", "No Care x No Care")) +
  labs(tag = "A")

###

#####

m_avg_nc_1 <- brm(avg_mass_res ~ female_F * male_F + 
                    carcass.s + female_size.s + male_size.s,
                  data = d_avg_NC, 
                  family = "gaussian",
                  prior = set_prior('normal(0,5)'),
                  chains = 4, cores = 4,
                  control = list(adapt_delta = 0.99, max_treedepth = 15),
                  save_all_pars = TRUE)

summary(m_avg_nc_1)
pp_check(m_avg_nc_1)


m_avg_nc_2 <- update(m_avg_nc_1,
                     formula. = avg_mass_res ~ male_F + female_F + 
                       carcass.s + female_size.s + male_size.s)

m_avg_nc_3 <- update(m_avg_nc_1,
                     formula. = avg_mass_res ~ female_F + 
                       carcass.s + female_size.s + male_size.s)

m_avg_nc_4 <- update(m_avg_nc_1,
                     formula. = avg_mass_res ~ male_F + 
                       carcass.s + female_size.s + male_size.s)

m_avg_nc_5 <- update(m_avg_nc_1,
                     formula. = avg_mass_res ~ carcass.s + female_size.s + male_size.s)

m_avg_nc_1 <- add_criterion(m_avg_nc_1, c("waic", "loo"))
m_avg_nc_2 <- add_criterion(m_avg_nc_2, c("waic", "loo"))
m_avg_nc_3 <- add_criterion(m_avg_nc_3, c("waic", "loo"))
m_avg_nc_4 <- add_criterion(m_avg_nc_4, c("waic", "loo"))
m_avg_nc_5 <- add_criterion(m_avg_nc_5, c("waic", "loo"))

w <- loo_compare(m_avg_nc_1, m_avg_nc_2, m_avg_nc_3, 
                 m_avg_nc_4, m_avg_nc_5, criterion = "loo")
print(w, simplify = F)
cbind(waic_diff = w[,1] * -2, se = w[,2] * 2)

model_weights(m_avg_nc_1, m_avg_nc_2, m_avg_nc_3, 
              m_avg_nc_4, m_avg_nc_5, weights = "loo")


Avg.carcass <- mean(d_avg_NC$carcass.s)
Avg.male <- mean(d_avg_NC$male_size.s)
Avg.female <- mean(d_avg_NC$female_size.s)

post2 <- d_avg_NC %>%
  data_grid(carcass.s = Avg.carcass,
            male_size.s = Avg.male,
            female_size.s = Avg.female,
            female_F = c(0, 1)) %>%
  add_fitted_draws(model = m_avg_nc_3,
                   re_formula = NA,
                   allow_new_levels = FALSE)

inters <- post2 %>%
  mean_hdci(.value, .width = 0.95)

# plot using predicted values from the model as boxplots
p2 <- ggplot(post2, aes(y = .value, x = as.factor(female_F), col = as.factor(female_F))) +
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(data = d_avg_NC, alpha = 0.3, 
              aes(y = avg_mass_res, x = as.factor(female_F), col = as.factor(female_F))) +
  scale_color_got(discrete = TRUE, option = "Tully", direction = 1, begin = 0) +
  theme_linedraw() +
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_y_continuous("Residual average larval mass") +
  scale_x_discrete("Female background", breaks = c(0,1), 
                   labels = c("No Care", "Full Care")) +
  labs(tag = "B")

grid.arrange(p1, p2, nrow = 1)

# 
# ##### residuals!
# Estimate   Est.Error          Q2.5         Q97.5
# [1,] -8.286201e-03 0.002067615 -1.233961e-02 -4.234910e-03
# [2,]  2.150848e-02 0.001996604  1.758349e-02  2.544773e-02
# [3,] -1.434866e-02 0.001996604 -1.827365e-02 -1.040941e-02
# [4,] -2.753947e-03 0.002565615 -7.804437e-03  2.284765e-03
# [5,]  7.738523e-03 0.001965650  3.877147e-03  1.159500e-02
# [6,]  1.766562e-02 0.001996604  1.374063e-02  2.160487e-02
# [7,] -1.134363e-02 0.001996054 -1.526569e-02 -7.390227e-03
# [8,]  4.148093e-04 0.002884264 -5.270885e-03  6.051451e-03
# [9,]  1.147785e-02 0.002097089  7.348020e-03  1.560150e-02
# [10,]  4.478217e-03 0.002001183  5.515908e-04  8.406029e-03
# [11,]  1.579316e-02 0.002043972  1.178120e-02  1.983211e-02
# [12,]  7.897848e-03 0.002097089  3.768020e-03  1.202150e-02
# [13,]  1.090019e-02 0.002114738  6.760852e-03  1.503540e-02
# [14,] -2.321767e-02 0.001946229 -2.702505e-02 -1.939773e-02
# [15,] -7.082152e-03 0.002097089 -1.121198e-02 -2.958503e-03
# [16,]  1.461785e-02 0.002097089  1.048802e-02  1.874150e-02
# [17,]  3.389512e-03 0.001932322 -4.095793e-04  7.181822e-03
# [18,] -4.253292e-02 0.001923798 -4.630654e-02 -3.875055e-02
# [19,]  1.477602e-02 0.001965650  1.091465e-02  1.863250e-02
# [20,]  3.001902e-02 0.002139517  2.581245e-02  3.422312e-02
# [21,]  3.773034e-03 0.001996054 -1.490215e-04  7.726440e-03
# [22,]  4.420129e-03 0.004266732 -3.935076e-03  1.276635e-02
# [23,]  6.300186e-03 0.002114738  2.160852e-03  1.043540e-02
# [24,]  1.883096e-02 0.002114738  1.469162e-02  2.296617e-02
# [25,]  1.559603e-02 0.002682698  1.031548e-02  2.087241e-02
# [26,]  3.762013e-02 0.004266732  2.926492e-02  4.596635e-02
# [27,]  4.312537e-02 0.002354272  3.849694e-02  4.774991e-02
# [28,]  2.109984e-02 0.002440767  1.629045e-02  2.589698e-02
# [29,]  1.109316e-02 0.002043972  7.081198e-03  1.513211e-02
# [30,]  9.425372e-03 0.002354272  4.796939e-03  1.404991e-02
# [31,]  5.619021e-03 0.002139517  1.412447e-03  9.823119e-03
# [32,]  2.658127e-02 0.002313888  2.202797e-02  3.112659e-02
# [33,] -5.160299e-03 0.001996054 -9.082355e-03 -1.206894e-03
# [34,] -9.280979e-03 0.002139517 -1.348755e-02 -5.076881e-03
# [35,]  1.085405e-02 0.002565615  5.803563e-03  1.589277e-02
# [36,]  1.175519e-02 0.002313888  7.201887e-03  1.630050e-02
# [37,]  2.528692e-02 0.003439755  1.854759e-02  3.201398e-02
# [38,]  3.971865e-03 0.001932322  1.727737e-04  7.764175e-03
# [39,] -1.902009e-02 0.001996604 -2.294508e-02 -1.508084e-02
# [40,] -2.787987e-02 0.004266732 -3.623508e-02 -1.953365e-02
# [41,] -6.305829e-04 0.002114738 -4.769917e-03  3.504631e-03
# [42,] -2.122565e-02 0.002085530 -2.534320e-02 -1.712397e-02
# [43,] -8.185841e-03 0.004343996 -1.674095e-02  3.121517e-04
# [44,] -4.342430e-03 0.002682698 -9.622977e-03  9.339461e-04
# [45,]  2.165508e-03 0.002014366 -1.783392e-03  6.117085e-03
# [46,] -1.250581e-02 0.001996604 -1.643080e-02 -8.566556e-03
# [47,]  1.018760e-02 0.001946229  6.380211e-03  1.400754e-02
# [48,]  2.742902e-02 0.002139517  2.322245e-02  3.163312e-02
# [49,]  1.540666e-02 0.002067615  1.135325e-02  1.945795e-02
# [50,]  5.206053e-03 0.002565615  1.555629e-04  1.024477e-02
# [51,]  9.103386e-03 0.001946229  5.296000e-03  1.292333e-02
# [52,] -3.285678e-02 0.002001183 -3.678341e-02 -2.892897e-02
# [53,] -1.013644e-03 0.004950232 -1.077128e-02  8.682170e-03
# [54,] -2.521814e-02 0.002799903 -3.073547e-02 -1.975535e-02
# [55,] -1.079178e-02 0.002001183 -1.471841e-02 -6.863971e-03
# [56,]  3.031785e-02 0.002097089  2.618802e-02  3.444150e-02
# [57,]  4.042370e-03 0.002067615 -1.104001e-05  8.093662e-03
# [58,] -1.191637e-02 0.001932322 -1.571546e-02 -8.124061e-03
# [59,] -5.237987e-02 0.004266732 -6.073508e-02 -4.403365e-02
# [60,]  3.275372e-03 0.002354272 -1.353061e-03  7.899908e-03
# [61,]  1.367067e-02 0.002440767  8.861282e-03  1.846781e-02
# [62,] -1.585977e-02 0.002155893 -2.009970e-02 -1.162814e-02
# [63,] -2.247987e-02 0.004266732 -3.083508e-02 -1.413365e-02
# [64,]  6.925372e-03 0.002354272  2.296939e-03  1.154991e-02
# [65,] -1.066793e-02 0.010565087 -3.140785e-02  1.006198e-02
# [66,] -1.219403e-02 0.001923798 -1.596765e-02 -8.411661e-03
# [67,]  1.065262e-03 0.002682698 -4.215285e-03  6.341638e-03
# [68,] -1.127723e-02 0.001996604 -1.520223e-02 -7.337984e-03
# [69,]  1.020571e-02 0.003560018  3.155642e-03  1.718226e-02
# [70,] -9.838275e-03 0.002114738 -1.397761e-02 -5.703062e-03
# [71,] -1.989152e-02 0.001996604 -2.381651e-02 -1.595227e-02
# [72,]  1.691217e-02 0.002014366  1.296327e-02  2.086375e-02
# [73,] -5.882916e-03 0.001923798 -9.656538e-03 -2.100550e-03
# [74,] -1.737090e-02 0.002313888 -2.192420e-02 -1.282559e-02
# [75,] -3.876783e-03 0.002001183 -7.803409e-03  5.102893e-05
# [76,] -7.145967e-03 0.002114738 -1.128530e-02 -3.010754e-03
# [77,] -4.449963e-02 0.002354272 -4.912806e-02 -3.987509e-02
# [78,] -3.944921e-04 0.002014366 -4.343392e-03  3.557085e-03
# [79,]  4.788632e-03 0.002085530  6.710817e-04  8.890316e-03
# [80,] -8.333089e-03 0.004950232 -1.809073e-02  1.362726e-03
# [81,]  2.694194e-03 0.001996604 -1.230797e-03  6.633444e-03
# [82,]  2.555973e-03 0.001923798 -1.217649e-03  6.338339e-03
# [83,] -6.343769e-04 0.001996604 -4.559368e-03  3.304873e-03
# [84,]  7.903860e-03 0.002146470  3.693298e-03  1.211918e-02
# [85,] -5.874628e-03 0.002354272 -1.050306e-02 -1.250092e-03
# [86,] -7.506783e-03 0.002001183 -1.143341e-02 -3.578971e-03
# [87,] -4.562948e-03 0.001996604 -8.487940e-03 -6.236984e-04
# [88,] -1.513781e-02 0.002146470 -1.934837e-02 -1.092249e-02
# [89,]  5.750840e-03 0.002313888  1.197539e-03  1.029615e-02
# [90,] -1.667463e-02 0.002354272 -2.130306e-02 -1.205009e-02
# [91,]  6.873034e-03 0.001996054  2.950979e-03  1.082644e-02
# [92,]  5.885922e-03 0.003066043 -1.550277e-04  1.189404e-02
# [93,]  2.038217e-03 0.002001183 -1.888409e-03  5.966029e-03
# [94,] -1.818620e-02 0.002067615 -2.223961e-02 -1.413491e-02
# [95,]  6.511302e-03 0.004343996 -2.043810e-03  1.500929e-02
# [96,]  8.980571e-03 0.005720203 -2.274702e-03  2.023270e-02
# [97,] -2.894215e-02 0.002097089 -3.307198e-02 -2.481850e-02
# [98,] -6.796136e-03 0.002155893 -1.103606e-02 -2.564508e-03
# [99,]  1.037991e-02 0.001996604  6.454918e-03  1.431916e-02
# [100,] -2.656783e-03 0.002001183 -6.583409e-03  1.271029e-03
# [101,]  9.749255e-03 0.003066043  3.708306e-03  1.575737e-02
# [102,] -1.635147e-02 0.002799903 -2.186881e-02 -1.088869e-02
# [103,]  1.181865e-03 0.002799903 -4.335473e-03  6.644647e-03
# [104,] -2.605806e-03 0.001996604 -6.530797e-03  1.333444e-03
# [105,] -3.844701e-03 0.002085530 -7.962252e-03  2.569823e-04
# [106,]  4.868445e-03 0.004343996 -3.686667e-03  1.336644e-02
# [107,] -1.961003e-02 0.002313888 -2.416333e-02 -1.506472e-02
# [108,] -1.572701e-02 0.003887620 -2.340636e-02 -8.097148e-03
# [109,] -8.277947e-03 0.002565615 -1.332844e-02 -3.239235e-03
# [110,]  7.041013e-03 0.002085530  2.923463e-03  1.114270e-02
# [111,]  2.076219e-02 0.002146470  1.655163e-02  2.497751e-02
# [112,]  2.302013e-02 0.004266732  1.466492e-02  3.136635e-02
# [113,] -3.600979e-03 0.002139517 -7.807553e-03  6.031195e-04
# [114,] -6.004328e-03 0.002440767 -1.081372e-02 -1.207189e-03
# [115,]  2.186816e-02 0.002043972  1.785620e-02  2.590711e-02
# [116,] -4.580979e-03 0.002139517 -8.787553e-03 -3.768805e-04
# [117,]  4.176023e-03 0.001965650  3.146466e-04  8.032496e-03
# [118,]  8.677848e-03 0.002097089  4.548020e-03  1.280150e-02
# [119,] -1.297411e-03 0.003066043 -7.338361e-03  4.710703e-03
# [120,]  1.624834e-02 0.001932322  1.244924e-02  2.004065e-02
# [121,]  2.592875e-02 0.002788855  2.045134e-02  3.140539e-02
# [122,]  7.059276e-03 0.002973379  1.203042e-03  1.286010e-02
# [123,]  3.123970e-02 0.001996054  2.731765e-02  3.519311e-02
# [124,]  1.281805e-02 0.002565615  7.767563e-03  1.785677e-02
# [125,]  1.387808e-02 0.002067615  9.824674e-03  1.792938e-02
# [126,]  7.039306e-03 0.001923798  3.265684e-03  1.082167e-02
# [127,] -9.377234e-03 0.001996604 -1.330223e-02 -5.437984e-03
# [128,] -1.415366e-02 0.002114738 -1.829299e-02 -1.001845e-02
# [129,] -7.948663e-03 0.001996604 -1.187365e-02 -4.009413e-03
# [130,]  1.177991e-02 0.001996604  7.854918e-03  1.571916e-02
# [131,] -7.079871e-03 0.004266732 -1.543508e-02  1.266345e-03
# [132,] -4.733495e-03 0.002440767 -9.542884e-03  6.364395e-05
# [133,] -1.417674e-02 0.002114738 -1.831607e-02 -1.004152e-02
# [134,]  1.022996e-05 0.002788855 -5.467179e-03  5.486867e-03
# [135,] -1.361308e-02 0.003439755 -2.035241e-02 -6.886021e-03
# [136,] -7.457856e-03 0.002313888 -1.201116e-02 -2.912545e-03
# [137,] -2.047398e-02 0.001965650 -2.433535e-02 -1.661750e-02
# [138,] -1.316178e-02 0.002001183 -1.708841e-02 -9.233971e-03
# [139,] -7.098724e-03 0.001932322 -1.089781e-02 -3.306414e-03
# [140,] -5.505947e-03 0.002565615 -1.055644e-02 -4.672347e-04
# [141,]  3.114853e-02 0.002799903  2.563119e-02  3.661131e-02
# [142,]  6.139701e-03 0.001996054  2.217645e-03  1.009311e-02
# [143,] -2.132341e-02 0.002155893 -2.556334e-02 -1.709178e-02
# [144,] -1.358225e-02 0.002192593 -1.790722e-02 -9.277896e-03
# [145,] -1.150763e-02 0.002067615 -1.556104e-02 -7.456338e-03



