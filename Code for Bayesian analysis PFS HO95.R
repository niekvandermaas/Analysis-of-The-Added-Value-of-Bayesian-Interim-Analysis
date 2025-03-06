# Bayesian interim analyses based on patient numbers HOVON95
# July 17 2023

# Total patients: n = 1197
# INTERIM 1 - 50% patients: n = 599
# INTERIM 2 - 75% patients: n = 898
# INTERIM 3 - 100% patients: n = 1197



set.seed(2024)
library("rjags")
library("tidyverse")
library("tidybayes")
library("survminer")
library("survival")
library('bayesplot')
library("ggthemes")
library(brms)

##############################################
#Bayesian cox model
##############################################
#Scenario: weakly informative priors for betas 
prior_1 <- c(set_prior("normal(0,10)", class = "b"))

#use a flat prior (non-informative)
bfit_1_95 <- brm(pfs_interim | cens(1-pfs_interim_i) ~ armR1c2, data = patient_interim_R1_599, 
                 family = cox(), iter = 50000, cores = 4, 
                 backend = "cmdstanr", prior = prior_1)
bfit_2_95 <- update(bfit_1_95, newdata = patient_interim_R1_898, cores = 4)
bfit_3_95 <- update(bfit_1_95, newdata = patient_interim_R1_1197, cores = 4)

bfit_1_95 %>% spread_draws(b_armR1c2) %>% mutate(b_armR1c2 = exp(b_armR1c2)) %>% median_hdi()
bfit_2_95 %>% spread_draws(b_armR1c2) %>% mutate(b_armR1c2 = exp(b_armR1c2)) %>% median_hdi()
bfit_3_95 %>% spread_draws(b_armR1c2) %>% mutate(b_armR1c2 = exp(b_armR1c2)) %>% median_hdi()

#get posterior draws to compile posterior distribution plots
posterior_draws_1 <- bfit_1_95 %>% spread_draws(b_armR1c2) %>% mutate(b_armR1c2 = exp(b_armR1c2)) %>% select(b_armR1c2) %>% rename(b_armR1c2_1 = b_armR1c2)
posterior_draws_2 <- bfit_2_95 %>% spread_draws(b_armR1c2) %>% mutate(b_armR1c2 = exp(b_armR1c2)) %>% select(b_armR1c2) %>% rename(b_armR1c2_2 = b_armR1c2)
posterior_draws_3 <- bfit_3_95 %>% spread_draws(b_armR1c2) %>% mutate(b_armR1c2 = exp(b_armR1c2)) %>% select(b_armR1c2) %>% rename(b_armR1c2_3 = b_armR1c2)


#plot posterior HR
#function for auc posterior
calc_auc_post <- function(data, threshold) {
  prob <- sum(data<threshold)/length(data)
  prb <- paste0(sapply(prob, round, 3)*100,"%")
  return(prb)
}

# results in one tibble for analysis
hr_1 <- list()
prob_1 <- list()
prob_2 <- list()
results <- cbind(posterior_draws_1, posterior_draws_2, posterior_draws_3)
for (i in 1:ncol(results)){
  hr_1[[i]] <- results[,i]
  
  prob_1[[i]] <- calc_auc_post(hr_1[[i]], 0.78)
  prob_2[[i]] <- calc_auc_post(hr_1[[i]], 1)
}


# Make some nice density plots
prob_1 <- unlist(prob_1)
prob_2 <- unlist(prob_2)
post <- data.frame(hr_1[[1]], hr_1[[2]], 
                   hr_1[[3]])

colnames(post) <- c("Interim analysis 1", "Interim analysis 2", "Interim analysis 3")

dat_text_1 <- data.frame(
  label = paste("HR < 0.78 =", prob_1),
  key   = c("Interim analysis 1", "Interim analysis 2", "Interim analysis 3")
)

dat_text_2 <- data.frame(
  label = paste("HR < 1 =", prob_2),
  key   = c("Interim analysis 1", "Interim analysis 2", "Interim analysis 3")
)

# plot posteriors with HDI
p <- post %>% 
  gather() %>% 
  
  mutate(key = factor(key)) %>% 
  
  ggplot() +
  geom_density(aes(x = value),
               fill = "grey92") +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  stat_pointinterval(aes(x = value, y = 0),
                     point_interval = mode_hdi, .width = c(.95, .5)) +
  scale_y_continuous(NULL, breaks = NULL) +
  xlab(NULL) +
  facet_wrap(~key, scales = "free", ncol = 2)

#extract values for fill
to_fill <- data_frame(
  x = ggplot_build(p)$data[[1]]$x,
  y = ggplot_build(p)$data[[1]]$y,
  key = factor(ggplot_build(p)$data[[1]]$PANEL, levels = c(1,2,3), 
               labels = c("Interim analysis 1", "Interim analysis 2", "Interim analysis 3"), ordered = T))

post %>% 
  gather() %>% 
  mutate(key = factor(key)) %>% 
  ggplot() +
  geom_density(aes(x = value), trim = TRUE,
               fill = "aliceblue", size = 0.8, alpha = 0.75) +
  geom_vline(aes(xintercept = 1), linetype = "dashed") +
  geom_vline(aes(xintercept = 0.78), linetype = "dashed") +
  scale_y_continuous(NULL, breaks = NULL) +
  scale_x_continuous(limits = c(0, 2), breaks = seq(0, 2, by = 0.5)) + 
  xlab(NULL) +
  #area 1
  geom_area(data = to_fill[to_fill$x <= 1, ], 
            aes(x=x, y=y), fill = "#01A2D9",
            alpha = 0.75) +
  #area 2
  geom_area(data = to_fill[to_fill$x <= 0.78, ], 
            aes(x=x, y=y), fill = "#014d64",
            alpha = 1) +
  stat_pointinterval(aes(x = value, y = 0),
                     point_interval = mode_hdi, .width = c(.95, .5)) +
  geom_text(
    data    = dat_text_1,
    mapping = aes(x = 1.5, y = 1, label = label)) +
  geom_text(
    data    = dat_text_2,
    mapping = aes(x = 1.5, y = 1.5, label = label)) +
  facet_wrap(~key, scales = "fixed", ncol = 1, strip.position = "left") + 
  theme_fivethirtyeight() + 
  theme(
    strip.text.y.left = element_text(angle = 90),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = "white"),
    axis.title.x =  element_text()
  ) +
  xlab("Hazard ratio") +
  scale_fill_economist()
