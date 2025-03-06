set.seed(123)
library(brms)
library(tidybayes)
library(ggthemes)
library(tidyverse)

##############################################
#Bayesian interim analysis progression free survival
##############################################
#prior 
prior_1 <- c(set_prior("normal(0,10)", class = "b"))

bfit_1_87 <- brm(pfs_interim | cens(1-pfs_interim_i) ~ arm, data = interim_analysis_1, family = cox(), iter = 50000, cores = 32, backend = "cmdstanr", prior = prior_1)
bfit_2_87 <- update(bfit_1_87, newdata = interim_analysis_2, cores = 4)
bfit_3_87 <- update(bfit_1_87, newdata = interim_analysis_3, cores = 4)

bfit_1_87 %>% spread_draws(b_arm) %>% mutate(b_arm = exp(b_arm)) %>% median_hdi()
bfit_2_87 %>% spread_draws(b_arm) %>% mutate(b_arm = exp(b_arm)) %>% median_hdi()
bfit_3_87 %>% spread_draws(b_arm) %>% mutate(b_arm = exp(b_arm)) %>% median_hdi()

mcmc_trace(bfit_1_87, pars = c("b_Intercept", "b_arm"))
mcmc_trace(bfit_2_87, pars = c("b_Intercept", "b_arm"))
mcmc_trace(bfit_3_87, pars = c("b_Intercept", "b_arm"))

#get posterior draws to compile posterior distribution plots
posterior_draws_1 <- bfit_1_87 %>% spread_draws(b_arm) %>% mutate(b_arm = exp(b_arm)) %>% select(b_arm) %>% rename(b_arm_1 = b_arm)
posterior_draws_2 <- bfit_1_87 %>% spread_draws(b_arm) %>% mutate(b_arm = exp(b_arm)) %>% select(b_arm) %>% rename(b_arm_2 = b_arm)
posterior_draws_3 <- bfit_1_87 %>% spread_draws(b_arm) %>% mutate(b_arm = exp(b_arm)) %>% select(b_arm) %>% rename(b_arm_3 = b_arm)

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
  
  prob_1[[i]] <- calc_auc_post(hr_1[[i]], 0.714)
  prob_2[[i]] <- calc_auc_post(hr_1[[i]], 1)
}

#make surival curves
ggsurvplot(
  survfit(Surv(pfs_interim, pfs_interim_i) ~ arm, interim_analysis_1),
  size = 1,                 # change line size
  risk.table = T,
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)


# Make some nice density plots
prob_1 <- unlist(prob_1)
prob_2 <- unlist(prob_2)
post <- data.frame(hr_1[[1]], hr_1[[2]], 
                   hr_1[[3]])

colnames(post) <- c("Interim analysis 1", "Interim analysis 2", "Interim analysis 3")

dat_text_1 <- data.frame(
  label = paste("HR < 0.714 =", prob_1),
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
  geom_vline(aes(xintercept = 0.714), linetype = "dashed") +
  scale_y_continuous(NULL, breaks = NULL) +
  scale_x_continuous(limits = c(0, 2), breaks = seq(0, 2, by = 0.5)) + 
  xlab(NULL) +
  #area 1
  geom_area(data = to_fill[to_fill$x <= 1, ], 
            aes(x=x, y=y), fill = "#01A2D9",
            alpha = 0.75) +
  #area 2
  geom_area(data = to_fill[to_fill$x <= 0.714, ], 
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
