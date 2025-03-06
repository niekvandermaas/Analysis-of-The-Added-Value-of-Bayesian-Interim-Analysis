#packages
library(readxl)
library(tidyverse)
library(lubridate)
library(survival)
library(survminer)
library(haven)
library(rjags)
library(tidybayes)
library(ggthemes)
library(bayesplot)

########################################################
# import data and exploratory analysis
# This can be skipped if you want to just run the model
# on some data which you already know
########################################################

#hovon 103
hovon_103_seli <- hovon_103 %>% filter(seli %in% 1)

######################################################################
#model specification Bayesian analysis selinexor complete remission
######################################################################

model_1 <- "model{

  y[1] ~ dbin(theta1, n)
  theta1 ~ dbeta(0.5, 0.5)
  y[2] ~ dbin(theta2, m)
  theta2 <- theta1 + diff
  diff ~ dnorm(0.0,1)
  p <- step(diff)

}"

interim_list <- list(
  data_50 = hovon_103_seli[1:51,] %>% mutate(arm = droplevels(arm)), 
  data_75 = hovon_103_seli[1:77,] %>% mutate(arm = droplevels(arm)),
  data_102 = hovon_103_seli %>% mutate(arm = droplevels(arm))
)

#get elapsed time
ia103_1 <- interim_list[[1]]
ia103_2 <- interim_list[[2]]
ia103_3 <- interim_list[[3]]

ss_103_1 <- min(ia103_1$dT0)
ss_103_2 <- min(ia103_2$dT0)
ss_103_3 <- min(ia103_3$dT0)

ia_103_1 <- max(ia103_1$ds1, na.rm = T)
ia_103_2 <- max(ia103_2$ds1, na.rm = T)
ia_103_3 <- max(ia103_3$ds1, na.rm = T)

#time elapsed since start of the study
(ss_103_1 %--% ia_103_1) %/% days(1) / 30.4368499
(ss_103_2 %--% ia_103_2) %/% days(1) / 30.4368499
(ss_103_3 %--% ia_103_3) %/% days(1) / 30.4368499

success <- lapply(interim_list, function(df) {
  with(df, tapply(df[["crind"]], arm, sum, na.rm = TRUE))
})

y1 <- sapply(success, "[[", 1)
y2 <- sapply(success, "[[", 2)
n <- sapply(lapply(interim_list, function(df) df %>% count(arm) %>% select(n)), "[[", "n")
data <- data.frame(y1 = y1, y2 = y2, n = n[1,], m = n[2,])

data_list <- lapply(1:nrow(data), function(i) list(y = c(data$y1[i], data$y2[i]), n = data$n[i], m = data$m[i]))
mcmc_samples <- lapply(data_list, function(dl) {
  model <- jags.model(textConnection(model_1), dl, n.chains = 3, n.adapt = 1000)
  update(model, 1000)
  coda.samples(model, c("diff", "theta1", "theta2", "p"), n.iter = 50000, inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 2022))
})

mcmc_samples[[1]] %>% spread_draws(diff) %>% median_hdi() #Interim 1
mcmc_samples[[2]] %>% spread_draws(diff) %>% median_hdi() #Interim 2
mcmc_samples[[3]] %>% spread_draws(diff) %>% median_hdi() #Interim 3

#traceplots
mcmc_trace(mcmc_samples[[1]], pars = c("diff"))
mcmc_trace(mcmc_samples[[2]], pars = c("diff"))
mcmc_trace(mcmc_samples[[3]], pars = c("diff"))

#Rhat
gelman.diag(mcmc_samples[[1]][,"diff"])
gelman.diag(mcmc_samples[[2]][,"diff"])
gelman.diag(mcmc_samples[[3]][,"diff"])


######################################################################
# posterior distribution selinexor complete remission
######################################################################
# results and treatment effect
theta_diff <- list()
prob <- list()
for (cr in 1:length(mcmc_samples)){
  theta_diff[[cr]] <- mcmc_samples[[cr]][[1]][,1]
  prob[[cr]] <- sum(theta_diff[[cr]]<0)/length(theta_diff[[cr]])
}

# Make density plots
prob <- unlist(prob)
post <- data.frame(theta_diff[[1]], theta_diff[[2]],
                   theta_diff[[3]])

colnames(post) <- c("Interim analysis 1", "Interim analysis 2", 
                    "Original interim analysis")
dat_text <- data.frame(
  label = paste0(sapply(prob[1:3], round, 2)*100, "%"),
  key   = c("Interim analysis 1", "Interim analysis 2", 
            "Original interim analysis"),
  mean = sapply(theta_diff, mean)
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
  xlab(NULL) + geom_text(
    data    = dat_text,
    mapping = aes(x = 0, y = 0.5, label = label),
    hjust = 1, vjust = -0.5) +
  facet_wrap(~key, scales = "free", ncol = 1)

#extract values for fill
to_fill <- data_frame(
  x = ggplot_build(p)$data[[1]]$x,
  y = ggplot_build(p)$data[[1]]$y,
  key = factor(ggplot_build(p)$data[[1]]$PANEL, levels = c(1,2,3), 
               labels = c("Interim analysis 1", "Interim analysis 2", 
                          "Original interim analysis")))

#plot
post %>% 
  gather() %>% 
  
  mutate(key = factor(key)) %>% 
  
  ggplot() +
  geom_density(aes(x = value), trim = TRUE,
               fill = "#01A2D9", size = 0.8, alpha = 1) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_area(data = to_fill[to_fill$x >= 0, ], 
            aes(x=x, y=y), fill = "#014d64") +
  stat_pointinterval(aes(x = value, y = 0),
                     point_interval = mode_hdi, .width = c(.95, .5)) +
  scale_x_continuous(breaks = c(-0.5,-0.25,0), labels = c("-50.0%", "-25.0%", "0.0%")) +
  scale_y_continuous(NULL, breaks = NULL) + 
  geom_text(
    data    = dat_text,
    mapping = aes(x = mean, y = 0, label = paste0(round(mean*100,1), "%")), vjust = -1) +
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
  xlab("Treatment difference") +
  scale_fill_economist()
