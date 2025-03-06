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


######################################################################
# Bayesian analysis AE 3-5 HO95
######################################################################
interim_list <- list(
  data_1 = patient_interim_R1_599 %>%
    left_join(ho95_ae35nohem, by = "pid") %>%
    mutate(ae_35 = ifelse(max(dR1) + dmonths(2) >= dstartae_35, 1, 0),
           trt = factor(armR1c2, labels = c("VMP", "HSCT"))), 
  data_2 = patient_interim_R1_898 %>%
    left_join(ho95_ae35nohem, by = "pid") %>%
    mutate(ae_35 = ifelse(max(dR1) + dmonths(2) >= dstartae_35, 1, 0),
           trt = factor(armR1c2, labels = c("VMP", "HSCT"))),
  data_3 = patient_interim_R1_1197 %>%
    left_join(ho95_ae35nohem, by = "pid") %>%
    mutate(ae_35 = ifelse(max(dR1) + dmonths(2) >= dstartae_35, 1, 0),
           trt = factor(armR1c2, labels = c("VMP", "HSCT")))
)

# Model
model_2 <- "model{ 
  y[1] ~ dbin(theta1, n)
  theta1 ~ dbeta(0.5, 0.5)
  y[2] ~ dbin(theta2, m)
  theta2 <- theta1 + diff
  diff ~ dnorm(0.0,1) #precision: 1/var
  p <- step(diff)
}"


success <- lapply(interim_list, function(df) {
  with(df, tapply(df[["ae_35"]], trt, sum, na.rm = TRUE))
})

y1 <- sapply(success, "[[", 1)
y2 <- sapply(success, "[[", 2)
n <- sapply(lapply(interim_list, function(df) df %>% count(trt) %>% select(n)), "[[", "n")
data <- data.frame(y1 = y1, y2 = y2, n = n[1,], m = n[2,])

data_list <- lapply(1:nrow(data), function(i) list(y = c(data$y1[i], data$y2[i]), n = data$n[i], m = data$m[i]))
mcmc_samples <- lapply(data_list, function(dl) {
  model <- jags.model(textConnection(model_2), dl, n.chains = 3, n.adapt = 1000)
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
# posterior distribution AE 3-5
######################################################################
# results and treatment effect
theta_diff <- list()
prob <- list()
for (cr in 1:length(mcmc_samples)){
  theta_diff[[cr]] <- mcmc_samples[[cr]][[1]][,1]
  prob[[cr]] <- sum(theta_diff[[cr]]>0)/length(theta_diff[[cr]])
}

# Make density plots
prob <- unlist(prob)
post <- data.frame(theta_diff[[1]], theta_diff[[2]],
                   theta_diff[[3]])

colnames(post) <- c("Interim analysis 1", "Interim analysis 2", "Interim analysis 3")
dat_text <- data.frame(
  label = paste0(sapply(prob[1:3], round, 2)*100, "%"),
  key   = c("Interim analysis 1", "Interim analysis 2", "Interim analysis 3"),
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
               labels = c("Interim analysis 1", "Interim analysis 2", "Interim analysis 3")))

#plot
post %>% 
  gather() %>% 
  
  mutate(key = factor(key)) %>% 
  
  ggplot() +
  geom_density(aes(x = value), trim = TRUE,
               fill = "#014d64", size = 0.8, alpha = 1) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_area(data = to_fill[to_fill$x >= 0, ], 
            aes(x=x, y=y), fill = "#01A2D9",
            alpha = 0.75) +
  stat_pointinterval(aes(x = value, y = 0),
                     point_interval = mode_hdi, .width = c(.95, .5)) +
  scale_x_continuous(breaks = c(-0.1, 0, 0.1,0.2), labels = c("-25.0%", "0.0%", "10.0%", "20.0%")) +
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

