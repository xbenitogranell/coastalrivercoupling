## This is to calculate first derivative from the HGAM model posteriors as a tool to detect regime shifts

#Clear workspace
rm(list=ls(all=TRUE))

# load functions used
library(tidyverse)
library(ggplot2)

# Read in HGAM S & I fitted values
fishes_plot_data <- read.csv("outputs/HGAM_fishes_data.csv", row.names=1)
head(fishes_plot_data)

##Derivatives and posterior distribution simulation
set.seed(10) #set a seed so this is repeatable
n_sims = 250

n_length = 44

years <- seq(min(fishes_plot_data$Year2),
             max(fishes_plot_data$Year2),
             length.out = n_length)

model <- model_gam_S
#model <- model_gam_I

pred <- fishes_modI_fit
#pred <- fishes_modI_fit

# Generate multivariate normal simulations of model coefficients
random_coefs <- t(rmvn(n_sims, mu = coef(model),V = vcov(model)))

confint_sims <- crossing(spp=unique(fishes_plot_data$spp),
                         Year2 = seq(min(fishes_plot_data$Year2),
                                      max(fishes_plot_data$Year2),
                                      length.out = n_length),
                         log_total_counts=0)

map_pred_sims <- predict(model,
                         confint_sims,
                         type = "lpmatrix") %*% random_coefs %>%
  as_data_frame() %>%
  bind_cols(confint_sims)%>%
  gather(key = simulation, value = pred, -Year2,-spp)


#specifying the step size for numerical derivative calculations
delta = 0.01

#calculating the predicted value for the current year plus delta
step_ahead_fits <- confint_sims %>%
  mutate(Year2 = Year2+delta)%>%
  predict(model, 
          ., type = "lpmatrix") %*% random_coefs 


#calculating the predicted value for the current year minus delta
step_behind_fits <- confint_sims %>%
  mutate(Year2 = Year2-delta)%>%
  predict(model,
          ., type = "lpmatrix") %*% random_coefs 



#Function for calculating first derivatives of time series given a before,
#after, and delta step size
calc_1st_deriv = function(fit_before, fit_after,delta) {
  (fit_after-fit_before)/(2*delta)
}

#using the predicted values for year plus and minus delta to calculate
#derivatives for each species for each simulation
derivs <- calc_1st_deriv(step_behind_fits,step_ahead_fits,delta = delta)%>%
  as_data_frame()%>%
  bind_cols(confint_sims)%>%
  gather(key = simulation,value = deriv, -spp,-Year2)

#Creating summaries of derivatives for each simulation for each year
deriv_summaries <- derivs %>%
  group_by(Year2,simulation)%>%
  summarize(deriv_mean = mean(deriv),
            deriv_sd = sd(deriv))%>%
  group_by(Year2)%>% #turning derivative summaries into 95% confidence intervals
  select(-simulation)%>%
  summarize_all(.funs = list(lower = ~quantile(.,probs = 0.025),
                             upper = ~quantile(.,probs = 0.975),
                             med   = ~quantile(.,probs = 0.5)))

#Plotting mean rate of change plus the 95% CI
mean_plot <- deriv_summaries %>%
  ggplot(aes(x = Year2, 
             y = deriv_mean_med, 
             ymin = deriv_mean_lower,
             ymax = deriv_mean_upper))+
  geom_ribbon(fill="grey")+
  geom_line()+
  geom_hline(yintercept = 0, linetype=2) +
  scale_y_continuous("")+
  xlab("Years") +
  theme_bw()
mean_plot 

#Plotting standard deviation of rate of change plus the 95% CI
sd_plot <- deriv_summaries %>%
  ggplot(aes(x = Year2, 
             y = deriv_sd_med, 
             ymin=deriv_sd_lower,
             ymax=deriv_sd_upper))+
  geom_ribbon(fill="grey")+
  geom_line()+
  geom_hline(yintercept = 0, linetype=2) +
  scale_y_continuous("")+
  xlab("Years") +
  theme_bw()
sd_plot


## save derivative summaries for later use
#write.csv(deriv_summaries, "outputs/.csv", row.names = FALSE)