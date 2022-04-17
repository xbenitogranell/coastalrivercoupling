## Random Forest to reduce the number of predictors

#Clear workspace
rm(list=ls(all=TRUE))

library(mgcv)
library(tidyverse)
library(ggplot2)
#library(gratia) # for visualizing GAM model outputs
library(party) #to perform Random Forest

# read in full data including river physico-chemical (from CHE, CAT), river flow (CHE), historical river flow (CHE), and mean fish catches (averaged across lagoons)
data_full <- read.csv("outputs/river_fisheries_lagoon_spp_data.csv",row.names=1)
str(data_full)

# read in breakpoints river environment
breakpoints_river <- read.csv("outputs/breakpoints_river.csv", row.names = 1)

# check how many levels per species and lagoon
levels(factor(data_full$Species))
levels(factor(data_full$Lagoon))

# This is for average monthly observations across Years to match with biomass catches
data_spp_env <- data_full %>% 
  filter(!Lagoon=="Platjola" & !Lagoon=="Clot - Baseta") %>% #filter out 
  droplevels() %>%
  group_by(Year) %>%
  dplyr::select(-Month, -Year) %>%
  # dplyr::select(Year, flow_che, SRP_che, mean_biomass, mean_biomass_log, Chl.Total, TOC, NO3, NO2, NH4,
  #               qmedmes, Wind_direction, Wind_speed, Sea_level_pressure, Precipitation, Mean_temperature, QMax, QMean) %>%
  summarise(across(everything(), mean, na.rm=TRUE)) %>%
  mutate(qmedmes_lag1=lag(qmedmes)) %>% #create a lagged variable for historical flow
  left_join(breakpoints_river, by="Year") %>% #join with river environmental dataset of breakpoints
  as.data.frame()

str(data_spp_env)

# Plot of possible transformations
pred.vars <- colnames(data_spp_env)[-1]

par(mfrow=c(3,2))
for (i in pred.vars) {
  x <- data_spp_env[ ,i]
  x = as.numeric(unlist(x))
  hist((x)) 
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log10(x+0.25))
  plot(log10(x+0.25))
}

# Make transformations
data_spp_env$alkalinity <- log10(data_spp_env$alkalinity + 0.25)
data_spp_env$ammonia_che <- log10(data_spp_env$ammonia_che + 0.25)
data_spp_env$flow_che <- sqrt(data_spp_env$flow_che)
data_spp_env$DBO_che <- log10(data_spp_env$flow_che + 0.25)
data_spp_env$SRP_che <- log10(data_spp_env$SRP_che + 0.25)
data_spp_env$total_phosphorous_che <- log10(data_spp_env$total_phosphorous_che + 0.25)
data_spp_env$SuspendedSolids_che <- log10(data_spp_env$SuspendedSolids_che + 0.25)
data_spp_env$total_nitrogen_che <- sqrt(data_spp_env$total_nitrogen_che)
data_spp_env$WaterT_che <- log10(data_spp_env$WaterT_che + 0.25)


# Run Random Forest
#use.dat <- na.omit(data_spp_env[,c(response.var,pred.vars)]) #get rid off NAs
use.dat <- data_spp_env %>% select(-Year, -mean_biomass)

set.seed(42)
catches.rf <- cforest(mean_biomass_log ~ ., data = use.dat, control = cforest_unbiased(mtry = 5, ntree = 200))
barplot(varimp(catches.rf))
