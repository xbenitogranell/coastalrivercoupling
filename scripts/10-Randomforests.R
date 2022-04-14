## Random Forest to reduce the number of predictors

#Clear workspace
rm(list=ls(all=TRUE))

library(mgcv)
library(tidyverse)
library(ggplot2)
#library(gratia) # for visualizing GAM model outputs

# read in full data including river physico-chemical (from CHE, CAT), river flow (CHE), historical river flow (CHE), and mean fish catches (averaged across lagoons)
data_full <- read.csv("outputs/river_fisheries_lagoon_spp_data.csv",row.names=1)
str(data_full)

# check how many levels per species and lagoon
levels(factor(data_full$Species))
levels(factor(data_full$Lagoon))

plot.ts(data_full$Year, data_full$qmedmes, type="l")
plot.ts(data_full$Year, data_full$mean_biomass, type="l")

# This is for average monthly observations across Years
data_spp_env <- data_full %>% 
  filter(!Lagoon=="Platjola" & !Lagoon=="Clot - Baseta") %>% #filter out 
  droplevels() %>%
  group_by(Year) %>%
  summarise_at(c("flow_che", "qmedmes"), mean, na.rm = TRUE)
  # dplyr::select(Year, flow_che, SRP_che, mean_biomass, mean_biomass_log, Chl.Total, TOC, NO3, NO2, NH4,
  #               qmedmes, Wind_direction, Wind_speed, Sea_level_pressure, Precipitation, Mean_temperature) %>%
  summarise(across(everything(), mean, na.rm=TRUE)) %>%
  as.data.frame()

