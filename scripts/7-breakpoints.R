### Long-term coastal lagoon fisheries landings
# Calculate breakpoints in environmental datasets

#Clear workspace and set wd
rm(list=ls(all=TRUE))
setwd("~/coastalrivercoupling")

#Install libraries 
#install.packages("strucchange")
library(strucchange)
library(changepoint)

library(tidyverse)
#remotes::install_github("lindeloev/mcp")
library(mcp)

# read in full data including river physico-chemical (from CHE, CAT), historical river flow (CHE),  fish catches (year, lagoon) and climatic data
data_full <- read.csv("outputs/river_fisheries_lagoon_spp_data.csv",row.names=1)
head(data_full)
str(data_full)

# This is for average variables across Year while holding lagoon
data_spp_env <- data_full %>% 
  filter(!Lagoon=="Platjola" & !Lagoon=="Clot - Baseta") %>% #filter out 
  droplevels() %>%
  group_by(Year, Lagoon) %>%
  dplyr::select(Year, Lagoon, flow_che, SRP_che, Lagoon, mean_biomass, mean_biomass_log, Chl.Total, TOC, NO3, NO2, NH4,
                qmedmes, Wind_direction, Wind_speed, Sea_level_pressure, Precipitation, Mean_temperature) %>%
  summarise(across(everything(), mean, na.rm=TRUE)) %>%
  as.data.frame()



## Read in environmental data
# CHE
data_che <- read.csv("data/data_che_full.csv", sep=";")
head(data_che)
str(data_che)

# transform to numeric
df1 <- data.frame(apply(apply(data_che[,7:ncol(data_che)], 2, gsub, patt=",", replace="."), 2, as.numeric))
data_che <- cbind(data_che[,1:6], df1) 
data_che_Tortosa <- data_che %>% filter(Sampling_Point==2) #filter cases Tortosa station

plot.ts(data_che_Tortosa$Caudal.en.superficie..m3.s.)

#Calculate mean annual flow Data CHE (1980-2004)
year_mean_che <- data_che_Tortosa %>%
  mutate(Year_f=factor(Year)) %>%
  rename(flow=Caudal.en.superficie..m3.s.) %>%
  #filter(!is.na(CaudalÂ.enÂ.superficieÂ..m3.s)) %>% #remove NANs
  group_by(Year_f) %>% 
  summarise_at(c("flow", "Fosfatos..mg.L.PO4."), mean, na.rm = TRUE) %>%
  filter(!is.na(Fosfatos..mg.L.PO4.)) #remove NANs

plot(cpt.meanvar(year_mean_che$flow, method="BinSeg",pen.value=0.01))
plot(cpt.meanvar(year_mean_che$Fósforo.Total..mg.L.P., method="BinSeg",pen.value=0.01))
plot(cpt.meanvar(year_mean_che$Fosfatos..mg.L.PO4., method="BinSeg",pen.value=0.01))

cpt.meanvar(year_mean_che$Fosfatos..mg.L.PO4., method="BinSeg",pen.value=0.01)
# where is the breakpoint?
year_mean_che[c(15),] #1995 is breakpoints for PO4
PO4_breakpoinnt <- ifelse(data_spp_env$Year==1995, 1, 0)

# CAT dataset
data_cat <- read.csv("data/data_cat_full.csv", sep=";")
head(data_cat)
str(data_cat)
data_cat <- data_cat[-1,] #remove first row which contains variable units

# transform to numeric
df2 <- data.frame(apply(apply(data_cat[,3:ncol(data_cat)], 2, gsub, patt=",", replace="."), 2, as.numeric))
data_cat <- cbind(data_cat[,1:2], df2)
plot.ts(data_cat$Chl.Total)

year_mean_chla <- data_cat %>%
  mutate(Year_f=factor(Year)) %>%
  group_by(Year_f) %>% 
  summarise_at(c("Chl.Total", "TOC"), mean, na.rm = T) %>%
  #mutate(mean_annual_chla=round(mean_annual_chla,2)) %>%
  filter(!is.na(Chl.Total)) #remove NANs

# calculates the optimal positioning and number of significant changepoints
plot(cpt.meanvar(year_mean_chla$Chl.Total, method="BinSeg",pen.value=0.01))
cpt.meanvar(year_mean_chla$Chl.Total, method="BinSeg",pen.value=0.01)

# where is the breakpoint?
year_mean_chla[c(6,12),] #1995 and 2001 are breakpoints for chla
# create dummy variable for Chla; see the input dataframe: either form 4-HGAM.R or 8-GAM-trends.R
chla_breakpoinnt <- ifelse(data_spp_env$Year==1995 | data_spp_env$Year==2001, 1, 0)

plot(cpt.meanvar(year_mean_chla$TOC, method="BinSeg",pen.value=0.01))
cpt.meanvar(year_mean_chla$TOC, method="BinSeg",pen.value=0.01)

# where is the breakpoint?
year_mean_chla[c(7,10),] #1996 and 1999 are breakpoints for TOC
# create dummy variable for TOC; see the input dataframe: either form 4-HGAM.R or 8-GAM-trends.R
TOC_breakpoinnt <- ifelse(data_spp_env$Year==1996 | data_spp_env$Year==1999, 1, 0)


# CEDEX historical flows (Datos mensuales de estaciones de aforo en río)
hist_flow <- read.csv("data/mensual_a_EBRO.csv", sep=";") %>%
  filter(indroea==9027) #filter cases Tortosa station

# this is to create year and month columns
test <- data.frame(str_split_fixed(hist_flow$anomes, "", n = 6))
hist_flow$year <- str_c(test$X1,'',test$X2,'',test$X3,'',test$X4)
hist_flow$month <- str_c(test$X5,'',test$X6)

head(hist_flow)
str(hist_flow)

plot.ts(hist_flow$qmedmes)

#Calculate mean annual flow Data CHE CEDEX (1912-2018)
range(hist_flow$year) #1912-2018
year_mean_flow <- hist_flow %>%
  mutate(Year_f=factor(year)) %>%
  rename(flow=qmedmes) %>%
  group_by(Year_f) %>% 
  summarise(mean_annual_flow=mean(flow, na.rm = T)) %>%
  mutate(mean_annual_flow=round(mean_annual_flow,2)) %>%
  filter(!is.na(mean_annual_flow)) #remove NANs
  
plot.ts(year_mean_flow$mean_annual_flow)

# calculates the optimal positioning and number of significant changepoints
plot(cpt.meanvar(year_mean_flow$mean_annual_flow, method="BinSeg",pen.value=0.01))
cpt.meanvar(year_mean_flow$mean_annual_flow, method="BinSeg",pen.value=0.01)
# where is the breakpoint?
year_mean_flow[52,] #1979
# create dummy variable for TOC; see the input dataframe: either form 4-HGAM.R or 8-GAM-trends.R
flow_breakpoinnt <- ifelse(data_spp_env$Year==1979, 1, 0)

# Merge breakpoint vectors here and with data spp_env
breakpoints_data <- data.frame(data_spp_env$Year, chla_breakpoinnt, flow_breakpoinnt, TOC_breakpoinnt, PO4_breakpoinnt)
colnames(breakpoints_data)[1] <- "Year"

# save dataset
write.csv(breakpoints_data, "outputs/breakpoints_river.csv")

data_spp_env <- cbind(data_spp_env, breakpoints_data)
  
  
# A plethora of change point analysis: https://lindeloev.github.io/mcp/articles/packages.html

library(EnvCpt)
fit_envcpt <- envcpt(year_mean_flow$mean_annual_flow)  # Fit all models at once
fit_envcpt$summary  # Show log-likelihoods
plot(fit_envcpt)

fit_envcpt$meancpt@cpts
