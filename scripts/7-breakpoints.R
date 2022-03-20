### Long-term coastal lagoon fisheries landings
# Calculate breakpoints in environmental datasets

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


## Read in environmental data
# CHE
data_che <- read.csv("data/data_che_full.csv", sep=";")
head(data_che)
str(data_che)

# transform to numeric
df1 <- data.frame(apply(apply(data_che[,7:ncol(data_che)], 2, gsub, patt=",", replace="."), 2, as.numeric))
data_che <- cbind(data_che[,1:6], df1) 
data_che_Tortosa <- data_che %>% filter(Sampling_Point==2) #filter cases Tortosa station

plot.ts(data_che_Tortosa$CaudalÂ.enÂ.superficieÂ..m3.s.)

#Calculate mean annual flow Data CHE (1980-2004)
year_mean_flow <- data_che_Tortosa %>%
  mutate(Year_f=factor(Year)) %>%
  rename(flow=CaudalÂ.enÂ.superficieÂ..m3.s.) %>%
  #filter(!is.na(CaudalÂ.enÂ.superficieÂ..m3.s)) %>% #remove NANs
  group_by(Year_f) %>% 
  summarise(mean_annual_flow=mean(flow, na.rm = T)) %>%
  mutate(mean_annual_flow=round(mean_annual_flow,2)) 

plot(cpt.meanvar(year_mean_flow$mean_annual_flow, method="BinSeg",pen.value=0.01))


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
  summarise(mean_annual_chla=mean(Chl.Total, na.rm = T)) %>%
  mutate(mean_annual_chla=round(mean_annual_chla,2)) %>%
  filter(!is.na(mean_annual_chla)) #remove NANs

plot(cpt.meanvar(year_mean_chla$mean_annual_chla, method="BinSeg",pen.value=0.01))
cpt.meanvar(year_mean_chla$mean_annual_chla, method="BinSeg",pen.value=0.01)

year_mean_chla[c(6,12),]

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
range(hist_flow$Year) #1912-2018
year_mean_flow <- hist_flow %>%
  mutate(Year_f=factor(Year)) %>%
  rename(flow=qmedmes) %>%
  group_by(Year_f) %>% 
  summarise(mean_annual_flow=mean(flow, na.rm = T)) %>%
  mutate(mean_annual_flow=round(mean_annual_flow,2)) %>%
  filter(!is.na(mean_annual_flow)) #remove NANs
  
plot.ts(year_mean_flow$mean_annual_flow)
plot(cpt.meanvar(year_mean_flow$mean_annual_flow, method="BinSeg",pen.value=0.01))
cpt.meanvar(year_mean_flow$mean_annual_flow, method="BinSeg",pen.value=0.01)

year_mean_flow[52,]

## A plethora of change point analysis: https://lindeloev.github.io/mcp/articles/packages.html

library(EnvCpt)
fit_envcpt <- envcpt(year_mean_flow$mean_annual_flow)  # Fit all models at once
fit_envcpt$summary  # Show log-likelihoods
plot(fit_envcpt)

fit_envcpt$meancpt@cpts
