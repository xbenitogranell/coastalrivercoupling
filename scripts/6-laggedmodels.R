## EXPLORATORY ANALYSIS OF THE EFFECT OF LAGGED RIVER DYNAMICS ON COASTAL FISHERIES 

#Clear workspace
rm(list=ls(all=TRUE))

#load lag functions
source("scripts/functions/functions_lags.R")

#Load functions used
library(tidyverse)
library(ggplot2)
library(cowplot)
library(nlme)

# Read in CHE data
data_che <- read.csv("data/data_che_full.csv", sep=";")
head(data_che)
str(data_che)

# transform to numeric
df1 <- data.frame(apply(apply(data_che[,7:ncol(data_che)], 2, gsub, patt=",", replace="."), 2, as.numeric))
data_che <- cbind(data_che[,1:6], df1) 
data_che_Tortosa <- data_che %>% filter(Sampling_Point==2) #filter cases Tortosa station

#Calculate mean annual flow Data CHE (1980-2004) to match with fish catches
year_mean_flow_che <- data_che_Tortosa %>%
  mutate(Year_f=factor(Year)) %>%
  rename(flow=CaudalÂ.enÂ.superficieÂ..m3.s.) %>%
  #filter(!is.na(CaudalÂ.enÂ.superficieÂ..m3.s)) %>% #remove NANs
  group_by(Year_f) %>% 
  summarise(mean_annual_flow=mean(flow, na.rm = T)) %>%
  mutate(mean_annual_flow=round(mean_annual_flow,2)) 

## CAT
data_cat <- read.csv("data/data_cat_full.csv", sep=";") %>%
  rename(Month=ï..Month)
head(data_cat)
str(data_cat)
data_cat <- data_cat[-1,] #remove first row which contains variable units

# transform to numeric
df2 <- data.frame(apply(apply(data_cat[,3:ncol(data_cat)], 2, gsub, patt=",", replace="."), 2, as.numeric))
data_cat <- cbind(data_cat[,1:2], df2)
plot.ts(data_cat$Chl.Total)

#Calculate mean annual chla data to match with fish catches
year_mean_chla <- data_cat %>%
  mutate(Year_f=factor(Year)) %>%
  group_by(Year_f) %>% 
  summarise(mean_annual_chla=mean(Chl.Total, na.rm = T)) %>%
  mutate(mean_annual_chla=round(mean_annual_chla,2)) %>%
  filter(!is.na(mean_annual_chla)) #remove NANs

## CEDEX historical flows (Datos mensuales de estaciones de aforo en río)
hist_flow <- read.csv("data/mensual_a_EBRO.csv", sep=";") %>%
  filter(indroea==9027) #filter cases Tortosa station

# this is to create year and month columns
test <- data.frame(str_split_fixed(hist_flow$anomes, "", n = 6))
hist_flow$Year <- str_c(test$X1,'',test$X2,'',test$X3,'',test$X4)
hist_flow$Month <- str_c(test$X5,'',test$X6)
hist_flow <- hist_flow %>% mutate(Year=as.numeric(Year)) %>%
  mutate(Month=as.numeric(Month))

head(hist_flow)
str(hist_flow)

#Read in catches lagoons clean data
catches_clean <- read.csv("data/catches_clean.csv")[-1]
str(catches_clean)

#Species to be included in the models (>30% average biomass)
include <- c("Anguilla anguilla", "Cyprinus carpio", "Mullet ind.", "Sparus aurata",
             "Atherina boyeri", "Dicentrarchus labrax", "Liza ramada", "Sprattus sprattus")

#Calculate average biomass per year of the most common fishes cached
fishes_spp <- catches_clean %>%
  filter(!is.na(Kgs)) %>% #remove NANs
  filter(Species %in% include) %>% #select species to be modeled
  group_by(Year2) %>%
  summarise(mean_biomass = mean(Kgs,na.rm=T)) %>%
  mutate(mean_biomass_log = log(mean_biomass)) %>%
  mutate(Year=Year2) %>%
  as.data.frame()


# subset river flow from CHE AND other physico-chemical variables, and JOIN with fish catches
data_full <- data_che_Tortosa %>% select(Day, Month, Year, CaudalÂ.enÂ.superficieÂ..m3.s.,
                                         FÃ³sforoÂ.TotalÂ..mg.LÂ.P.,FosfatosÂ..mg.LÂ.PO4.,
                                         TemperaturaÂ.delÂ.aguaÂ..ÂºC.,
                                         OxÃ.genoÂ.disueltoÂ..mg.LÂ.O2.,
                                         MateriasÂ.enÂ.suspensiÃ³nÂ..mg.L.) %>%
  rename(flow_che=CaudalÂ.enÂ.superficieÂ..m3.s.) %>%
  rename(Total_phosphorous=FÃ³sforoÂ.TotalÂ..mg.LÂ.P.) %>%
  rename(SRP=FosfatosÂ..mg.LÂ.PO4.) %>%
  rename(WaterT=TemperaturaÂ.delÂ.aguaÂ..ÂºC.) %>%
  rename(Oxygen=OxÃ.genoÂ.disueltoÂ..mg.LÂ.O2.) %>%
  rename(SuspendedSolids=MateriasÂ.enÂ.suspensiÃ³nÂ..mg.L.) %>%
  full_join(fishes_spp, by = 'Year') %>% #here join with mean fishes catches across lagoons
  left_join(data_cat[c("Year", "Month", "Chl.Total")], by=c("Year", "Month")) %>% #here join chl data from CAT dataset
  full_join(hist_flow[c("Year", "Month", "qmedmes")], by=c('Year', 'Month')) #here historical flow from CHE dataset

  
# plot the data
#plotting the data
ggplot(data=gather(data_full, variable, value, -Year, -Year2, -Day, -Month), 
       aes(x=Year, 
           y=value, 
           group=variable)) + 
  geom_line() + 
  facet_wrap("variable",scales = "free_y", ncol = 2) +
  xlab("Years") +
  ylab("") +
  ggtitle("") +
  theme_bw()


#save the dataset
write.csv(data_full, "outputs/river_fisheries_data.csv")

## Generate lagged datasets
# Prepare data

# set the driver variable
"qmedmes" #river flow
"Chl.Total" #Chla
"SRP" #phosphorous

env <- data_full %>% group_by(Year) %>%
  summarise_at(c("Chl.Total", "SRP", "Total_phosphorous","flow","WaterT","Oxygen","SuspendedSolids"), mean, na.rm = TRUE) %>%
  #summarise(flow=mean(SRP, na.rm=TRUE)) %>%
  filter(!is.na(environment)) %>% #remove NANs
  as.data.frame()

# set the response variable
catches <- data_full %>% group_by(Year) %>%
  summarise(biomass=mean(mean_biomass)) %>%
  filter(!is.na(biomass)) %>% #remove NANs
  as.data.frame() %>%
  left_join(env, by="Year") %>%
  filter(!is.na(Chl.Total))

#plotting the data
ggplot(data=gather(catches, variable, value, -Year), 
       aes(x=Year, 
           y=value, 
           group=variable)) + 
  geom_line() + 
  facet_wrap("variable",scales = "free_y") +
  xlab("Years") +
  ylab("") +
  ggtitle("") +
  theme_bw()


# create two vectors with year~response and year~driver
env <- catches[,c(1,6)]
catches <- catches[,c(1,2)]

## Fix the non regular time series 
env <- catches[,c(1,3)]
env$Year <- seq(1, length(env$Year), by=1)
catches <- catches[,c(1,2)]
catches$Year <- seq(1, length(catches$Year), by=1)


## Backward lags
lags<-1:10 #change length of lags depending on each driver

#backward dataset 
#to assess the effect of “past” environment (e.g. flow; data.to.lag) on fish catches (reference.data)
lag.data.backward <- backwardLags(
  lags=lags, 
  reference.data=catches, 
  data.to.lag=env
)

#preparing plotting data (4 lags only)
temp.backward <- lag.data.backward
temp.backward$lag <- as.numeric(temp.backward$lag)
temp.backward <- temp.backward[temp.backward$lag%in% c(1, 3, 6, 10),]

plot.past <- ggplot(data=temp.backward, aes(x=environment, y=biomass, group=lag)) + 
  geom_point(shape=21, fill="gray50", color="black", size=2, alpha=0.5) +
  facet_wrap("lag", ncol=1) +
  xlab("") +
  ylab("") +
  ggtitle("") +
  theme(text=element_text(size=12),
        plot.title=element_text(size = 16), legend.position="none") +
  geom_smooth(method = lm, size=2, color="red4", se=FALSE, aes(alpha=0.5))

plot.past


#fitting a GLS model per lag on backward datasets
backward.results <- modelLagData(
  model.formula="biomass ~ environment", 
  lagged.data=lag.data.backward
)

backward.results$value <- round(backward.results$value, 2)

# fitting a null model
backward.results.random <- modelRandomLagData(
  lagged.data=lag.data.backward, 
  model.formula="biomass ~ environment", 
  iterations=1000
)

backward.results.random$value <- round(backward.results.random$value, 2)


# plot model results with own code
#axes limits
max.lag <- max(c(backward.results$lag))
max.coefficient <- round(max(c(backward.results[backward.results$variable=="Coefficient", "value"], backward.results.random[backward.results.random$variable=="Coefficient", "upper"])) + 0.1, 1)
min.coefficient <- round(min(c(backward.results[backward.results$variable=="Coefficient", "value"], backward.results.random[backward.results.random$variable=="Coefficient", "lower"])) - 0.1, 1)
max.R2 <- round(max(c(backward.results[backward.results$variable=="R2", "value"])), 1)

library(viridis)
viridis.colors <- viridis(10, option="D")

backward.plot.coefficient <- ggplot(data=subset(backward.results, variable=="Coefficient"), aes(x=lag, y=value)) +
  geom_ribbon(data=subset(backward.results.random, variable=="Coefficient"), aes(ymin=lower, ymax=upper), alpha=0.3, fill="light grey") +
  geom_line(data=subset(backward.results.random, variable=="Coefficient"), aes(x=lag, y=value), alpha=0.6, color="light grey", size=1) +
  geom_hline(yintercept=0, color="black", linetype=2) +
  geom_ribbon(aes(ymin=lower,ymax=upper), alpha=0.3, fill=viridis.colors[2]) +
  geom_line(size=1.5, color=viridis.colors[1]) +
  ggtitle(expression("River flow" %->% "Catches (biomass)")) +
  theme(legend.position="none") +
  xlab("") +
  ylab("Standardized coefficient") +
  scale_y_continuous(breaks=seq(min.coefficient, max.coefficient, by=0.8)) +
  scale_x_reverse()+
  theme(axis.text = element_text(size=12),
        plot.title = element_text(size=14),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x = element_blank()) +
  theme_classic()
#coord_cartesian(ylim = c(min.coefficient, max.coefficient + 0.5))

backward.plot.coefficient


backward.plot.R2 <- ggplot(data=subset(backward.results, variable=="R2"), aes(x=lag, y=value, group=variable)) +
  geom_ribbon(data=subset(backward.results.random, variable=="R2"), aes(ymin=lower, ymax=upper), alpha=0.3, fill="light grey") +
  geom_line(data=subset(backward.results.random, variable=="R2"), aes(x=lag, y=value), alpha=0.6, color="light grey", size=1.5) +
  geom_line(size=1.5, color=viridis.colors[2]) +
  theme(legend.position="none") +
  xlab("Years (before fish catches)") +
  ylab("Pseudo R squared") +
  scale_y_continuous(breaks=seq(0, max.R2, by=0.1)) +
  scale_x_reverse()+
  theme(axis.text = element_text(size=12),
        axis.text.x = element_text(size=12),
        plot.title = element_text(size = 14),
        plot.margin = unit(c(0.2, 0.5, 0, 0), "cm")) +
  theme_classic() +
  coord_cartesian(ylim = c(0, max.R2 + 0.05))

backward.plot.R2

#Combine plots
plot_composite <- plot_grid(backward.plot.coefficient, 
                            backward.plot.R2, ncol = 1, 
                            rel_heights = c(1, 1), align="v") + 
  theme(plot.margin = unit(c(0.5, -1, 0.5, 0.5), "cm"))
plot_composite

# save plot
ggsave("outputs/riverflow_catches_asycnchronousModel.png",
      plot = plot_composite,
      width=8,
      height=6,
      units="in",
      dpi = 400)

### 
# Run Lagged models with ecological memory (https://github.com/BlasBenito/memoria)
###
library(memoria)

# Prepare time lagged data
biomass.env.lagged <- prepareLaggedData(
  input.data = catches,
  response = "biomass",
  drivers = c("Chl.Total", "SRP","flow"),
  time = "Year",
  oldest.sample = "last",
  lags = seq(1, 5, by=1),
  time.zoom=NULL,
  scale=FALSE
)

#computing memory
memory.output <- computeMemory(
  lagged.data = biomass.env.lagged,
  drivers = c("Chl.Total", "SRP","flow"),
  response = "Response",
  add.random = TRUE,
  random.mode = "white.noise",
  repetitions = 100,
)

plotMemory(memory.output)
