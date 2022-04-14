## EXPLORATORY ANALYSIS OF THE EFFECT OF LAGGED RIVER DYNAMICS ON COASTAL FISHERIES 

#Clear workspace
rm(list=ls(all=TRUE))

#load lag functions
source("scripts/functions/functions_lags.R") #from Blas Benito

#Load functions used
library(tidyverse)
library(ggplot2)
library(cowplot)
library(nlme)

# Read in data full
data_full <- read.csv("outputs/river_fisheries_lagoon_spp_data.csv")
str(data_full)

### Generate lagged datasets
# Prepare data
# set the driver variable
"qmedmes" #river flow
"Chl.Total" #Chla
"SRP" #phosphorous

env <- data_full %>% group_by(Year) %>%
  summarise_at(c("Chl.Total", "SRP_che", "total_phosphorous_che","flow_che", "qmedmes", "WaterT_che","SuspendedSolids_che", "TOC",
                 "NO3", "NO2", "NH4", "Wind_speed", "Mean_temperature", "Precipitation"), mean, na.rm = TRUE) %>%
  filter(!is.na(environment)) %>% #remove NANs
  as.data.frame()

# set the response variable
catches <- data_full %>% group_by(Year) %>%
  summarise(biomass=mean(mean_biomass)) %>%
  filter(!is.na(biomass)) %>% #remove NANs
  as.data.frame() %>%
  left_join(env, by="Year") 
  #filter(!is.na(Chl.Total))

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
lags<-1:30 #change length of lags depending on each driver

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
