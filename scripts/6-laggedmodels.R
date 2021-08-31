## EXPLORATORY ANALYSIS OF THE EFFECT OF LAGGED RIVER DYNAMICS ON COASTAL FISHERIES 

source("scripts/functions_lagged.R")

#Clear workspace
rm(list=ls(all=TRUE))

#Load functions used
library(tidyverse)
library(ggplot2)
library(cowplot)

# Read in CHE data
data_che_raw <- read.csv("data/data_che_full.csv", sep=";", dec = ".")
str(data_che_raw)

#replace commas with dots for decimals, replace NAs with 0, round to 2 decimals, and bind columns
data_che <- as.data.frame(apply(apply(data_che_raw[,7:ncol(data_che_raw)], 2, gsub, patt=",", replace="."), 2, as.numeric))
data_che <- data_che %>% mutate_all(~replace(., is.na(.), 0)) %>% mutate(across(where(is.numeric), round, 2))
data_che <- cbind(data_che_raw[,c(2:6)], data_che) %>%
  rename_with(~ tolower(gsub(".", "_", .x, fixed = TRUE))) 


#Read in CAT data
data_cat_raw <- read.csv("data/data_cat_full.csv", sep=";", dec = ".")
str(data_cat_raw)

# Remove 2nd row (units)
data_cat_raw <- data_cat_raw[-1,]

#replace commas with dots for decimals, replace NAs with 0, round to 2 decimals, and bind columns
data_cat <- as.data.frame(apply(apply(data_cat_raw[,3:ncol(data_cat_raw)], 2, gsub, patt=",", replace="."), 2, as.numeric))
data_cat <- data_cat %>% mutate_all(~replace(., is.na(.), 0)) %>% mutate(across(where(is.numeric), round, 2))
data_cat <- cbind(data_cat_raw[,c(1:2)], data_cat) %>%
  mutate(year=Year)


#Prepare fish catches data
#Read in catches lagoons clean data
catches_clean <- read.csv("data/catches_clean.csv")[-1]
str(catches_clean)

#Species to be included in the models (>30% average biomass)
include <- c("Anguilla anguilla", "Cyprinus carpio", "Mullet ind.", "Sparus aurata",
             "Atherina boyeri", "Dicentrarchus labrax", "Liza ramada", "Sprattus sprattus")

#Calculate average biomass per year of the most common fishes cacthed
fishes_spp <- catches_clean %>%
  filter(!is.na(Kgs)) %>% #remove NANs
  filter(Species %in% include) %>% #select species to be modeled
  group_by(Year2) %>%
  summarise(mean_biomass = mean(Kgs,na.rm=T)) %>%
  mutate(mean_biomass_log = log(mean_biomass)) %>%
  mutate(year=Year2) %>%
  as.data.frame()


# subset river flow AND other physico-chemical variables, and JOIN with fish catches
data_full <- data_che %>% select(day, month, year, caudalâ_enâ_superficieâ__m3_s_,
                                 fã³sforoâ_totalâ__mg_lâ_p_,fosfatosâ__mg_lâ_po4_) %>%
  rename(flow=caudalâ_enâ_superficieâ__m3_s_) %>%
  rename(Total_phosphorous=fã³sforoâ_totalâ__mg_lâ_p_) %>%
  rename(SRP=fosfatosâ__mg_lâ_po4_) %>%
  right_join(fishes_spp, by = 'year') %>%
  right_join(data_cat[c("year","Chl.Total")], by='year')   #here join chl data from CAT dataset
  
#save the dataset
#write.csv(data_full, "outputs/river_fisheries_data.csv")

## Generate lagged datasets
# Prepare data
env <- data_full %>% group_by(year) %>%
  summarise(environment=mean(flow)) %>%
  as.data.frame()

catches <- data_full %>% group_by(year) %>%
  summarise(biomass=mean(mean_biomass_log)) %>%
  as.data.frame()

## Backward lags
lags<-1:10

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
#ggsave("outputs/agropastoralism_diatoms_asycnchronousModel.png",
#       plot = plot_composite,
#       width=8,
#       height=6,
#       units="in",
#       dpi = 400)

