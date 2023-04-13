##Long-term coastal lagoon fisheries landings

#Load libraries
library(tidyverse)
library(rfishbase)

# Read in the data of species biomass per year and lagoon (1965-2010)
catches_raw <- read.csv("data/catches_lagoons_years_raw.csv", sep = ";") #from the original file "Dades Captures.xlsx"
str(catches_raw)
range(catches_raw$Year2)

# Updated fish catches dataset (1965-2021)
#catches_raw <- read.csv("data/catches_lagoons_years_raw_1965_2020.csv", sep = ";") #from the original file "Dades Captures.xlsx"

#Rename some columns
colnames(catches_raw)[1] <- "Name"

#transform relevant columns to factors or numeric
catches_raw$Name <- as.factor(catches_raw$Name)
catches_raw$Common_name <- as.factor(catches_raw$Common_name)
catches_raw$Stage <- as.factor(catches_raw$Stage)
catches_raw$Group <- as.factor(catches_raw$Group)
catches_raw$Order <- as.factor(catches_raw$Order)
catches_raw$Family <- as.factor(catches_raw$Family)
catches_raw$Species <- as.factor(catches_raw$Species)
catches_raw$Lagoon <- as.factor(catches_raw$Lagoon)
#catches_raw$Year <- as.factor(catches_raw$Year)
#catches_raw$Year2 <- as.factor(catches_raw$Year2) #not sure if leave it as integer
catches_raw$Kgs <- as.numeric(gsub(",", ".", gsub("\\.", "", catches_raw$Kgs))) #replace commas with dots for decimals

#Rename some columns
colnames(catches_raw)[1] <- "Name"
colnames(catches_raw)[2] <- "commonName"

# extract traits using Cano-Barbacil 2019 FBW and average values across databases
traitsDB <- read.csv("data/Fish_trait_raw_data_SI.csv", sep = ";") %>%
  group_by(Species) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

# rfishbase function
# traits <- species(unique(catches_raw$Species))
# colnames(traits)
#fields=c("Species", "Resilience", "StockDefs"))

#Save dataset for posterior analysis
#write.csv(catches_raw, "outputs/catches_clean.csv")

#### Read in FQ data
## CHE data
data_che <- read.csv("data/data_che_full.csv", sep=";")
head(data_che)
str(data_che)

# transform to numeric
df1 <- data.frame(apply(apply(data_che[,7:ncol(data_che)], 2, gsub, patt=",", replace="."), 2, as.numeric))
data_che <- cbind(data_che[,1:6], df1) 
data_che_Tortosa <- data_che %>% filter(Sampling_Point==2) %>% #filter cases Tortosa station
  dplyr::select(c(3,4,7,8,12,15,18,19,21,22,23,24,32)) %>%
  rename(alkalinity=3) %>%
  rename(ammonia_che=4) %>%
  rename(flow_che=5) %>%
  rename(DBO_che=6) %>%
  rename(SRP_che=7) %>%
  rename(total_phosphorous_che=8) %>%
  rename(SuspendedSolids_che=9) %>%
  rename(Nitrates_che=10) %>%
  rename(Nitrites_che=11) %>%
  rename(total_nitrogen_che=12) %>%
  rename(WaterT_che=13)

#Calculate mean annual flow Data CHE (1980-2004) to match with fish catches
year_mean_flow_che <- data_che_Tortosa %>%
  mutate(Year_f=factor(Year)) %>%
  rename(flow=12) %>% #12 is the column position 
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

plot(hist_flow$Year, hist_flow$qmedmes)


## Read in climatic data (Tortosa station) 1910-2013
climatic_data_1910_2013 <- read.csv("data/data_climatic_tortosa_1910_2013.csv")[-1] %>% #remove first column (date)
  rename(Day=DAY) %>%
  rename(Month=MONTH) %>%
  rename(Year=YEAR) %>%
  rename(Wind_direction=5) %>%
  rename(Wind_speed=6) %>%
  rename(Sea_level_pressure=9) %>%
  rename(Precipitation=10) %>%
  rename(Mean_temperature=12)%>%
  rename(Min_temperature=13) %>%
  rename(Max_temperature=14) %>%
  dplyr::select(-c(4,7,8,11))

str(climatic_data_1910_2013)  
summary(climatic_data_1910_2013)
plot(climatic_data_1910_2013$Year, climatic_data_1910_2013$Mean_temperature)

## Read in climatic data (Tortosa station) 2000-2020
climatic_data_2000_2020 <- read.csv("data/data_climatic_tortosa_2000_2020.csv", sep = ";") %>%
  rename(Month=month) %>%
  rename(Year=year) %>%
  dplyr::select(1,2,14,20,26,27,35,26)

str(climatic_data_2000_2020)

# transform to numeric
df1 <- data.frame(apply(apply(climatic_data_2000_2020[,4:ncol(climatic_data_2000_2020)], 2, gsub, patt=",", replace="."), 2, as.numeric))
climatic_data_2000_2020 <- cbind(climatic_data_2000_2020[,1:3], df1) 
summary(climatic_data_2000_2020)
plot(climatic_data_2000_2020$Year, climatic_data_2000_2020$p_mes)

## Read in flow CHE data (1980-2007)
flow_che_1980_2007 <- read.csv("data/flow_che.csv", sep = ";") %>%
  rename(Month=ï..Month)
str(flow_che_1980_2007)

# transform to numeric
df1 <- data.frame(apply(apply(flow_che_1980_2007[,3:ncol(flow_che_1980_2007)], 2, gsub, patt=",", replace="."), 2, as.numeric))
flow_che_1980_2007 <- cbind(flow_che_1980_2007[,1:2], df1) 
summary(flow_che_1980_2007)
plot(flow_che_1980_2007$Year, flow_che_1980_2007$QMean)

## Read in flow CHE data (1997-2020)
flow_che_1997_2020 <- read.csv("data/flow_che_1997_2020.csv", sep = ";") %>%
  rename(Month=ï..month) %>%
  rename(Year=year)


#Species to be included in dataset (>30% average biomass)
include <- c("Anguilla anguilla", "Cyprinus carpio", "Mullet ind.", "Sparus aurata",
             "Atherina boyeri", "Dicentrarchus labrax", "Liza ramada", "Sprattus sprattus", "Callinectes sapidus")

#Calculate average biomass per year of the most common fishes catched
fishes_spp <- catches_raw %>%
  filter(!is.na(Kgs)) %>% #remove NANs
  filter(Species %in% include) %>% #select species 
  group_by(Year2, Species, Lagoon) %>% #group by species, lagoon and year
  summarise(mean_biomass = mean(Kgs,na.rm=T)) %>%
  mutate(mean_biomass_log = log(mean_biomass)) %>%
  mutate(Year=Year2) %>%
  as.data.frame() %>%
  ungroup()

# subset river flow from CHE AND other physico-chemical variables (CAT) and climatic data, and JOIN with fish catches
data_full <- fishes_spp %>%
  left_join(data_che_Tortosa, by = "Year") %>%
  left_join(hist_flow[c("Year","qmedmes")], by=c('Year')) %>% #here historical flow from CHE dataset
  left_join(data_cat[c("Year", "Month", "Chl.Total", "TOC", "NO3", "NO2", "NH4")], by=c("Year", "Month")) %>% #here join chl data from CAT dataset
  left_join(climatic_data_1910_2013[c("Year", 'Month', 'Wind_direction', 'Wind_speed', 'Sea_level_pressure', 'Precipitation', 'Mean_temperature', 'Min_temperature', 'Max_temperature')], by=c('Year', 'Month')) %>% #here climatic data (Tortosa station)
  left_join(flow_che_1980_2007, by=c("Year", "Month")) %>%
  filter(Year>=1964)

# subset river flow from CHE AND other physico-chemical variables (CAT) and climatic data
data_FQ_climatic <- data_che_Tortosa %>%
  left_join(hist_flow[c("Year","qmedmes")], by=c('Year')) %>% #here historical flow from CHE dataset
  left_join(data_cat[c("Year", "Month", "Chl.Total", "TOC", "NO3", "NO2", "NH4")], by=c("Year", "Month")) %>% #here join chl data from CAT dataset
  left_join(climatic_data_1910_2013[c("Year", 'Month', 'Wind_direction', 'Wind_speed', 'Sea_level_pressure', 'Precipitation', 'Mean_temperature', 'Min_temperature', 'Max_temperature')], by=c('Year', 'Month')) %>% #here climatic data (Tortosa station)
  left_join(flow_che_1980_2007, by=c("Year", "Month")) %>%
  filter(Year>=1964)

str(data_full)
range(data_full$Year)
range(data_full$Month)

var_plt <- c("Year", "flow_che", "Chl.Total", "total_phosphorous_che")

## quick plot the data
ggplot(data=gather(data_FQ_climatic, variable, value, -Year) %>%
         filter(variable %in% var_plt), 
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
#write.csv(data_full, "outputs/river_fisheries_data.csv")
#write.csv(data_full, "outputs/river_fisheries_lagoon_spp_data.csv")



