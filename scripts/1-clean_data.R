##Long-term coastal lagoon fisheries landings

#Load libraries
library(tidyverse)

# Read in the data of species biomass per year and lagoon
catches_raw <- read.csv("data/catches_lagoons_years_raw.csv", sep = ";") #from the original file "Dades Captures.xlsx"
str(catches_raw)

#transform relevant columns to factors or numeric
catches_raw$ï..Name <- as.factor(catches_raw$ï..Name)
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

#Save dataset for posterior analysis
#write.csv(catches_raw, "data/catches_clean.csv")

# Assign new name to the dataframe
catches_clean <- catches_raw

#### Read in environmental data
## CHE data
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

## Read in climatic data (Tortosa station)
climatic_data <- read.csv("data/data_climatic_tortosa.csv")[-1] %>% #remove first column (date)
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

str(climatic_data)  

## Read in flow data (CHE)
flow_che <- read.csv("data/flow_che.csv", sep = ";") %>%
  rename(Month=ï..Month)
str(flow_che)

# transform to numeric
df1 <- data.frame(apply(apply(flow_che[,3:ncol(flow_che)], 2, gsub, patt=",", replace="."), 2, as.numeric))
flow_che <- cbind(flow_che[,1:2], df1) 

#Read in catches lagoons clean data
catches_clean <- read.csv("data/catches_clean.csv")[-1]
str(catches_clean)

#Species to be included in dataset (>30% average biomass)
include <- c("Anguilla anguilla", "Cyprinus carpio", "Mullet ind.", "Sparus aurata",
             "Atherina boyeri", "Dicentrarchus labrax", "Liza ramada", "Sprattus sprattus")

#Calculate average biomass per year of the most common fishes catched
fishes_spp <- catches_clean %>%
  filter(!is.na(Kgs)) %>% #remove NANs
  filter(Species %in% include) %>% #select species to be modeled
  group_by(Year2, Species, Lagoon) %>% #group by species, lagoon and year
  #group_by(Year2) %>% #group by species, lagoon and year
  summarise(mean_biomass = mean(Kgs,na.rm=T)) %>%
  mutate(mean_biomass_log = log(mean_biomass)) %>%
  mutate(Year=Year2) %>%
  as.data.frame()

# subset river flow from CHE AND other physico-chemical variables (CAT) and climatic data, and JOIN with fish catches
data_full <- data_che_Tortosa %>% dplyr::select(c(2,3,4,7,8,12,15,18,19,21,22,23,24,32)) %>%
  rename(alkalinity=4) %>%
  rename(ammonia_che=5) %>%
  rename(flow_che=6) %>%
  rename(DBO_che=7) %>%
  rename(SRP_che=8) %>%
  rename(total_phosphorous_che=9) %>%
  rename(SuspendedSolids_che=10) %>%
  rename(Nitrates_che=11) %>%
  rename(Nitrites_che=12) %>%
  rename(total_nitrogen_che=13) %>%
  rename(WaterT_che=14) %>%
  full_join(fishes_spp, by = 'Year') %>% #here join with mean fishes catches across lagoons
  left_join(data_cat[c("Year", "Month", "Chl.Total", "TOC", "NO3", "NO2", "NH4")], by=c("Year", "Month")) %>% #here join chl data from CAT dataset
  full_join(hist_flow[c("Year", "Month", "qmedmes")], by=c('Year', 'Month')) %>% #here historical flow from CHE dataset
  full_join(climatic_data[c("Year", 'Month', 'Wind_direction', 'Wind_speed', 'Sea_level_pressure', 'Precipitation', 'Mean_temperature', 'Min_temperature', 'Max_temperature')], by=c('Year', 'Month')) %>% #here climatic data (Tortosa station)
  full_join(flow_che, by=c("Year", "Month"))

head(data_full)

## quick plot the data
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
#write.csv(data_full, "outputs/river_fisheries_data.csv")
write.csv(data_full, "outputs/river_fisheries_lagoon_spp_data.csv")



