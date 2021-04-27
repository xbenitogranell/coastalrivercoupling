## Exploratory analysis: time series plots

#Load functions used
library(tidyverse)
library(ggplot2)

#Read in catches lagoons clean data
catches_clean <- read.csv("data/catches_clean.csv")[-1]
str(catches_clean)

#Calculate average biomass per species per year and lagoon
year_mean_biomass <- catches_clean %>%
  filter(!is.na(Kgs)) %>% #remove NANs
  filter(!str_detect(Species, "Ind")) %>% #drop Undetermined species
  group_by(Species,Year2,Lagoon) %>% 
  summarise(average_biomass=mean(Kgs,na.rm = T)) %>%
  mutate(average_biomass=round(average_biomass,2)) %>% #place 2 decimals for biomass
  ungroup() %>%
  group_by(Year2, Lagoon) %>%
  mutate(percentage = average_biomass / sum(average_biomass) * 100)#calculate mean percentage


#Plot temporal trends in biomass per lagoon
plt <- ggplot(year_mean_biomass, aes(x = Year2, y = average_biomass)) +
  geom_line() +
  facet_wrap(Lagoon~., scales = "free") +
  theme_bw()+ theme(legend.position = "none")+
  labs(y = "", x = "", title = "")
plt

#subset spp having mean average percentage > 20
subset <- year_mean_biomass %>% filter(percentage>20)

plt_spp <- ggplot(subset, aes(x=Year2, y=average_biomass, colour=Species))+
  geom_line()+
  facet_wrap(Lagoon~., scales = "free") +
  theme_bw()+ theme(legend.position = "bottom")+
  labs(y = "", x = "", title = "")
plt_spp

#Plot temporal trend in % of total assemblage of fishes per lagoon
ggplot(subset, aes(fill=Species, y=average_biomass, x=Year2)) + 
  facet_wrap(Lagoon~., scales = "free") +
  geom_bar(position="fill", stat="identity")+
  theme_bw() +
  labs(y="% total assemblage", x="Years")



