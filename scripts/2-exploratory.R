## Exploratory analysis: time series plots

#Clear workspace
rm(list=ls(all=TRUE))

#Load functions used
library(tidyverse)
library(ggplot2)

#Read in catches lagoons clean data
catches_clean <- read.csv("outputs/catches_clean_1965_2010.csv")[-1]
str(catches_clean)

range(catches_clean$Year2)

#Calculate average biomass per year and lagoon
year_mean_biomass <- catches_clean %>%
  filter(!is.na(Kgs)) %>% #remove NANs
  filter(!str_detect(Species, "Ind")) %>% #drop Undetermined species
  group_by(Year2,Lagoon) %>% 
  summarise(average_biomass=mean(Kgs,na.rm = T)) %>%
  mutate(average_biomass=round(average_biomass,2)) 

#Plot temporal trends in biomass per lagoon
plt_lagoon <- ggplot(year_mean_biomass, aes(x = Year2, y = average_biomass)) +
  geom_line() +
  facet_wrap(Lagoon~., scales = "free") +
  theme_bw()+ theme(legend.position = "none")+
  labs(y = "Biomass (Kgs)", x = "Years", title = "")
plt_lagoon

ggsave("outputs/lagoon_mean_biomass.png",
       plot = plt_lagoon,
       width=8,
       height=6,
       units="in",
       dpi = 400)

#subset spp having mean average percentage > 20
species_mean_biomass <- catches_clean %>%
  filter(!is.na(Kgs)) %>% #remove NANs
  filter(!str_detect(Species, "Ind")) %>% #drop Undetermined species
  group_by(Species, Year2,Lagoon) %>% 
  summarise(average_biomass=mean(Kgs,na.rm = T)) %>%
  mutate(average_biomass=round(average_biomass,2)) %>% #place 2 decimals for biomass
  ungroup() %>%
  group_by(Year2, Lagoon) %>%
  mutate(percentage = average_biomass / sum(average_biomass) * 100) %>% #calculate mean percentage 
  filter(percentage>20)

plt_spp <- ggplot(species_mean_biomass, aes(x=Year2, y=average_biomass, colour=Species))+
  geom_line()+
  facet_wrap(Lagoon~., scales = "free") +
  theme_bw()+ theme(legend.position = "bottom")+
  labs(y = "Biomass (Kgs)", x = "Years", title = "")
plt_spp

ggsave("outputs/species_mean_biomass.png",
       plot = plt_spp,
       width=8,
       height=6,
       units="in",
       dpi = 400)

#Plot temporal trend in % of total assemblage of fishes per lagoon
plt_proportion <- ggplot(species_mean_biomass, aes(fill=Species, y=average_biomass, x=Year2)) + 
  facet_wrap(Lagoon~., scales = "free") +
  geom_bar(position="fill", stat="identity")+
  theme_bw() +
  labs(y="% total assemblage", x="Years")
plt_proportion

ggsave("outputs/species_proportion.png",
       plot = plt_proportion,
       width=8,
       height=6,
       units="in",
       dpi = 400)


