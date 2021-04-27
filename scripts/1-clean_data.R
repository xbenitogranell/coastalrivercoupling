##Long-term coastal lagoon fisheries landings

#Read raw data and clean for intermediate data analyses

#Load libraries
library(tidyverse)

# Read in the data of species biomass per year and lagoon
catches_raw <- read.csv("data/catches_lagoons_years_clean.csv", sep = ";") #from the original file "Dades Captures.xlsx"
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
write.csv(catches_raw, "data/catches_clean.csv")


