## Model temporal trends of fish biomass

#Clear workspace
rm(list=ls(all=TRUE))

# load functions used
library(tidyverse)
library(ggplot2)
library(mgcv) #to perform GAMs
library(gratia) # to perform diagnostic plots of GAM objects
library(car)

#Read in catches lagoons clean data
catches_clean <- read.csv("data/catches_clean.csv")[-1] 
str(catches_clean)

#Species to be included in the models (>30% average biomass)
include <- c("Anguilla anguilla", "Cyprinus carpio", "Mullet ind.", "Sparus aurata",
             "Atherina boyeri", "Dicentrarchus labrax", "Liza ramada", "Sprattus sprattus")
exclude <- c("Platjola", "Clot - Baseta")

## This is to fit a GAM on temporal total biomass
#Calculate average biomass per year and lagoon
fishes_comm <- catches_clean %>%
  filter(!is.na(Kgs)) %>% #remove NANs
  filter(!str_detect(Species, "Ind")) %>% #drop Undetermined species
  filter(Species %in% include) %>% #select species to be modeled
  filter(!Lagoon %in% exclude) %>%
  group_by(Year2) %>%
  summarise(mean_biomass = mean(Kgs,na.rm=T)) %>%
  mutate(mean_biomass_log = log(mean_biomass)) %>%
  mutate(year_f=factor(Year2)) %>%
  mutate(Year=Year2) %>%
  as.data.frame()
  
# Global model
set.seed(10) #allow replication of results
model_gam <- gam(mean_biomass ~ s(Year2, k=40),
                 data=fishes_comm, family = Gamma(link = "log"),  
                 method = "REML", select=TRUE)

# Example of model with a random factor ("bs="re") (Not run)
mod1 <- gam(mean_biomass ~ s(Year2, k=40) + s(Lagoon, bs="re"),
                  data=catches_clean, family = Gamma(link = "log"),  
                  method = "REML", select=TRUE)


#check gam outputs
appraise(model_gam)
gam.check(model_gam)
appraise(model_gam)
draw(model_gam)
summary(model_gam)

pacf(residuals(model_gam)) #check autocorrelation


#controlling for temporal autocorrelation
library(lme4)
model_gamm <- gamm(mean_biomass_log ~ s(Year2, k=20),
                  data=fishes_comm, family = gaussian(link = "identity"),  
                  method = "REML", correlation=corCAR1(form = ~Year2))

#Compare AIC of the two models
AIC_table <- AIC(model_gam, model_gamm$lme)%>%
  rownames_to_column(var= "Model")%>%
  mutate(data_source = rep(c("diatom_data")))%>%
  group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  dplyr::select(-data_source)%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0)))
AIC_table


#Create list of GAM fit and observation data
full_model_data <- list(fishes_comm)
fits <- list(model_gam)
names(fits) <- c("global_model")

# Fit GAM model on temporal biomass for each lagoon separately
#Calculate average biomass per year and lagoon
fishes_comm_lagoon <- catches_clean %>%
  filter(!is.na(Kgs)) %>% #remove NANs
  filter(!str_detect(Species, "Ind")) %>% #drop Undetermined species
  filter(Species %in% include) %>% #select species to be modeled
  group_by(Year2, Lagoon) %>%
  summarise(mean_biomass = mean(Kgs,na.rm=T)) %>%
  mutate(mean_biomass_log = log(mean_biomass)) %>%
  rename(lake=Lagoon) %>%
  filter(!lake=="Platjola" & !lake=="Clot - Baseta") %>% 
  as.data.frame()

levels(fishes_comm_lagoon$lake)

# Create list of lagoon datasets
LagoonData <- split(fishes_comm_lagoon, fishes_comm_lagoon$Lagoon)

## fit mean biomass across lagoons  
fitGam <- function(i, data, kk, select=TRUE) {
  Ldata <- data[[i]]
  k <- kk[i]
  fit <- gam(mean_biomass_log ~ s(Year2, k = k),
             data = Ldata, method = "REML",
             family = gaussian(link = "identity"))
}

lagoons <- c("Encanyissada", "Tancada", "Canal vell", "Goleta")
k <- c(rep(30,4))
fits <- lapply(seq_along(LagoonData[lagoons]), fitGam, data = LagoonData[lagoons], kk = k)
names(fits) <- names(LagoonData[lagoons])


#Create a function the extract GAM results
GAM_results <- function(i, models){
  model <- models[[i]]
  formula <- paste("", format(model$formula))
  n <- as.numeric(anova(model)["n"]) #take n samples
  family.function <- model$family[1] #take family
  family.link <- model$family[2] #take link
  fit.deviance <- model$deviance #take model deviance
  null.deviance <- model$null.deviance #take null model deviance
  r.sq <- anova(model)["r.sq"] #take r square 
  s.table <- anova(model)["s.table"] #take summary table of smooth terms 
  
  table <- cbind.data.frame(formula,n,family.function,family.link,fit.deviance,null.deviance,r.sq,s.table) #combine extracted columns
  colnames(table) <- c("formula","n","family.function", "family.link", "deviance", "null.deviance","r.sq",
  "s.edf", "s.ref.df", "p") #
  return(table)

}

summaryGam <- lapply(seq_along(fits), GAM_results, models=fits)
names(summaryGam) <- names(fits)

#extract dataframes from the summary GAM list and save table
summary_Gam_tables <- plyr::ldply(summaryGam, data.frame) 
write.table(summary_Gam_tables, file = "outputs/summary_gam_lagoons_table.txt", 
            sep = "\t", na = "", col.names = TRUE, row.names = FALSE)


#Load function to calculate first derivative on fish community biomass
source("scripts/functions/Deriv.R")

## Apply derivative function to global model (list)
derivs <- lapply(fits, Deriv, n = 200)
Term <- "Year2"
confints <- lapply(derivs, confint, term = Term)

#Create dataset to predict model fitted values
makePredData <- function(data, n = 200) {
  data.frame(Year2 = seq(min(fishes_comm$Year2), max(fishes_comm$Year2), length.out = n))
}

#wrap-up the function
predData <- lapply(fits, makePredData)

#Predict results of gam model  
modelPreds <- function(i, models, newdata, se.fit = TRUE) {
  predict(models[[i]], newdata = newdata[[i]], se.fit = se.fit, type="response")
}

#wrap-up the function
predTrends <- lapply(seq_along(fits), modelPreds,
                     models = fits, newdata = predData)

names(predTrends) <- names(fits)

#Create a function for calculating first derivatives
signifWrap <- function(i, data, derivObj, term, ciObj) {
  signifD(data[[i]]$fit,
          d = derivObj[[i]][[Term]]$deriv,
          ciObj[[i]][[Term]]$upper,
          ciObj[[i]][[Term]]$lower)
}

#wrap-up the function
signifCores <- lapply(seq_along(fits), signifWrap, data = predTrends,
                      derivObj = derivs, term = Term, ciObj = confints)

names(signifCores) <- names(fits)

# Function to plot GAM derivative trends
plotTrends <- function(i, trends, dates, signifs, obs, cex = 0.8) {
  ptitle <- names(dates)[i]
  trends <- trends[[i]]
  dates <- dates[[i]]$Year2
  signifs <- signifs[[i]]
  fit <- trends$fit
  uci <- fit + (1.96 * trends$se.fit)
  lci <- fit - (1.96 * trends$se.fit)
  mean_biomass <- obs[[i]]$mean_biomass_log
  incr <- unlist(signifs$incr)
  decr <- unlist(signifs$decr)
  incr.col <- "blue"
  decr.col <- "red"
  ylim <- range(mean_biomass, fit, uci, lci) 
  plot(dates, fit, ylim = ylim, type = "l",
       ylab="Log(Biomass (Kgs))", xlab="Years", main=ptitle)
  points(obs[[i]]$Year2, mean_biomass, pch = 16, cex = cex)
  lines(dates, uci, lty = "dashed")
  lines(dates, lci, lty = "dashed")
  lines(dates, incr, col = incr.col, lwd = 3)
  lines(dates, decr, col = decr.col, lwd = 3)
}

pdf("outputs/catches_lagoons_gam_fits_derivatives.pdf",
    height = 8, width = 12,
    pointsize = 12)

layout(1)
#par(mfrow=c(5,1),mar=c(0,0,1,0),oma=c(3,3,2,2))

par(mfrow=c(length(fits),1),mar=c(0,4,1,2),oma=c(3,3,2,2))
lapply(seq_along(fits[lagoons]), plotTrends, trends = predTrends,
       dates = predData, signifs = signifCores, obs = LagoonData[lagoons], cex = 1)
dev.off()

# For global model
lapply(seq_along(fits), plotTrends, trends = predTrends,
       dates = predData, signifs = signifCores, obs = full_model_data, cex = 1)



#--------------#

#extract dataframes from lists
df1 <- plyr::ldply(predData, data.frame) 
df2 <- plyr::ldply(predTrends, data.frame)[,-1] #this is to remove first column .id to avouid duplicates
df3 <- plyr::ldply(signifCores, data.frame) [,-1] 

gamderivAll <- cbind(df1, df2, df3) 

colnames(gamderivAll) <- c("lake", "Year", "fit", "se.fit", "incr", "decr") 

#Calculate standard error
gamderivAll <- mutate(gamderivAll, upper = fit + (2 * se.fit),
                      lower = fit - (2 * se.fit))

# this is to manually sort lakes
gamderivAll$lake = factor(gamderivAll$lake, levels=c("Encanyissada", "Tancada", "Canal vell", "Goleta"))

# Plot
deriv_gam_plot <- ggplot(gamderivAll) + 
  facet_grid(lake~., scales = "free_y") +
  geom_line(aes(x = Year, y = fit)) +
  geom_line(aes(x = Year, y = incr), color="blue", size=1.5) +
  geom_line(aes(x = Year, y = decr), color="red", size=1.5) +
  geom_ribbon(aes(x=Year, 
                  ymin = lower,
                  ymax = upper), alpha=0.1)+  
  geom_point(data=fishes_comm_lagoon, aes(x = Year2, y = mean_biomass_log), size=1,color="black") +
  labs(y = "Log(Total biomass (Kgs))", x = "Years") +
  ggtitle("Ebro Delta Coastal Lagoons Fish catches trend") +
  theme_bw()
deriv_gam_plot

# save plot
ggsave("outputs/gam_derivative_globalmodel_catches.png",
       plot = deriv_gam_plot,
       width=8,
       height=6,
       units="in",
       dpi = 400)



