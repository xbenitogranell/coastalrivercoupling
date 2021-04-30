## Model temporal trends of fish biomass

#Clear workspace
rm(list=ls(all=TRUE))

# load functions used
library(tidyverse)
library(ggplot2)
library(mgcv)
library(gratia)
library(car)

#Read in catches lagoons clean data
catches_clean <- read.csv("data/catches_clean.csv")[-1]
str(catches_clean)

#Species to be included in the models (>30% average biomass)
include <- c("Anguilla anguilla", "Cyprinus carpio", "Mullet ind.", "Sparus aurata",
             "Atherina boyeri", "Dicentrarchus labrax", "Liza ramada", "Sprattus sprattus")


## This is to fit a GAM on temporal total biomass
#Calculate average biomass per year and lagoon
fishes_comm <- catches_clean %>%
  filter(!is.na(Kgs)) %>% #remove NANs
  filter(!str_detect(Species, "Ind")) %>% #drop Undetermined species
  #filter(Species %in% include) %>% #select species to be modeled
  group_by(Year2) %>%
  summarise(mean_biomass = mean(Kgs,na.rm=T)) %>%
  mutate(mean_biomass_log = log(mean_biomass)) %>%
  mutate(year_f=factor(Year2)) %>%
  as.data.frame()
  

# Global model
set.seed(10)
model_gam <- gam(mean_biomass ~ s(Year2, k=40),
                 data=fishes_comm, family = Gamma(link = "log"),  
                 method = "REML", select=TRUE)
appraise(model_gam)

model_gam_log <- gam(mean_biomass_log ~ s(Year2, k=30),
                 data=fishes_comm, family = gaussian(link="identity"),  
                 method = "REML", select=TRUE)
appraise(model_gam_log)

#Compare different model fits using AIC
AIC(model_gam, model_gam_log)

model_gam <- model_gam_log

pacf(residuals(model_gam)) #check autocorrelation

#check gam outputs
gam.check(model_gam)
appraise(model_gam)
draw(model_gam)
summary(model_gam)


#controlling for temporal autocorrelation
library(lme4)
model_gamm <- gamm(mean_biomass_log ~ s(Year2, k=20),
                  data=fishes_comm, family = gaussian(link = "identity"),  
                  method = "REML", correlation=corCAR1(form = ~Year2))

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
  as.data.frame()


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
k <- c(rep(40,4))
fits <- lapply(seq_along(LagoonData[lagoons]), fitGam, data = LagoonData[lagoons], kk = k)
names(fits) <- names(LagoonData[lagoons])

#lapply(fits, summary)

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
  
  #this is for extract F value and F df
  Null.n <- as.numeric(model$df.null + 1)
  Fit.n <- as.numeric(model$df.null + 1)
  Null.res.df <- as.numeric(model$df.null)
  Fit.res.df <- as.numeric(model$df.residual) 
  Null.fac.df <- as.numeric(Null.n - Null.res.df)
  Fit.fac.df <- as.numeric(Fit.n - Fit.res.df)
  F.df1 <- as.numeric(Fit.fac.df - Null.fac.df)
  F.df2 <- as.numeric(Fit.n - Fit.fac.df)
  Null.MS <- as.numeric((null.deviance - fit.deviance) / (Fit.fac.df - Null.fac.df))
  Fit.MS <- as.numeric(fit.deviance / (Fit.n - Fit.fac.df))	  
  F.value <- as.numeric(Null.MS / Fit.MS)
  
  table <- cbind.data.frame(formula,n,family.function,family.link,fit.deviance,null.deviance,r.sq,
                            F.df1,F.value,s.table) #combine extracted columns
  colnames(table) <- c("formula","n","family.function", "family.link", "deviance", "null.deviance","r.sq",
  "s.edf", "s.ref.df", NA, "F.df1", "F.value", "p") #
  return(table)

}

summaryGam <- lapply(seq_along(fits), GAM_results, models=fits)
names(summaryGam) <- names(fits)

#extract dataframes from the summary GAM list and save table
summary_Gam_tables <- plyr::ldply(summaryGam, data.frame) 
write.table(summary_Gam_tables, file = "outputs/summary_gam_lagoons_table.txt", 
            sep = "\t", na = "", col.names = TRUE, row.names = FALSE)


#Load function to calculate first derivative on fish community biomass
source("scripts/Deriv.R")

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

deriv_gam_plot <- ggplot(gamderivAll) + 
  geom_line(aes(x = Year, y = fit)) +
  geom_line(aes(x = Year, y = incr), color="blue", size=1.5) +
  geom_line(aes(x = Year, y = decr), color="red", size=1.5) +
  geom_ribbon(aes(x=Year, 
                  ymin = lower,
                  ymax = upper), alpha=0.1)+  
  geom_point(data=fishes_comm, aes(x = Year2, y = mean_biomass_log), size=1,color="black") +
  labs(y = "Log(Total biomass (Kgs))", x = "Years") +
  ggtitle("Ebro Delta Coastal Lagoons Fish catches trend")
deriv_gam_plot

ggsave("outputs/gam_derivative_globalmodel_catches.png",
       plot = deriv_gam_plot,
       width=8,
       height=6,
       units="in",
       dpi = 400)

# this is to manually sort lakes
gamderivAll$Lake_f = factor(gamderivAll$lake, levels=c("Encanyissada", "Tancada", "Canal vell", "Goleta"))

deriv_plot <- ggplot(gamderivAll) +
  facet_grid(lake~., scales = "free_y") +
  geom_line(aes(x = Year, y = fit)) +
  geom_line(aes(x = Year, y = incr), color="blue", size=1.5) +
  geom_line(aes(x = Year, y = decr), color="red", size=1.5) +
  geom_ribbon(aes(x=Year, 
                  ymin = lower,
                  ymax = upper), alpha=0.1)+  
  geom_point(data=fishes_comm_lagoon, aes(x = Year2, y = mean_biomass_log), size=1,color="black") +
  labs(y = "Total biomass (Kgs)", x = "Years") 
deriv_plot


#ggsave("planktic_benthic_ratio_gam_plot.png", deriv_plot, height = 8, width = 10)




## This is to model temporal trends in biomass between spp 
#model S HGAM : similar smootheness between groups (spp) without global smooth 

#Species to be included in the models (>30% average biomass)
include <- c("Anguilla anguilla", "Cyprinus carpio", "Mullet ind.", "Sparus aurata",
             "Atherina boyeri", "Dicentrarchus labrax", "Liza ramada")

#Calculate average biomass per year and lagoon
fishes_spp <- catches_clean %>%
  filter(!is.na(Kgs)) %>% #remove NANs
  filter(Species %in% include) %>% #select species to be modeled
  group_by(Year2, Species) %>%
  summarise(mean_biomass = mean(Kgs,na.rm=T)) %>%
  mutate(mean_biomass_log = log(mean_biomass)) %>%
  mutate(year_f=factor(Year2)) %>%
  mutate(spp=factor(Species)) %>%
  as.data.frame()

levels(fishes_spp$spp)

set.seed(10) #set a seed so this is repeatable

model_gam_S <- gam(mean_biomass ~ s(Year2, spp, k=20, bs="fs"),
                    data=fishes_spp, family = Gamma(link = "log"), 
                    method = "REML", select = TRUE)

gam.check(model_gam_S)
draw(model_gam_S)
summary(model_gam_S)

#model I HGAM: different smootheness for each taxa without global smooth
model_gam_I<- gam(mean_biomass ~ s(Year2, by=spp, k=30, bs="fs") +
                     s(spp, bs="re"),
                   data=fishes_spp, family = Gamma(link = "log"),
                   method = "REML", select=TRUE)

gam.check(model_gam_I)
draw(model_gam_I)
summary(model_gam_I)


#Compare different model fits using AIC
AIC_table <- AIC(model_gam_S, model_gam_I)%>%
  rownames_to_column(var= "Model")%>%
  mutate(data_source = rep(c("diatom_data")))%>%
  group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  dplyr::select(-data_source)%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0)))

AIC_table


#Create synthetic data to predict over a range of ages
fishes_plot_data <- with(fishes_spp, as_tibble(expand.grid(Year2 = seq(min(fishes_spp$Year2), max(fishes_spp$Year2)),
                                                        spp = factor(levels(fishes_spp$spp)))))

fishes_modS_fit <- predict(model_gam_S, 
                         newdata = fishes_plot_data,
                         se.fit = TRUE)

fishes_modI_fit <- predict(model_gam_I,
                         newdata = fishes_plot_data,
                         se.fit = TRUE)

#non-shared trends
fishes_plot_data$modS_fit <- as.numeric(fishes_modS_fit$fit)
fishes_plot_data$modI_fit <- as.numeric(fishes_modI_fit$fit)

# comparing non-shared trends
fishes_plot_data <- gather(fishes_plot_data, key=model, value=fit, modS_fit, modI_fit)

fishes_plot_data <- mutate(fishes_plot_data, se= c(as.numeric(fishes_modS_fit$se.fit),
                                               as.numeric(fishes_modI_fit$se.fit)),
                         upper = exp(fit + (2 * se)),
                         lower = exp(fit - (2 * se)),
                         fit   = exp(fit))


#Plot the model output for non-shared trends, with means plus standard deviations for each model.
fishes_plot_model_labels <- paste("Model", c("S", "I"))
fishes_plot_model_labels <- factor(fishes_plot_model_labels, levels = fishes_plot_model_labels)

#non-shared trends
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

fishes_plot <- ggplot(fishes_plot_data) +
  facet_wrap(~spp, nrow = 4,scales = "free_y")+
  geom_ribbon(aes(x=Year2,
                  ymin = lower,
                  ymax = upper,
                  fill = model),
              alpha=0.2)+
  geom_point(data= fishes_spp, aes(x = Year2, y = mean_biomass), size=0.06) +
  #scale_y_continuous(trans = "log10")+
  geom_line(aes(x = Year2, y = fit, color = model))+
  labs(y = "Biomass (Kgs)", x = "Years") +
  scale_fill_brewer(name = "", palette = "Dark2",
                    labels = fishes_plot_model_labels) +
  scale_colour_brewer(name = "",
                      palette = "Dark2", labels = fishes_plot_model_labels)+
  theme(legend.position = "top",
        strip.text = element_text(size=10))

fishes_plot




#This chunk will compare individual species dynamics among lagoons
taxa <- fishes %>% filter(Species=="Anguilla anguilla") %>%
  mutate(lagoon_f=factor(Lagoon)) %>%
  filter(!lagoon_f %in% c("Clot - Baseta", "Platjola")) #few occurrences
  
taxa$lagoon_f <- droplevels(taxa$lagoon_f)                      

taxa_data <- taxa
str(taxa_data)

# Model G HGAM: global smoother
taxa_modG <- gam(Kgs ~ s(lagoon_f, Year2, bs="re") + s(lagoon_f, bs="re"),
                     data=taxa_data,
                     family=Gamma(link = "log"),
                     method="REML", select=TRUE,
                     drop.unused.levels = FALSE)


# Model GS: A global smoother plus group-level smoothers that have the same wiggliness
taxa_modGS <- gam(Kgs ~ s(Year2, lagoon_f, bs="re") +
                        s(Year2, k=10), data=taxa_data, 
                        family=Gamma(link = "log"), method="REML",select=TRUE,
                        drop.unused.levels = FALSE)


# Model GI: A global smoother plus group-level smoothers with differing wiggliness
taxa_modGI <- gam(Kgs ~ s(Year2, by=lagoon_f, k=10)+
                        s(lagoon_f, bs="re") + 
                        s(lagoon_f, Year2,bs="re"),
                      data=taxa_data,
                      family=Gamma(link = "log"),
                      method="REML",select=TRUE, drop.unused.levels = FALSE)

#qqplot, using gratia's qq_plot function, with simulated confidence intervals
pltG <- qq_plot(taxa_modG, method = "simulate")+
  labs(subtitle = NULL, title=NULL)
pltGS <- qq_plot(taxa_modGS, method = "simulate")+
  labs(subtitle = NULL, title=NULL, y=NULL)
pltGI <- qq_plot(taxa_modGI, method = "simulate")+
  labs(subtitle = NULL, title=NULL, y=NULL)

plot_grid(pltG, pltGS,pltGI, ncol = 3, align = "hv", axis = "lrtb",labels=c("G","GS","GI"))



#Compare different model fits using AIC
AIC_table <- AIC(taxa_modG, taxa_modGS, taxa_modGI)%>%
  rownames_to_column(var= "Model")%>%
  mutate(data_source = rep(c("taxa_data")))%>%
  group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  dplyr::select(-data_source)%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0)))

AIC_table


#Create synthetic data to use to compare predictions
taxa_plot_data <- expand.grid(lagoon_f = levels(taxa_data$lagoon_f),
                              Year2 = seq(min(taxa_data$Year2), max(taxa_data$Year2)))

#extract predicted values and standard errors for the three models
taxa_modG_fit <- predict(taxa_modG, 
                         newdata = taxa_plot_data, 
                         se.fit = TRUE)

taxa_modGS_fit <- predict(taxa_modGS, 
                          newdata = taxa_plot_data, 
                          se.fit = TRUE)


taxa_modGI_fit <- predict(taxa_modGI, 
                          newdata = taxa_plot_data, 
                          se.fit = TRUE)

taxa_plot_data$modG_fit <- as.numeric(taxa_modG_fit$fit)
taxa_plot_data$modGS_fit <- as.numeric(taxa_modGS_fit$fit)
taxa_plot_data$modGI_fit <- as.numeric(taxa_modGI_fit$fit)

taxa_plot_data <- gather(taxa_plot_data, 
                         key = model, 
                         value = fit, 
                         modG_fit, 
                         modGS_fit, 
                         modGI_fit)

taxa_plot_data <- mutate(taxa_plot_data, 
                         se = c(as.numeric(taxa_modG_fit$se.fit),
                                as.numeric(taxa_modGS_fit$se.fit),
                                as.numeric(taxa_modGI_fit$se.fit)),
                         upper = exp(fit + (2 * se)),
                         lower = exp(fit - (2 * se)),
                         fit   = exp(fit))

taxa_plot_model_labels = paste("Model", c("G","GS","GI"))
taxa_plot_model_labels = factor(taxa_plot_model_labels, 
                                levels= taxa_plot_model_labels)

# Plot

#define title
title <- taxa_data$Species[1]

taxa_plot <- ggplot(taxa_plot_data, aes(x=Year2))+
  facet_wrap(~lagoon_f, nrow = 2)+
  geom_ribbon(aes(x = Year2, ymin = lower, ymax = upper, fill = model), 
              alpha = 0.2) +
  geom_point(data= taxa_data, 
             aes(x = Year2, 
                 y = Kgs),
             size=0.06)+
  geom_line(aes(x = Year2, y = fit, colour = model)) +
  labs(y = "Biomass (Kgs)", 
       x = "Year") +
  scale_x_continuous(expand = c(0,0))+
  scale_fill_brewer(name = "", 
                    palette = "Dark2",
                    labels = taxa_plot_model_labels) +
  scale_colour_brewer(name = "",
                      palette = "Dark2", 
                      labels = taxa_plot_model_labels)+
  theme(legend.position = "top") +
  ggtitle(title)

taxa_plot

# save plot
# Save plot
ggsave("outputs/taxa_HGAM.png",
       plot = taxa_plot,
       width =10,
       height=8,
       units="in",
       dpi = 400)
