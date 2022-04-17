## this code is for running a GAM with an IT approach

#Clear workspace
rm(list=ls(all=TRUE))
#setwd("~/coastalrivercoupling")

# load functions used
library(FSSgam)
library(mgcv)
library(tidyverse)
library(MuMIn)

# read in full data including river physico-chemical (from CHE, CAT), historical river flow (CHE),  fish catches (year, lagoon) and climatic data
data_full <- read.csv("outputs/river_fisheries_lagoon_spp_data.csv",row.names=1)
head(data_full)
str(data_full)
data_full$Species

# read in breakpoints river environment
breakpoints_river <- read.csv("outputs/breakpoints_river.csv", row.names = 1)

# taxa <- "Anguilla|Dicentrarchus labrax|Mullet|Sparus"
# taxa <- "Anguilla"

# This is for average variables across Year while holding lagoon for random effects
data_spp_env <- data_full %>% 
  filter(!Lagoon=="Platjola" & !Lagoon=="Clot - Baseta") %>% #filter out 
  droplevels() %>%
  group_by(Year, Lagoon) %>%
  dplyr::select(Year, flow_che, SRP_che, mean_biomass, mean_biomass_log, Chl.Total, TOC, NO3, NO2, NH4,
                qmedmes, Wind_direction, Wind_speed, Sea_level_pressure, Precipitation, Mean_temperature, QMax, QMean) %>%
  summarise(across(everything(), mean, na.rm=TRUE)) %>%
  mutate(qmedmes_lag1=lag(qmedmes)) %>% #create a lagged variable for historical flow
  left_join(breakpoints_river, by="Year") %>% #join with river environmental dataset of breakpoints
  as.data.frame()

names(data_spp_env)

#pred.vars <- colnames(data_spp_env[,c(4,5,6,7,9,13,14)])

# Plot of possible transformations
par(mfrow=c(3,2))
for (i in pred.vars) {
  x <- dat_full[ ,i]
  x = as.numeric(unlist(x))
  hist((x)) 
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

# Run the full subset model selection----
loop.vars <- unique(as.character(data_spp_env$Species)) #loop through species
loop.vars <- unique(as.character(data_spp_env$Lagoon)) #loop through lagoons

#resp.vars <- resp.vars[-9] 
setwd("~/coastalrivercoupling/outputs") #Set wd for example outputs - will differ on your computer

# specify predictors
#cyclic.vars <- c("Year")
factor.vars <- c("Lagoon")
pred.vars <- c("Year", "flow_che", "Chl.Total")
response.var <- c("mean_biomass_log")

use.dat <- na.omit(data_spp_env[,c(factor.vars,response.var,pred.vars, "Species")]) #get rid off NAs
use.dat <- na.omit(data_spp_env[,c(factor.vars,response.var,pred.vars, "Lagoon")]) #get rid off NAs

out.all <- list()
var.imp <- list()
fss.all <- list()
top.all <- list()
i=1

# Loop through the FSS function for each Taxa----
for(i in 1:length(loop.vars)){
  use.dat <- use.dat[which(use.dat$Species==resp.vars[i]),] #loop through species
  use.dat <- use.dat[which(use.dat$Lagoon==lagoon.vars[i]),] #loop through lagoons
  
  Model1 <- gam(mean_biomass_log ~ s(Year, k=10, bs="cc"),
             data=use.dat, method = "REML")
  model.set <- generate.model.set(use.dat=use.dat,
                               test.fit=Model1,
                               pred.vars.cont=pred.vars, 
                               max.predictors = 3)
  out.list=fit.model.set(model.set)
  #names(out.list)
  # examine the list of failed models
  #out.list$failed.models
  #out.list$success.models
  fss.all=c(fss.all,list(out.list))
  mod.table=out.list$mod.data.out
  mod.table=mod.table[order(mod.table$AICc),]
  out.i=mod.table
  out.all=c(out.all,list(out.i))
  var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw))
  all.less.2AICc=mod.table[which(mod.table$delta.AICc<2),]
  top.all=c(top.all,list(all.less.2AICc))
  
  # plot the all best models
  par(oma=c(1,1,4,1))
  for(r in 1:nrow(all.less.2AICc)){
    best.model.name=as.character(all.less.2AICc$modname[r])
    best.model=out.list$success.models[[best.model.name]]
    if(best.model.name!="null"){
      plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
      mtext(side=3,text=resp.vars[i],outer=T)}
  }
}
dev.off()
  
names(out.all)=resp.vars
names(var.imp)=resp.vars
names(top.all)=resp.vars
names(fss.all)=resp.vars

all.mod.fits=do.call("rbind",out.all)
all.var.imp=do.call("rbind",var.imp)
top.mod.fits=do.call("rbind",top.all)

require(car)
require(doBy)
require(gplots)
require(RColorBrewer)

pdf(file="var_importance_heatmap_functional_biomass.pdf",height=5,width=7,pointsize=10)

heatmap.2(all.var.imp,notecex=0.4,  dendrogram ="none",
          col=colorRampPalette(c("white","yellow","orange","red"))(30),
          trace="none",key.title = "",keysize=2,
          notecol="black",key=T,
          sepcolor = "black",margins=c(12,14), lhei=c(3,10),lwid=c(3,10),
          Rowv=FALSE,Colv=FALSE)
dev.off()

# write.csv(all.mod.fits[,-2],"all_model_fits_functional_biomass.csv")
# write.csv(top.mod.fits[,-2],"top_model_fits_functional_biomass.csv")
# write.csv(model.set$predictor.correlations,"predictor_correlations.csv")


## Individual models for each lagoon
use.dat <- data_spp_env 

# specify predictors
pred.vars <- c("Year", "flow_che", "Chl.Total", "qmedmes", "SRP_che", "TOC", "NO3", "Wind_speed", "Sea_level_pressure",
               "QMax", "qmedmes_lag1", 'chla_breakpoinnt', 'flow_breakpoinnt', 'PO4_breakpoinnt')
response.var <- c("mean_biomass_log")

use.dat <- na.omit(data_spp_env[,c(response.var,pred.vars, "Lagoon")]) #get rid off NAs

start.fit <- gam(mean_biomass_log ~ s(Year, bs="cc", k=10) +
                s(Lagoon, bs="re"),
              family= gaussian(link = "identity"),
              data=use.dat, method="REML", select=TRUE)

model.set <- generate.model.set(use.dat=use.dat,
                             test.fit=start.fit,
                             pred.vars.cont=pred.vars,
                             max.predictors=3)

out.list <- fit.model.set(model.set,parallel=T)
names(out.list)

#write.csv(out.list$predictor.correlations,"predictor_correlations.csv")
out.list$predictor.correlations

# examine the list of failed models
length(out.list$failed.models)
length(out.list$success.models)

# look at the model selection table
mod.table=out.list$mod.data.out
mod.table=mod.table[order(mod.table$AICc),]
head(mod.table)
#write.csv(mod.table[,-2],"modfits.csv")

barplot(out.list$variable.importance$bic$variable.weights.raw,las=2,
        ylab="Relative variable importance")

write.csv(model.set$predictor.correlations,"predictor_correlations.csv")

# extract the best model
mod.table=mod.table[order(mod.table$AIC),]
head(mod.table)

best.model=out.list$success.models[[as.character(mod.table$modname[1])]]
plot(best.model,all.terms=T,pages=1)

gam.check(best.model)
draw(best.model, residuals=TRUE)
summary(best.model)

