## this code is for running a GAM with an IT approach

#Clear workspace
rm(list=ls(all=TRUE))
#setwd("~/coastalrivercoupling")

## incorporate a dummy variable of presence of river regime shift or breakpoints of predictors

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

taxa <- "Anguilla|Dicentrarchus labrax|Mullet|Sparus"
taxa <- "Anguilla"

# This is for average variables across Year while holding lagoon for random effects
data_spp_env <- data_full %>% 
  filter(!Lagoon=="Platjola" & !Lagoon=="Clot - Baseta") %>% #filter out 
  droplevels() %>%
  filter(str_detect (Species, taxa)) %>% #select species to be modeled
  dplyr::select(Year, flow_che, SRP, Lagoon, Species, mean_biomass, mean_biomass_log, Chl.Total, TOC, NO3, NO2, NH4,
                qmedmes, Wind_direction, Wind_speed, Sea_level_pressure, Precipitation, Mean_temperature) %>%
  group_by(Year, Lagoon, Species) %>%
  summarise(across(everything(), mean, na.rm=TRUE)) %>%
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
resp.vars <- unique(as.character(data_spp_env$Species))
#resp.vars <- resp.vars[-9] 
setwd("~/coastalrivercoupling/outputs") #Set wd for example outputs - will differ on your computer

# specify predictors
cyclic.vars <- c("Year")
factor.vars <- c("Lagoon")
pred.vars <- c("Year", "flow_che", "Chl.Total")
response.var <- c("mean_biomass_log")

use.dat <- na.omit(data_spp_env[,c(factor.vars,response.var,pred.vars, "Species")]) #get rid off NAs
out.all <- list()
var.imp <- list()
fss.all <- list()
top.all <- list()
i=1

# Loop through the FSS function for each Taxa----
for(i in 1:length(resp.vars)){
  
  use.dat <- use.dat[which(use.dat$Species==resp.vars[i]),]
  
  Model1 <- gam(mean_biomass_log ~ s(Year, k=10, bs="cc"),
             data=use.dat, method = "REML")
  
  model.set <- generate.model.set(use.dat=use.dat,
                               test.fit=Model1,
                               factor.smooth.interactions = TRUE,
                               cyclic.vars=cyclic.vars,
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

