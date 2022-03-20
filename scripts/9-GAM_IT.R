## this code is for running a GAM with an IT approach

#Clear workspace
rm(list=ls(all=TRUE))
setwd("~/coastalrivercoupling")

## incorporate a dummy variable of presence of river regime shift or breakpoints of predictors

# load functions used
library(FSSgam)
library(mgcv)
library(tidyverse)

# read in full data including river physico-chemical (from CHE, CAT), historical river flow (CHE),  fish catches (year, lagoon) and climatic data
dat_full <- read.csv("outputs/river_fisheries_lagoon_spp_data.csv",row.names=1)
head(dat_full)
str(dat_full)

pred.vars <- colnames(dat_full[,c(4,5,6,7,9,13,14)])

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
resp.vars <- unique(as.character(dat_full$Species))
resp.vars <- resp.vars[-9] 
setwd("~/coastalrivercoupling/outputs") #Set wd for example outputs - will differ on your computer

# specify predictors
cyclic.vars <- c("Year")
factor.vars <- c("Lagoon")
pred.vars <- c("Year", "flow_che")
response.var <- c("mean_biomass_log")


use.dat <- na.omit(dat_full[,c(factor.vars,response.var,pred.vars)]) #get rid off NAs
out.all <- list()
var.imp <- list()

# Loop through the FSS function for each Taxa----
for(i in 1:length(resp.vars)){
  use.dat <- use.dat[which(use.dat$Species==resp.vars[i]),]
  
  Model1 <- gam(mean_biomass_log ~ s(Year) + s(Lagoon,bs="re"),
             data=use.dat)
  
  model.set <- generate.model.set(use.dat=use.dat,
                               test.fit=Model1,
                               cyclic.vars=cyclic.vars,
                               pred.vars.cont=pred.vars,
                               k=3)
  
  out.list=fit.model.set(model.set,
                         max.models=600,
                         parallel=T)
  names(out.list)
  
  out.list$failed.models # examine the list of failed models
  mod.table=out.list$mod.data.out  # look at the model selection table
  mod.table=mod.table[order(mod.table$AICc),]
  mod.table$cumsum.wi=cumsum(mod.table$wi.AICc)
  out.i=mod.table[which(mod.table$delta.AICc<=3),]
  out.all=c(out.all,list(out.i))
  # var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Either raw importance score
  var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Or importance score weighted by r2
  
  # plot the best models
  for(m in 1:nrow(out.i)){
    best.model.name=as.character(out.i$modname[m])
    
    png(file=paste(name,m,resp.vars[i],"mod_fits.png",sep="_"))
    if(best.model.name!="null"){
      par(mfrow=c(3,1),mar=c(9,4,3,1))
      best.model=out.list$success.models[[best.model.name]]
      plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
      mtext(side=2,text=resp.vars[i],outer=F)}  
    dev.off()
  }
}

# Model fits and importance---
names(out.all)=resp.vars
names(var.imp)=resp.vars
all.mod.fits=do.call("rbind",out.all)
all.var.imp=do.call("rbind",var.imp)
write.csv(all.mod.fits[,-2],file=paste(name,"all.mod.fits.csv",sep="_"))
write.csv(all.var.imp,file=paste(name,"all.var.imp.csv",sep="_"))

# Generic importance plots-
heatmap.2(all.var.imp,notecex=0.4,  dendrogram ="none",
          col=colorRampPalette(c("white","yellow","red"))(10),
          trace="none",key.title = "",keysize=2,
          notecol="black",key=T,
          sepcolor = "black",margins=c(12,8), lhei=c(4,15),Rowv=FALSE,Colv=FALSE)

