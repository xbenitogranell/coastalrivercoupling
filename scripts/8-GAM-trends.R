## Model river FQ and flow trends: data-driven GAM approach

library(mgcv)
library(tidyverse)
library(ggplot2)
library(gratia)

#https://github.com/SwampThingPaul/SWFL_RTHab/blob/main/src/LOSOM_CREHAB_v2.R   
#https://github.com/tetratech/baytrends   

# read in full data including river physico-chemical (from CHE, CAT), river flow (CHE), historical river flow (CHE), and mean fish catches (averaged across lagoons)
river_fisheries_data <- read.csv("outputs/river_fisheries_data.csv",row.names=1)
head(river_fisheries_data)
str(river_fisheries_data)

# plot the data
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


# Model mean biomass
m1 <- gam(mean_biomass_log ~ s(Year,bs="tp"),
          data=data_full, method = 'REML')

summary(m1)
appraise(m1) #check model residuals
draw(m1, residuals = TRUE) #plot partial effects 

# Model river flow with linear trend
m2a <- gam(log10(qmedmes) ~  Year + s(Year,bs="tp") + s(Month,bs="cc",k=12) + ti(Month,Year,bs=c("cc","tp")),
          data=data_full, method = 'REML')

summary(m2a)
appraise(m2) #check model residuals
draw(m2a, residuals = TRUE) #plot partial effects 

# Model river flow without linear trend
m2b <- gam(log10(flow_che) ~ s(Year,bs="tp") + s(Month,bs="cc",k=12) + ti(Month,Year,bs=c("cc","tp")),
           data=data_full, method = 'REML')

summary(m2b)
draw(m2b, residuals = TRUE) #plot partial effects 


#Compare different model fits using AIC
AIC_table <- AIC(m2a, m2b)%>%
  rownames_to_column(var= "Model")%>%
  mutate(data_source = rep(c("model_data")))%>%
  group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  dplyr::select(-data_source)%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0)))
AIC_table


# Model some river FQ variables
m3 <- gam(log10(Total_phosphorous) ~ s(Month,bs="cc",k=4) + s(Year,bs="tp") + ti(Month,Year,bs=c("cc","tp")),
          data=data_full, method = 'REML')

summary(m3)
appraise(m3) #check model residuals
draw(m3, residuals=TRUE)

