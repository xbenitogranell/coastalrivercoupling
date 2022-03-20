## Model river FQ and flow trends: data-driven GAM approach

library(mgcv)
library(tidyverse)
library(ggplot2)
library(gratia) # for visualizing GAM model outputs

#https://github.com/SwampThingPaul/SWFL_RTHab/blob/main/src/LOSOM_CREHAB_v2.R   
#https://github.com/tetratech/baytrends   

# read in full data including river physico-chemical (from CHE, CAT), river flow (CHE), historical river flow (CHE), and mean fish catches (averaged across lagoons)
river_fisheries_data <- read.csv("outputs/river_fisheries_data.csv",row.names=1)
head(river_fisheries_data)
str(river_fisheries_data)

# Model mean biomass
m1 <- gam(mean_biomass_log ~ s(Year,bs="tp"),
          data=river_fisheries_data, method = 'REML')

summary(m1)
appraise(m1) #check model residuals
draw(m1, residuals = TRUE) #plot partial effects 

# Model river flow with linear trend
m2a <- gam(log10(qmedmes) ~  Year + s(Year,bs="tp") + s(Month,bs="cc",k=12) + ti(Month,Year,bs=c("cc","tp")),
           data=data_full, method = 'REML')

summary(m2a)
appraise(m2a) #check model residuals
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

## GAM Covariates model: test for temporal contributions of each covariate on fish catch trends over time
# red in full data including river physico-chemical (from CHE, CAT), historical river flow (CHE), species biomass data per lagoon, and climatic data
data_full <- read.csv("outputs/river_fisheries_lagoon_spp_data.csv",row.names=1)
head(data_full)

# Create lagged variables
# bellver$y.1 <- c(NA,bellver$y[1:(nrow(bellver)-1)]) #y is the variable to be lagged

# check how many levels per species and lagoobn
levels(data_full$Species)
levels(data_full$Lagoon)

taxa <- "Anguilla"

# This is for average variables across Year while holding lagoon for random effects
data_spp_env <- data_full %>% 
  filter(!Lagoon=="Platjola" & !Lagoon=="Clot - Baseta") %>% #filter out 
  droplevels() %>%
  filter(str_detect (Species, taxa)) %>% #select species to be modeled
  dplyr::select(Year, flow_che, SRP, Lagoon, mean_biomass, mean_biomass_log, Chl.Total, TOC, NO3, NO2, NH4,
         qmedmes, Wind_direction, Wind_speed, Sea_level_pressure, Precipitation, Mean_temperature) %>%
  group_by(Year, Lagoon) %>%
  summarise(across(everything(), mean, na.rm=TRUE)) 
  
# plot the data to check everything looks correct
ggplot(data=gather(data_spp, variable, value, -Year, -Lagoon),
       aes(x=Year,
           y=value,
           group=variable)) +
  geom_line() +
  geom_smooth() +
  facet_wrap("variable",scales = "free_y", ncol=2) +
  xlab("Years") +
  ylab("") +
  ggtitle("") +
  theme_bw()

# This chunk run a GAM model      
str(data_spp_env)     

mod1 <- gam(mean_biomass ~ s(Year, k=20) + s(flow_che, k=15, bs="ad") + s(Lagoon, bs="re"),
            data = data_spp_env, method = "REML", 
            select = TRUE, family = gaussian(link="log"),
            na.action = na.omit) 

pacf(residuals(mod1)) # indicates non AR1
plot(mod1, page=1, scale=0)
gam.check(mod1)
appraise(mod1)
draw(mod1)
summary(mod1)

#accounting for temporal autocorrelation
mod1.car <- gam(mean_biomass ~ s(Year, k=20) + s(flow_che, k=15, bs="ad") + s(Lagoon, bs="re"),
            data = data_spp_env, method = "REML", 
            select = TRUE, family = gaussian(link="log"),
            na.action = na.omit,
            correlation = corCAR1(form = ~ Year)) 

pacf(residuals(mod1.car$lme)) # indicates AR1
plot(mod1.car$gam, page=1, scale = 0)
summary(mod1.car$gam)
draw(mod1.car)

#Compare different model fits using AIC
AIC_table <- AIC(mod1, mod1.car$lme)%>%
  rownames_to_column(var= "Model")%>%
  mutate(data_source = rep(c("diatom_data")))%>%
  group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  dplyr::select(-data_source)%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0)))

AIC_table

# Predict GAM
## data to predict at
#Check NAs and remove rows
row.has.na <- apply(data_spp_env, 1, function(x){any(is.na(x))})
sum(row.has.na)
data_spp_env <- data_spp_env[!row.has.na,]

## Predict
predGam <- cbind(data_spp_env, 
                 data.frame(predict.gam(mod1, data_spp_env, 
                                        type = "terms" , se.fit = TRUE)))
#plot
var <- predGam$
se.var <- predGam$

predGamPlt <- ggplot(predGam, aes(x = Age, y = var)) +
  geom_line() +
  geom_ribbon(aes(ymin = var + (2 * se.var), ymax = var - (2 * se.var)),alpha=0.4) +
  geom_point()+
  #scale_x_reverse() +
  labs(y = "Agropastoralism", x = "cal years BP", title = "")+
  theme(legend.position = "none")+
  geom_hline(yintercept=0, linetype="dashed")+
  theme_bw()+ theme(legend.position = "none")+
  #ggtitle(expression("Agropastoralism" %->% "Diatoms")) +
  theme(strip.text.x = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=12),
        axis.title.y=element_text(size=14))
