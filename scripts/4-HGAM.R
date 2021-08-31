## This is to model temporal trends in biomass between spp using Hierarchical GAMs

#Clear workspace
rm(list=ls(all=TRUE))

# load functions used
library(tidyverse)
library(ggplot2)
library(mgcv)
library(gratia)
library(car)
library(cowplot)

#Read in catches lagoons clean data
catches_clean <- read.csv("data/catches_clean.csv")[-1]
str(catches_clean)

#Species to be included in the models (>30% average biomass)
include <- c("Anguilla anguilla", "Cyprinus carpio", "Mullet ind.", "Sparus aurata",
             "Atherina boyeri", "Dicentrarchus labrax", "Liza ramada", "Sprattus sprattus")


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

levels(fishes_spp$spp) #check number of species (levels)


#model S HGAM: similar wiggliness between groups (spp) without a global smooth 
set.seed(10) #set a seed so this is repeatable

model_gam_S <- gam(mean_biomass ~ s(Year2, spp, k=30, bs="fs"),
                   data=fishes_spp, family = Gamma(link = "log"), 
                   method = "REML", select = TRUE)

gam.check(model_gam_S)
appraise(model_gam_S)
draw(model_gam_S)
summary(model_gam_S)

#model I HGAM: different wiggliness for each taxa without a global smooth
model_gam_I<- gam(mean_biomass ~ s(Year2, by=spp, k=30, bs="fs") +
                    s(spp, bs="re"),
                  data=fishes_spp, family = Gamma(link = "log"),
                  method = "REML", select=TRUE)

gam.check(model_gam_I)
appraise(model_gam_I)
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

### Save results
## Save fitted HGAM values

write.csv(fishes_plot_data, "outputs/HGAM_fishes_data.csv")

## Save plot
ggsave("outputs/taxa_HGAM.png",
       plot = taxa_plot,
       width =10,
       height=8,
       units="in",
       dpi = 400)

## This chunk will model individual species among lagoons
taxa <- "Anguilla anguilla" 

taxa_data <- catches_clean %>%
  filter(!is.na(Kgs)) %>% #remove NANs
  group_by(Year2, Species, Lagoon) %>%
  summarise(mean_biomass = mean(Kgs,na.rm=T)) %>%
  mutate(mean_biomass_log = log(mean_biomass)) %>%
  mutate(year_f=factor(Year2)) %>%
  mutate(spp=factor(Species)) %>%
  filter(Species==taxa) %>%
  mutate(lagoon_f=factor(Lagoon)) %>%
  filter(!lagoon_f %in% c("Clot - Baseta", "Platjola")) %>% #few occurrences
  as.data.frame()

taxa_data$lagoon_f <- droplevels(taxa_data$lagoon_f)                      
str(taxa_data)

# Model G HGAM: global smoother
taxa_modG <- gam(mean_biomass ~ s(lagoon_f, Year2, bs="re") + s(lagoon_f, bs="re"),
                 data=taxa_data,
                 family=Gamma(link = "log"),
                 method="REML", select=TRUE,
                 drop.unused.levels = FALSE)


# Model GS: A global smoother plus group-level smoothers that have the same wiggliness
taxa_modGS <- gam(mean_biomass ~ s(Year2, lagoon_f, bs="fs") +
                    s(Year2, k=10), data=taxa_data, 
                  family=Gamma(link = "log"), method="REML",select=TRUE,
                  drop.unused.levels = FALSE)


# Model GI: A global smoother plus group-level smoothers with differing wiggliness
taxa_modGI <- gam(mean_biomass ~ s(Year2, by=lagoon_f, k=10)+
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
                 y = mean_biomass),
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

# Save plot
ggsave("outputs/taxa_HGAM.png",
       plot = taxa_plot,
       width =10,
       height=8,
       units="in",
       dpi = 400)




