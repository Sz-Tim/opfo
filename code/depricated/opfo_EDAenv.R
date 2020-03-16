# this script explores the environmental data with a main goal of error checking

library(tidyverse); source("code/lc_cols.R")
theme_set(theme_bw() + theme(panel.grid.minor=element_blank()))

write_data <- FALSE

plot_data <- readxl::read_xlsx("data/opfo_2019_envData.xlsx", 1) %>%
  mutate(BDM=as.numeric(BDM), 
         Sample_date=lubridate::as_date(Date, tz="GMT", format="%Y%m%d"))
site_data <- read_csv("~/Documents/unil/gis/data/site_info.csv")


summarise_plots <- function(x.df) {
  x.df %>% 
    summarise(n_plots=n(),
              SoilTemp_mean=mean(SoilTemp, na.rm=T),
              SoilTemp_sd=sd(SoilTemp, na.rm=T),
              Grass_mean=mean(Grass, na.rm=T),
              Grass_sd=sd(Grass, na.rm=T),
              Forb_mean=mean(Forb, na.rm=T),
              Forb_sd=sd(Forb, na.rm=T),
              Shrub_mean=mean(Shrub, na.rm=T),
              Shrub_sd=sd(Shrub, na.rm=T),
              Bare_mean=mean(Bare, na.rm=T),
              Bare_sd=sd(Bare, na.rm=T),
              Litter_mean=mean(Litter, na.rm=T),
              Litter_sd=sd(Litter, na.rm=T),
              Moss_mean=mean(Moss, na.rm=T),
              Moss_sd=sd(Moss, na.rm=T))
}

env_data <- full_join(select(site_data, -c(BDM_id, source)), 
                      plot_data, 
                      by="BDM")
site_summary <- plot_data %>% group_by(BDM_id, BDM, Sample_date) %>%
  summarise_plots %>%
  full_join(select(site_data, -c(BDM_id, source)), ., by="BDM")
lc_summary <- plot_data %>% group_by(BDM_id, BDM, Sample_date, Categorie) %>%
  summarise_plots %>%
  full_join(select(site_data, -c(BDM_id, source)), ., by="BDM")

if(write_data) {
  write_csv(env_data, "data/opfo_envDataProcessed.csv")
  write_csv(site_summary, "data/opfo_siteSummaryProcessed.csv")
  write_csv(lc_summary, "data/opfo_lcSummaryProcessed.csv")
}


########
## general checks
########
str(env_data)
summary(env_data)



########
## graphical checks
########

# Soil temperature
ggplot(env_data, aes(x=MAT_mean/10, y=SoilTemp, colour=el_mean)) + 
  stat_smooth(method="loess", span=2) + geom_point() + 
  labs(x="Mean Annual Temp (ºC)", y="Soil Temp (ºC)") +
  scale_color_viridis_c("Elevation (m)")
ggplot(env_data, aes(x=MAT_mean/10, y=SoilTemp, colour=Sample_date)) + 
  stat_smooth(method="loess", span=2) + geom_point() + 
  labs(x="Mean Annual Temp (ºC)", y="Soil Temp (ºC)") 
ggplot(env_data, aes(x=Sample_date, y=el_mean, colour=SoilTemp)) + 
  stat_smooth(method="loess", span=2) + geom_point() + 
  labs(x="Sample Date", y="Elevation") 
ggplot(env_data, aes(x=el_mean, y=SoilTemp, colour=Sample_date)) + 
  geom_point() + labs(x="Elevation", y="Soil Temp")
ggplot(env_data, aes(x=Sample_date, y=SoilTemp, colour=el_mean)) + 
  geom_point() + labs(x="Date", y="Soil Temp")

ggplot(env_data, aes(x=el_mean, y=SoilTemp, colour=Categorie)) + 
  geom_jitter(alpha=0.7) + labs(x="Elevation", y="Soil Temp") +
  scale_colour_manual(values=lc_cols) + 
  stat_smooth(se=F, span=2) 
ggplot(env_data, aes(x=el_mean, y=Grass, colour=Categorie)) + 
  geom_jitter(alpha=0.7) + labs(x="Elevation", y="Grass") +
  scale_colour_manual(values=lc_cols) + 
  stat_smooth(se=F, span=2) + 
  facet_wrap(~Categorie)
ggplot(env_data, aes(x=el_mean, y=Forb, colour=Categorie)) + 
  geom_jitter(alpha=0.7) + labs(x="Elevation", y="Forb") +
  scale_colour_manual(values=lc_cols) + 
  stat_smooth(se=F, span=2) + 
  facet_wrap(~Categorie)
ggplot(env_data, aes(x=el_mean, y=Shrub, colour=Categorie)) + 
  geom_jitter(alpha=0.7) + labs(x="Elevation", y="Shrub") +
  scale_colour_manual(values=lc_cols) + 
  stat_smooth(se=F, span=2) + 
  facet_wrap(~Categorie)
ggplot(env_data, aes(x=el_mean, y=Bare, colour=Categorie)) + 
  geom_jitter(alpha=0.7) + labs(x="Elevation", y="Bare") +
  scale_colour_manual(values=lc_cols) + 
  stat_smooth(se=F, span=2) + 
  facet_wrap(~Categorie)
ggplot(env_data, aes(x=el_mean, y=Litter, colour=Categorie)) + 
  geom_jitter(alpha=0.7) + labs(x="Elevation", y="Litter") +
  scale_colour_manual(values=lc_cols) + 
  stat_smooth(se=F, span=2) + 
  facet_wrap(~Categorie)
ggplot(env_data, aes(x=el_mean, y=Moss, colour=Categorie)) + 
  geom_jitter(alpha=0.7) + labs(x="Elevation", y="Moss") +
  scale_colour_manual(values=lc_cols) + 
  stat_smooth(se=F, span=2) + 
  facet_wrap(~Categorie)


# Mean and sd by elevation
ggplot(site_summary, aes(x=el_mean, y=SoilTemp_mean, colour=Sample_date)) +
  geom_pointrange(aes(ymin=SoilTemp_mean-SoilTemp_sd, 
                      ymax=SoilTemp_mean+SoilTemp_sd)) +
  labs(x="Elevation", y="Soil Temp (mean ± SD)")
ggplot(site_summary, aes(x=el_mean, y=SoilTemp_sd, colour=Sample_date)) +
  geom_point(aes()) +
  labs(x="Elevation", y="Soil Temp SD")

ggplot(site_summary, aes(x=el_mean, y=Grass_mean, colour=Sample_date)) +
  geom_pointrange(aes(ymin=Grass_mean-Grass_sd, 
                      ymax=Grass_mean+Grass_sd)) +
  labs(x="Elevation", y="Grass (mean ± SD)")
ggplot(site_summary, aes(x=el_mean, y=Forb_mean, colour=Sample_date)) +
  geom_pointrange(aes(ymin=Forb_mean-Forb_sd, 
                      ymax=Forb_mean+Forb_sd)) +
  labs(x="Elevation", y="Forb (mean ± SD)")
ggplot(site_summary, aes(x=el_mean, y=Shrub_mean, colour=Sample_date)) +
  geom_pointrange(aes(ymin=Shrub_mean-Shrub_sd, 
                      ymax=Shrub_mean+Shrub_sd)) +
  labs(x="Elevation", y="Shrub (mean ± SD)")
ggplot(site_summary, aes(x=el_mean, y=Bare_mean, colour=Sample_date)) +
  geom_pointrange(aes(ymin=Bare_mean-Bare_sd, 
                      ymax=Bare_mean+Bare_sd)) +
  labs(x="Elevation", y="Bare (mean ± SD)")
ggplot(site_summary, aes(x=el_mean, y=Litter_mean, colour=Sample_date)) +
  geom_pointrange(aes(ymin=Litter_mean-Litter_sd, 
                      ymax=Litter_mean+Litter_sd)) +
  labs(x="Elevation", y="Litter (mean ± SD)")
ggplot(site_summary, aes(x=el_mean, y=Moss_mean, colour=Sample_date)) +
  geom_pointrange(aes(ymin=Moss_mean-Moss_sd, 
                      ymax=Moss_mean+Moss_sd)) +
  labs(x="Elevation", y="Moss (mean ± SD)")


# by categorie
ggplot(lc_summary, aes(x=el_mean, y=SoilTemp_mean, colour=Categorie)) +
  geom_point() + facet_wrap(~Categorie) +
  geom_linerange(aes(ymin=SoilTemp_mean-SoilTemp_sd, 
                      ymax=SoilTemp_mean+SoilTemp_sd)) +
  scale_colour_manual(values=lc_cols) + 
  labs(x="Elevation", y="Soil Temp (mean ± SD)")
ggplot(lc_summary, aes(x=el_mean, y=SoilTemp_sd, colour=Categorie)) +
  geom_point(aes()) + facet_wrap(~Categorie) +
  scale_colour_manual(values=lc_cols) + 
  labs(x="Elevation", y="Soil Temp SD")

ggplot(lc_summary, aes(x=el_mean, y=Grass_mean, colour=Categorie)) +
  geom_point() + facet_wrap(~Categorie) +
  geom_linerange(aes(ymin=Grass_mean-Grass_sd, 
                      ymax=Grass_mean+Grass_sd)) +
  scale_colour_manual(values=lc_cols) +
  labs(x="Elevation", y="Grass (mean ± SD)")
ggplot(lc_summary, aes(x=el_mean, y=Forb_mean, colour=Categorie)) +
  geom_point() + facet_wrap(~Categorie) +
  geom_linerange(aes(ymin=Forb_mean-Forb_sd, 
                      ymax=Forb_mean+Forb_sd)) +
  scale_colour_manual(values=lc_cols) +
  labs(x="Elevation", y="Forb (mean ± SD)")
ggplot(lc_summary, aes(x=el_mean, y=Shrub_mean, colour=Categorie)) +
  geom_point() + facet_wrap(~Categorie) +
  geom_linerange(aes(ymin=Shrub_mean-Shrub_sd, 
                      ymax=Shrub_mean+Shrub_sd)) +
  scale_colour_manual(values=lc_cols) +
  labs(x="Elevation", y="Shrub (mean ± SD)")
ggplot(lc_summary, aes(x=el_mean, y=Bare_mean, colour=Categorie)) +
  geom_point() + facet_wrap(~Categorie) +
  geom_linerange(aes(ymin=Bare_mean-Bare_sd, 
                      ymax=Bare_mean+Bare_sd)) +
  scale_colour_manual(values=lc_cols) +
  labs(x="Elevation", y="Bare (mean ± SD)")
ggplot(lc_summary, aes(x=el_mean, y=Litter_mean, colour=Categorie)) +
  geom_point() + facet_wrap(~Categorie) +
  geom_linerange(aes(ymin=Litter_mean-Litter_sd, 
                      ymax=Litter_mean+Litter_sd)) +
  scale_colour_manual(values=lc_cols) +
  labs(x="Elevation", y="Litter (mean ± SD)")
ggplot(lc_summary, aes(x=el_mean, y=Moss_mean, colour=Categorie)) +
  geom_point() + facet_wrap(~Categorie) +
  geom_linerange(aes(ymin=Moss_mean-Moss_sd, 
                      ymax=Moss_mean+Moss_sd)) +
  scale_colour_manual(values=lc_cols) +
  labs(x="Elevation", y="Moss (mean ± SD)")





