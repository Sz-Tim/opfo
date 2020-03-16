# This script is primarily for detecting abnormalities for quality control in
# the opération fourmis structured sampling ant data. Specifically, it checks:
#
# Quality checks:
# 1. General location of points
#   - are all points within Vaud?
#   - are all 999xxxx points within sampling squares?
# 2. Point type distribution
#   - do the number of soil/transect/tree points align with Tanja's datasheet?
# 3. Point-habitat alignment
#   - does the GIS-extracted habitat align with the label entered with the tube?
#
#
# General patterns:
# 1. Site-level -- colony density (soil) vs.
#   - elevation
#   - MAT
#   - TAR
#   - soil temperature (mean)
#   - habitat proportions
# 2. Site-level -- colony density (transect) vs.
#   - elevation
#   - MAT
#   - TAR
#   - soil temperature (mean)
#   - habitat proportions








########
## set up
########

#--- libraries and options
library(tidyverse); library(ggspatial); library(sf); library(googlesheets)
library(viridis); theme_set(theme_bw()); source("code/lc_cols.R")
fonts <- theme(axis.text=element_text(size=14),
               axis.title=element_text(size=16),
               legend.text=element_text(size=12),
               legend.title=element_text(size=14),
               strip.text=element_text(size=14))
scales::trans_new("log_p1", 
                  function(x) {log(x+1)},
                  function(x) {exp(x)-1})
overwrite <- FALSE  # just a failsafe



#--- site, land cover, & plot info
site_summary <- read_csv("data/opfo_siteSummaryProcessed.csv")
lc_summary <- read_csv("data/opfo_lcSummaryProcessed.csv")
plot_i <- read_csv("data/opfo_envDataProcessed.csv")



#--- spatial data

## environmental
gis.dir <- "~/Documents/unil/GIS/data/"
VD_raw <- st_read(paste0(gis.dir, "Vaud_boundaries.shp")) 
VD <- st_union(VD_raw)
dem <- raster::raster(paste0(gis.dir, "aster/CH_21781.tif"))
slope <- raster::raster(paste0(gis.dir, "aster/CH_slope_21781.tif"))
dem_VD <- raster::crop(dem, VD_raw) %>% raster::mask(., st_zm(VD_raw))
dem_agg <- raster::aggregate(dem_VD, fact=42.5) # for ~1km cells

## site and plot
sites.sf <- st_read(paste0(gis.dir, "Fourmis_data_shp/MilieuxBDMdiss.shp")) %>%
  st_transform(3857)
site.box.sf <- st_read(paste0(gis.dir, "Fourmis_data_shp/BDMplus.shp")) %>%
  st_transform(st_crs(sites.sf)) %>% mutate(BDM=Id)
plot.sf <- st_read("data/gis/opfo_soil_25per.shp") %>% 
  st_transform(st_crs(sites.sf)) %>%
  mutate(Plot_id=str_remove_all(plot_id, "[[[:space:]]_]")) %>%
  select(Plot_id, propArea) %>% left_join(., plot_i, by="Plot_id") %>% 
  group_by(BDM_id) %>% mutate(nPlots=n()) %>% ungroup %>%
  mutate(el_plot=raster::extract(dem, .),
         slope_plot=raster::extract(slope, .))


#--- ant data
ant.df <- gs_read(ss=gs_title("Scientific Ant Inventory VD"), ws="Transfer") %>%
  filter(!is.na(lon)) %>% rename(TubeNo=CATALOGUENUMBER, Plot_id=PLOT) %>%
  mutate(SampleDate=lubridate::dmy(SampleDate)) %>%
  left_join(., select(plot.sf, -c(3,21:30)), by="Plot_id")


# convert veg metrics to midpoints
veg_lookup <- data.frame(score=c(0, 0.1, 1, 2, 3, 4, 5),
                         mp=c(0, 0.5, 2.5, 15, 37.5, 62.5, 87.5))
ant.df$Grass <- veg_lookup$mp[match(ant.df$Grass, veg_lookup$score)]
ant.df$Forb <- veg_lookup$mp[match(ant.df$Forb, veg_lookup$score)]
ant.df$Shrub <- veg_lookup$mp[match(ant.df$Shrub, veg_lookup$score)]
ant.df$Litter <- veg_lookup$mp[match(ant.df$Litter, veg_lookup$score)]
ant.df$Bare <- veg_lookup$mp[match(ant.df$Bare, veg_lookup$score)]
ant.df$Moss <- veg_lookup$mp[match(ant.df$Moss, veg_lookup$score)]
plot.sf$Grass <- veg_lookup$mp[match(plot.sf$Grass, veg_lookup$score)]
plot.sf$Forb <- veg_lookup$mp[match(plot.sf$Forb, veg_lookup$score)]
plot.sf$Shrub <- veg_lookup$mp[match(plot.sf$Shrub, veg_lookup$score)]
plot.sf$Litter <- veg_lookup$mp[match(plot.sf$Litter, veg_lookup$score)]
plot.sf$Bare <- veg_lookup$mp[match(plot.sf$Bare, veg_lookup$score)]
plot.sf$Moss <- veg_lookup$mp[match(plot.sf$Moss, veg_lookup$score)]

ant.df <- ant.df %>% group_by(BDM) %>% 
  mutate(elAnomaly=mnt25-mean(mnt25, na.rm=T))

ggplot(ant.df, aes(SampleDate, elAnomaly)) + geom_point(alpha=0.7)
ggplot(ant.df, aes(elAnomaly, SoilTempAnomaly)) + geom_point(alpha=0.7)
ggplot(ant.df, aes(south, SoilTempAnomaly)) + geom_point(alpha=0.7)
ggplot(ant.df, aes(topo, SoilTempAnomaly)) + geom_point(alpha=0.7)


ant.df %>% group_by(BDM) %>% 
  summarise(SoilTempSD=sd(SoilTempAnomaly, na.rm=T), 
            date=first(SampleDate), 
            mnEl=mean(mnt25, na.rm=T),
            sdEl=sd(mnt25, na.rm=T)) %>% 
  ggplot(aes(date, SoilTempSD, colour=mnEl)) + stat_smooth(se=F, span=2) + 
  geom_point(size=3) + scale_colour_viridis_c()
ant.df %>% group_by(BDM) %>% 
  summarise(SoilTempSD=sd(SoilTempAnomaly, na.rm=T), 
            date=first(SampleDate), 
            mnEl=mean(mnt25, na.rm=T),
            sdEl=sd(mnt25, na.rm=T)) %>% 
  ggplot(aes(date, SoilTempSD, colour=sdEl)) + stat_smooth(se=F, span=2) + 
  geom_point(size=3) + scale_colour_viridis_c()
ant.df %>% group_by(BDM) %>% 
  summarise(SoilTempSD=sd(SoilTempAnomaly, na.rm=T), 
            date=first(SampleDate), 
            mnEl=mean(mnt25, na.rm=T),
            sdEl=sd(mnt25, na.rm=T)) %>% 
  ggplot(aes(sdEl, SoilTempSD, colour=mnEl)) + stat_smooth(se=F, span=2) + 
  geom_point(size=3) + scale_colour_viridis_c()
ant.df %>% group_by(BDM) %>% 
  summarise(SoilTempSD=sd(SoilTempAnomaly, na.rm=T), 
            date=first(SampleDate), 
            mnEl=mean(mnt25, na.rm=T),
            sdEl=sd(mnt25, na.rm=T)) %>% 
  ggplot(aes(mnEl, SoilTempSD, colour=sdEl)) + stat_smooth(se=F, span=2) + 
  geom_point(size=3) + scale_colour_viridis_c()

ggplot(ant.df, aes(SampleDate, SoilTempAnomaly, colour=region)) + 
  geom_hline(yintercept=0, colour="gray") + geom_rug(alpha=0.25) +
  geom_point(size=3, shape=1) + #stat_smooth(se=F, span=2, size=0.5) +
  scale_colour_viridis_d() + facet_wrap(~Categorie) +
  theme(panel.grid.minor=element_blank())

ggplot(ant.df, aes(mnt25, SoilTempAnomaly, colour=region)) + 
  geom_hline(yintercept=0, colour="gray") + geom_rug(alpha=0.25) +
  geom_point(size=3, shape=1) + stat_smooth(se=F, method="lm", size=0.5) +
  scale_colour_viridis_d() + facet_wrap(~Categorie) +
  theme(panel.grid.minor=element_blank())

ggplot(ant.df, aes(mnt25, slope25, colour=region)) + geom_rug(alpha=0.25) +
  geom_point(size=3, shape=1) + #stat_smooth(se=F, span=2, size=0.5) +
  scale_colour_viridis_d() + facet_wrap(~Categorie) +
  theme(panel.grid.minor=element_blank())

ggplot(ant.df, aes(mnt25, topo, colour=region)) + 
  geom_hline(yintercept=0, colour="gray") + geom_rug(alpha=0.25) +
  geom_point(size=3, shape=1) +# stat_smooth(se=F, span=2, size=0.5) +
  scale_colour_viridis_d() + facet_wrap(~Categorie) +
  theme(panel.grid.minor=element_blank())

ggplot(ant.df, aes(mnt25, south, colour=region)) + geom_rug(alpha=0.25) +
  geom_point(size=3, shape=1) + 
  scale_colour_viridis_d() + facet_wrap(~Categorie) +
  theme(panel.grid.minor=element_blank())




ggplot(filter(ant.df, Categorie=="Autre"), aes(mnt25, TypeOfOpen)) + 
  geom_point(alpha=0.5)




# Pulling open types from samples
if(overwrite) {
  open_types_id <- ant.df %>% filter(Categorie=="Autre") %>% 
    select(Plot_id, TypeOfOpen) %>% filter(str_sub(Plot_id,5) != "00") %>% 
    arrange(Plot_id, TypeOfOpen) %>% group_by(Plot_id) %>% 
    summarise(Type=first(TypeOfOpen))
  plot_i$TypeOfOpen <- open_types_id$Type[match(plot_i$Plot_id, open_types_id$Plot_id)]
  write_csv(select(plot_i, Plot_id, TypeOfOpen), "~/Desktop/opfo_qc/plot_openTypes.csv")
  
  # need to re-pad 0's for BDM_id, Plot_id, Plot_id_orig because excel sucks
  # reload plot_i
  plot_i <- read_csv("data/opfo_envDataProcessed.csv") %>%
    mutate(Plot_id=str_pad(Plot_id, 6, "left", 0),
           Plot_id_orig=str_pad(Plot_id, 6, "left", 0),
           BDM_id=str_pad(BDM_id, 2, "left", 0))
  write_csv(plot_i, "data/opfo_envDataProcessed.csv") 
}




plot_i %>% filter(Categorie=="Autre") %>% group_by(BDM, TypeOfOpen, region) %>% 
  summarise(nOpenType=n(), mnEl=mean(el_mean)) %>%
  ungroup %>% group_by(BDM) %>% mutate(propOpenType=nOpenType/sum(nOpenType)) %>%
  ggplot(aes(mnEl, propOpenType, colour=TypeOfOpen)) + geom_point() + 
  stat_smooth(se=F) + facet_wrap(~region)


ant_plots.df <- ant.df %>% group_by(Plot_id) %>%
  select(c(3,21:56,58,88)) %>%
  summarise_all(first)
ant_plots.df <- ant.df %>% group_by(Plot_id) %>%
  select(c(3,21:56,58,88)) %>%
  summarise_all(first)

plot_ant_sum <- ant.df %>% group_by(Plot_id) %>%
  summarise(nTubes=n(), nGen=n_distinct(GENUSID, na.rm=T), 
            nSp=n_distinct(SPECIESID, na.rm=T)) %>%
  left_join(plot.sf, ., by="Plot_id") %>%
  mutate(nTubes=replace_na(nTubes, 0),
         nGen=replace_na(nGen, 0),
         nSp=replace_na(nSp, 0)) %>% 
  left_join(., ant_plots.df, by="Plot_id")
ggplot(plot_ant_sum, aes(el_plot, nTubes)) + 
  geom_point(alpha=0.5) + stat_smooth(se=F, method="lm") +
  facet_wrap(~Categorie)
ggplot(plot_ant_sum, aes(el_plot, nGen)) + 
  geom_point(alpha=0.5) + stat_smooth(se=F, method="lm") +
  facet_wrap(~Categorie)
ggplot(plot_ant_sum, aes(el_plot, nSp)) + 
  geom_point(alpha=0.5) + stat_smooth(se=F, method="lm") +
  facet_wrap(~Categorie)

site_ant_sum <- ant.df %>% filter(TypeOfSample=="soil") %>% group_by(BDM_id) %>%
  summarise(nTubes=n(), nGen=n_distinct(GENUSID, na.rm=T), 
            nSp=n_distinct(SPECIESID, na.rm=T)) %>%
  left_join(site_summary, ., by="BDM_id") %>%
  mutate(nTubes=replace_na(nTubes, 0),
         nGen=replace_na(nGen, 0),
         nSp=replace_na(nSp, 0)) %>%
  left_join(., site.box.sf %>% select(BDM), by="BDM") %>% 
  st_as_sf
  

ggplot(site_ant_sum, aes(el_mean, nTubes/n_plots, colour=region)) + 
  stat_smooth(method="lm") + geom_point() 
ggplot(site_ant_sum, aes(el_mean, nGen/n_plots, colour=region)) + 
  stat_smooth(method="lm") + geom_point() 
ggplot(site_ant_sum, aes(el_mean, nSp/n_plots, colour=region)) + 
  stat_smooth(method="lm") + geom_point() 


gen_ant_sum <- ant.df %>% group_by(GENUSID) %>% 
  summarise(nPlots=n_distinct(Plot_id), nBDM=n_distinct(BDM_id))
sp_ant_sum <- ant.df %>% group_by(SPECIESID) %>% 
  summarise(nPlots=n_distinct(Plot_id), nBDM=n_distinct(BDM_id))


ggplot(site_ant_sum, aes(MAT_mean, nSp)) + geom_point() + stat_smooth(method="lm")

ggplot(plot_ant_sum, aes(Grass, nGen, group=BDM)) + facet_wrap(~region) +
  geom_point(alpha=0.1) +
  geom_line(stat="smooth", method="lm", se=F, colour="black", alpha=0.5)
ggplot(plot_ant_sum, aes(Grass+Forb, nGen, group=BDM)) + facet_wrap(~region) +
  geom_point(alpha=0.1) +
  geom_line(stat="smooth", method="lm", se=F, colour="black", alpha=0.5)
ggplot(plot_ant_sum, aes(SoilTemp, nGen, group=BDM)) + facet_wrap(~region) +
  geom_point(alpha=0.1) +
  geom_line(stat="smooth", method="lm", se=F, colour="black", alpha=0.5)
ggplot(plot_ant_sum, aes(Bare, nGen, group=BDM)) + facet_wrap(~region) +
  geom_point(alpha=0.1) +
  geom_line(stat="smooth", method="lm", se=F, colour="black", alpha=0.5)
ggplot(plot_ant_sum, aes(Litter, nGen, group=BDM)) + facet_wrap(~region) +
  geom_point(alpha=0.1) +
  geom_line(stat="smooth", method="lm", se=F, colour="black", alpha=0.5)
ggplot(plot_ant_sum, aes(SoilTempAnomaly, nGen, group=BDM)) + facet_wrap(~region) +
  geom_point(alpha=0.1) +
  geom_line(stat="smooth", method="lm", se=F, colour="black", alpha=0.5)
ggplot(plot_ant_sum, aes(topo, nGen, group=BDM)) + facet_wrap(~region) +
  geom_point(alpha=0.1) +
  geom_line(stat="smooth", method="lm", se=F, colour="black", alpha=0.5)


# clim.f <- dir(paste0(gis.dir, "chelsa"))
# clim.key <- read_csv(paste0(gis.dir, "bioclim_variable_names.csv"))
# clim <- vector("list", length(clim.f)) %>%
#   setNames(str_sub(clim.f, start=8, end=15))
# for(i in 1:length(clim)) {
#   clim[[i]] <- raster::raster(dir(paste0(gis.dir, "chelsa"), full.names=T)[i]) %>%
#     raster::projectRaster(., dem_VD) %>%
#     raster::crop(., dem_VD)
# }




# for(i in 1:length(clim)) {
#   plot.sf <- plot.sf %>% mutate(NEW=raster::extract(clim[[i]], .))
#   names(plot.sf)[length(plot.sf)] <- paste0(names(clim)[i], "_plot")
# }

# Environmental data

site_habitats <- plot_i %>% 
  group_by(BDM_id, BDM, Categorie) %>%
  summarise(n_plots=n()) %>%
  spread(Categorie, n_plots)
site_habitats[is.na(site_habitats)] <- 0
site.box.sf$region <- plot_i$region[match(site.box.sf$BDM, plot_i$BDM)]

p <- ggplot() + layer_spatial(data=dem_VD) + 
  geom_sf(data=VD, fill=NA, colour="gray60", size=0.1) +
  geom_sf(data=site.box.sf, aes(colour=region), size=3) + 
  scale_fill_gradient("Elevation (m)", low="black", high="white") +
  scale_colour_brewer("Region", type="qual", palette=2)
ggsave("eda/regions.pdf", p, width=9, height=7)


# Ant data
ant.df <- gs_read(ss=gs_title("Scientific Ant Inventory VD"), ws="Transfer") %>%
  filter(!is.na(lon)) %>% rename(TubeNo=CATALOGUENUMBER, Plot_id=PLOT) %>%
  left_join(., select(plot_i, -c(1,19)), by="Plot_id")
ant.sf <- st_as_sf(ant.df, coords=c("lon", "lat"))


#-- spatial data
ants.sf <- st_read("data/gis/echantillons.shp") %>%  # issue with crs
  st_transform(st_crs(sites.sf)) %>%
  st_join(., site.box.sf, join=st_within) %>%
  filter(!is.na(BDM) & TubeNo >= 9990000) %>%
  mutate(BDM_id=site_summary$BDM_id[match(.$Id, site_summary$BDM)],
         el_tube=raster::extract(dem, .),
         slope_tube=raster::extract(slope, .)) %>%
  rename(TUBE=TubeNo) %>%
  select(TUBE, TypeEchant, MilieuOuve, el_tube, slope_tube)
for(i in 1:length(clim)) {
  ants.sf <- ants.sf %>% mutate(NEW=raster::extract(clim[[i]], .)) 
  names(ants.sf)[length(ants.sf)] <- paste0(names(clim)[i], "_tube")
}
ants.df <- read_csv("data/opfo_antSamplesProcessed.csv") %>%
  rename(Plot_id=PLOTCORRECTED) %>%
  select(-c(BDMCORRECTED, BDM)) %>%
  full_join(., st_set_geometry(ants.sf, NULL), by="TUBE") %>%
  full_join(., plot_i, by="Plot_id")


plot_summary <- ants.df %>% 
  mutate(Type=(as.numeric(Plot_id) %% 100)>0) %>%
  group_by(BDM, Plot_id) %>%
  summarise(nTubeSoil=sum(!is.na(TUBECORRECTED) & Type==1), 
            nTubeTransect=sum(!is.na(TUBECORRECTED) & Type==0)) %>%
  full_join(., plot_i, by=c("BDM", "Plot_id")) %>%
  filter(!is.na(nTubeSoil))
plot_bioclim <- plot_summary %>% ungroup %>% mutate_at(34:51, scale) %>% 
  gather("bioclimVar", "bioclimValue", 34:51) %>% 
  mutate(bioclimName=clim.key$name[match(.$bioclimVar, 
                                         paste0(clim.key$code, "_plot"))])

site_ant_summary <- plot_summary %>%
  ungroup %>% group_by(BDM) %>%
  summarise(propAnts=sum(nTubeSoil>0)/n(),
            meanAnts=mean(nTubeSoil)) %>%
  full_join(., site_summary, by="BDM")
  
ggplot(filter(plot_bioclim, nTubeSoil>0), 
       aes(x=el_plot, y=bioclimValue, colour=region)) + 
  geom_point(alpha=0.5) + 
  facet_wrap(~bioclimName,  
             labeller=labeller(bioclimName=label_wrap_gen(20))) + 
  scale_colour_brewer(type="qual", palette=2) + 
  labs(x="Elevation (m)", y="Normalized Value")
ggsave("eda/bioclim_region.pdf", width=9, height=7)

ggplot(filter(plot_bioclim, nTubeSoil>0 & 
                bioclimVar %in% c("bio10_02_plot", "bio10_03_plot", 
                                  "bio10_04_plot", "bio10_07_plot")), 
       aes(x=el_plot, y=bioclimValue, colour=region)) + 
  geom_point(alpha=0.8) + 
  facet_wrap(~bioclimName, scales="free_y",
             labeller=labeller(bioclimName=label_wrap_gen(20))) + 
  scale_colour_brewer("Region", type="qual", palette=2) + 
  labs(x="Elevation (m)", y="Normalized Value")
ggsave("eda/bioclim_variability_region_unscaled.pdf", width=5.5, height=4)



ggplot(filter(plot_bioclim, nTubeSoil>0), 
       aes(x=bioclimValue, y=nTubeSoil, colour=region)) + 
  geom_point(alpha=0.5) + stat_smooth(se=F, span=2) + facet_wrap(~bioclimName)


ggplot(plot_summary, aes(x=el_plot, colour=Categorie)) + 
  geom_density(size=1) + scale_colour_manual(values=lc_cols) + ylim(0, 0.005)


ggplot(plot_summary, aes(x=el_plot, y=log(nTubeSoil+1))) +
  geom_point(alpha=0.5) + stat_smooth(se=F, span=2, method="loess") + 
  facet_wrap(~Categorie, scales="free_y") + ylim(0, NA) + 
  labs(x="Plot elevation (m)", y="Number of colonies")
ggplot(plot_summary, aes(x=bio10_02_plot, y=log(nTubeSoil+1))) +
  geom_point(alpha=0.5) + stat_smooth(se=F, span=2, method="loess") + 
  facet_wrap(~Categorie, scales="free_y") + ylim(0, NA) + 
  labs(x="Mean diurnal temp range", y="Number of colonies")
ggplot(plot_summary, aes(x=bio04_plot, y=log(nTubeSoil+1))) +
  geom_point(alpha=0.5) + stat_smooth(se=F, span=2, method="loess") + 
  facet_wrap(~Categorie, scales="free_y") + ylim(0, NA) + 
  labs(x="Temp seasonality", y="Number of colonies")
ggplot(plot_summary, aes(x=bio07_plot, y=log(nTubeSoil+1))) +
  geom_point(alpha=0.5) + stat_smooth(se=F, span=2, method="loess") + 
  facet_wrap(~Categorie, scales="free_y") + ylim(0, NA) + 
  labs(x="Temp annual range", y="Number of colonies")
ggplot(plot_summary, aes(x=bio15_plot, y=log(nTubeSoil+1))) +
  geom_point(alpha=0.5) + stat_smooth(se=F, span=2, method="loess") + 
  facet_wrap(~Categorie, scales="free_y") + ylim(0, NA) + 
  labs(x="Precip seasonality", y="Number of colonies")
ggplot(plot_summary, aes(x=slope_plot, y=log(nTubeSoil+1))) +
  geom_point(alpha=0.5) + stat_smooth(se=F, span=2, method="loess") + 
  facet_wrap(~Categorie, scales="free_y") + ylim(0, NA) + 
  labs(x="Plot slope", y="Number of colonies")
ggplot(plot_summary, aes(x=SoilTemp, y=log(nTubeSoil+1))) +
  geom_point(alpha=0.5, shape=1) + stat_smooth(se=F, span=2, method="loess") + 
  geom_rug(alpha=0.1) + facet_wrap(~Categorie) + 
  labs(x="Soil temperature (ºC)", y="Number of colonies")
ggplot(plot_summary, aes(x=Grass, y=nTubeSoil, colour=region)) +
  geom_point(alpha=0.5) + stat_smooth(se=F, span=2, method="loess") + 
  facet_wrap(~Categorie)
ggplot(plot_summary, aes(x=Forb, y=nTubeSoil, colour=region)) +
  geom_point(alpha=0.5) + stat_smooth(se=F, span=2, method="loess") + 
  facet_wrap(~Categorie)
ggplot(plot_summary, aes(x=Shrub, y=nTubeSoil, colour=region)) +
  geom_point(alpha=0.5) + stat_smooth(se=F, span=2, method="loess") + 
  facet_wrap(~Categorie)
ggplot(plot_summary, aes(x=Litter, y=nTubeSoil, colour=region)) +
  geom_point(alpha=0.5) + stat_smooth(se=F, span=2, method="loess") + 
  facet_wrap(~Categorie)
ggplot(plot_summary, aes(x=Bare, y=nTubeSoil, colour=region)) +
  geom_point(alpha=0.5) + stat_smooth(se=F, span=2, method="loess") + 
  facet_wrap(~Categorie)
ggplot(plot_summary, aes(x=Moss, y=nTubeSoil, colour=region)) +
  geom_point(alpha=0.5) + stat_smooth(se=F, span=2, method="loess") + 
  facet_wrap(~Categorie)
ggplot(plot_summary, aes(x=Grass+Forb+Shrub+Moss, y=nTubeSoil, colour=region)) +
  geom_point(alpha=0.5) + stat_smooth(se=F, span=2, method="loess") + 
  facet_wrap(~Categorie)


ggplot() + 
  geom_point(data=filter(plot_summary, nTubeSoil==0), 
             aes(x=el_plot, y=SoilTemp), size=2, shape=1, colour="gray80") +
  geom_point(data=filter(plot_summary, nTubeSoil>0), 
             aes(x=el_plot, y=SoilTemp, colour=nTubeSoil), size=2, alpha=0.8) +
  scale_colour_viridis_c() + facet_wrap(~region)

ggplot() + 
  geom_point(data=filter(plot_summary, nTubeSoil==0), 
             aes(x=Sample_date, y=SoilTemp), size=2, shape=1, colour="gray80") +
  geom_point(data=filter(plot_summary, nTubeSoil>0), 
             aes(x=Sample_date, y=SoilTemp, colour=nTubeSoil), size=2, alpha=0.8) +
  scale_colour_viridis_c() + facet_wrap(~region)


plot_summary <- plot_summary %>% group_by(BDM) %>% 
  mutate(mnSoilTemp=mean(SoilTemp, na.rm=T),
         SoilTemp_dev=SoilTemp-mnSoilTemp) %>%
  ungroup

ggplot(plot_summary, aes(x=as.factor(nTubeSoil), y=SoilTemp_dev)) + 
  geom_hline(yintercept=0, colour="gray60", linetype=2) +
  geom_point(aes(fill=Categorie), shape=21,
             alpha=0.5, position=position_jitter(width=0.1, height=0)) +
  geom_violin(fill=NA) +
  scale_fill_manual(values=lc_cols) +
  facet_wrap(~Categorie) + 
  labs(x="Number of Colonies", y="Plot Soil Temp Anomaly within Site")

ggplot(plot_summary, aes(group=nTubeSoil, fill=nTubeSoil, x=SoilTemp_dev)) + 
  geom_vline(xintercept=0, colour="gray60", linetype=2) +
  geom_density(alpha=0.2) + geom_rug(alpha=0.4) +
  scale_fill_viridis() +
  labs(x="Plot Soil Temp Anomaly within Site") +
  facet_grid(nTubeSoil~.) + theme(legend.position="none")

ggplot(plot_summary, aes(group=Categorie, fill=Categorie, x=SoilTemp_dev)) + 
  geom_vline(xintercept=0, colour="gray60", linetype=2) +
  geom_density(alpha=0.2) + geom_rug(alpha=0.4) +
  scale_fill_manual(values=lc_cols) +
  labs(x="Plot Soil Temp Anomaly within Site") +
  facet_grid(nTubeSoil~., scales="free_y")



ants.soil.lc <- ants.sf %>% st_set_geometry(NULL) %>% 
  filter(TypeEchant=="soil") %>% group_by(Categorie) %>% 
  summarise(ct_soil=n(), prop=ct_soil/nrow(ants.sf)) %>%
  full_join(., lc_summary %>% group_by(Categorie) %>%
              summarise(n_plots_lc=sum(n_plots)), by="Categorie")
ants.soil.site.lc <- ants.sf %>% st_set_geometry(NULL) %>% 
  filter(TypeEchant=="soil") %>% group_by(BDM_id, Categorie) %>% 
  summarise(ct_soil=n()) %>%
  full_join(., lc_summary, by=c("BDM_id", "Categorie"))
ants.soil.site.lc$ct_soil[is.na(ants.soil.site.lc$ct_soil)] <- 0

ants.tran.site <- ants.sf %>% st_set_geometry(NULL) %>% 
  filter(TypeEchant=="transect") %>% group_by(BDM_id) %>% 
  summarise(ct_tran=n()) %>%
  full_join(., site_summary, by="BDM_id")

ants.type <- ants.sf %>% st_set_geometry(NULL) %>% 
  group_by(TypeEchant, Categorie) %>% 
  summarise(ct_ants=n()) %>%
  ungroup %>% group_by(TypeEchant) %>% 
  mutate(prop_ants=ct_ants/sum(ct_ants)) 
ants.type$lc_prop <- prop.totals[match(ants.type$Categorie, names(prop.totals))]








########
## Maps
########

# Regions
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=region), size=4, shape=21) +
  scale_fill_viridis_d("Region", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map__regions.jpg", width=7, height=6, dpi=400)

# Colony density (soil)
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=ct_soil/n_plots/12), 
          size=4, shape=21) +
  scale_fill_viridis_c("Colonies/L", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_colonies_per_L.jpg", width=7, height=6, dpi=400)

# Colony density (transect)
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=ct_tran/(2*n_plots/25)), 
          size=4, shape=21) +
  scale_fill_viridis_c("Colonies/km", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_colonies_per_km.jpg", width=7, height=6, dpi=400)

# Colony density (transect)
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=ct_tran/(2*n_plots/25)), 
          size=4, shape=21) +
  scale_fill_viridis_c("Colonies/km", option="D", trans="log1p") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_colonies_per_km_logCol.jpg", width=7, height=6, dpi=400)

# Mean Annual Temp
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=MAT_mean/10), 
          size=4, shape=21) +
  scale_fill_viridis_c("MAT (ºC)", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_MAT.jpg", width=7, height=6, dpi=400)

# Temp Annual Range
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=TAR_mean/10), 
          size=4, shape=21) +
  scale_fill_viridis_c("TAR (ºC)", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_TAR.jpg", width=7, height=6, dpi=400)

# Elevation (mean)
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=el_mean), 
          size=4, shape=21) +
  scale_fill_viridis_c("Elevation\nmean (m)", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_el_mean.jpg", width=7, height=6, dpi=400)

# Elevation (sd)
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=el_stdev), 
          size=4, shape=21) +
  scale_fill_viridis_c("Elevation\nSD (m)", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_el_sd.jpg", width=7, height=6, dpi=400)

# Soil Temp (mean)
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=SoilTemp_mean), 
          size=4, shape=21) +
  scale_fill_viridis_c("Soil Temp\nmean (ºC)", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_SoilTemp_mean.jpg", width=7, height=6, dpi=400)

# Soil Temp (sd)
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=SoilTemp_sd), 
          size=4, shape=21) +
  scale_fill_viridis_c("Soil Temp\nsd (ºC)", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_SoilTemp_sd.jpg", width=7, height=6, dpi=400)

# Grass (mean)
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=Grass_mean), 
          size=4, shape=21) +
  scale_fill_viridis_c("Grass\nmean", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_Grass_mean.jpg", width=7, height=6, dpi=400)

# Forb (mean)
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=Forb_mean), 
          size=4, shape=21) +
  scale_fill_viridis_c("Forb\nmean", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_Forb_mean.jpg", width=7, height=6, dpi=400)

# Shrub (mean)
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=Shrub_mean), 
          size=4, shape=21) +
  scale_fill_viridis_c("Shrub\nmean", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_Shrub_mean.jpg", width=7, height=6, dpi=400)

# Bare (mean)
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=Bare_mean), 
          size=4, shape=21) +
  scale_fill_viridis_c("Bare\nmean", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_Bare_mean.jpg", width=7, height=6, dpi=400)

# Litter (mean)
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=Litter_mean), 
          size=4, shape=21) +
  scale_fill_viridis_c("Litter\nmean", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_Litter_mean.jpg", width=7, height=6, dpi=400)

# Moss (mean)
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=Moss_mean), 
          size=4, shape=21) +
  scale_fill_viridis_c("Moss\nmean", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_Moss_mean.jpg", width=7, height=6, dpi=400)

# Autre 
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=Autre/n_plots), 
          size=4, shape=21) +
  scale_fill_viridis_c("Autre", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_Autre.jpg", width=7, height=6, dpi=400)

# Foret Conifere 
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=ForetConifere/n_plots), 
          size=4, shape=21) +
  scale_fill_viridis_c("Conifer\nForest", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_ForetConifere.jpg", width=7, height=6, dpi=400)

# Foret Feuillus 
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=ForetFeuillus/n_plots), 
          size=4, shape=21) +
  scale_fill_viridis_c("Broadleaf\nForest", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_ForetFeuillus.jpg", width=7, height=6, dpi=400)

# Foret Mixe 
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, aes(fill=ForetMixe/n_plots), 
          size=4, shape=21) +
  scale_fill_viridis_c("Mixed\nForest", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_ForetMixe.jpg", width=7, height=6, dpi=400)

# Foret (all) 
ggplot() + geom_sf(data=VD, fill="gray95", colour="gray80") +
  geom_sf(data=site.ants.centroids, size=4, shape=21, 
          aes(fill=(ForetConifere+ForetFeuillus+ForetMixe)/n_plots)) +
  scale_fill_viridis_c("All\nForest", option="D") + 
  theme(panel.grid=element_blank()) + fonts
ggsave("eda/map_ForetAll.jpg", width=7, height=6, dpi=400)






########
## General patterns
########

##---- 1. Site-level -- colony density (soil) vs. ...
# Elevation
ggplot(ants.soil.site, aes(x=el_mean, y=ct_soil/n_plots/12)) + 
  ylim(0,NA) + stat_smooth(method="loess", span=2) +
  geom_point() + labs(x="Elevation (m)", y="Colonies/L") + fonts
ggsave("eda/el_colonies_per_L.jpg", width=7, height=7, dpi=400)

# Mean Annual Temperature
ggplot(ants.soil.site, aes(x=MAT_mean/10, y=ct_soil/n_plots/12)) + 
  ylim(0,NA) + stat_smooth(method="loess", span=2) +
  geom_point() + labs(x="Mean Annual Temp (ºC)", y="Colonies/L") + fonts
ggsave("eda/MAT_colonies_per_L.jpg", width=7, height=7, dpi=400)

# Temperature Annual Range
ggplot(ants.soil.site, aes(x=TAR_mean/10, y=ct_soil/n_plots/12)) + 
  ylim(0,NA) + stat_smooth(method="loess", span=2) +
  geom_point() + labs(x="Temp Annual Range (ºC)", y="Colonies/L") + fonts
ggsave("eda/TAR_colonies_per_L.jpg", width=7, height=7, dpi=400)

# Soil Temperature (mean)
ggplot(ants.soil.site, aes(x=SoilTemp_mean, y=ct_soil/n_plots/12)) + 
  ylim(0,NA) + stat_smooth(method="loess", span=2) +
  geom_point() + labs(x="Mean Soil Temp (ºC)", y="Colonies/L") + fonts
ggsave("eda/SoilTempMean_colonies_per_L.jpg", width=7, height=7, dpi=400)

# Habitat Proportions
ggplot(site.ants.long, aes(x=n_plots_habitat/n_plots, colour=el_mean,
                           y=ct_soil/n_plots/12)) +
  geom_point(alpha=0.8, size=2) + ylim(0,.2) + xlim(0,NA) +
  scale_colour_viridis_c("Elevation (m)") +
  facet_wrap(~Habitat, scales="free_x") +
  fonts + theme(legend.position=c(0.9,0.1),
               axis.text=element_text(size=10)) +
  labs(x="Proportion habitat", y="Colonies/L")
ggsave("eda/habProp_colonies_per_L.jpg", width=10, height=9, dpi=400)








ggplot(site.ants.sf, aes(SoilTemp_sd, ct_soil/n_plots/12, fill=region)) + 
  geom_point(size=3, shape=21, alpha=0.8) + scale_fill_viridis_d() +
  facet_wrap(~region)
ggplot(site.ants.sf, aes(SoilTemp_mean, ct_soil/n_plots/12, fill=region)) + 
  geom_point(size=3, shape=21, alpha=0.8) + scale_fill_viridis_d() +
  facet_wrap(~region)
ggplot(site.ants.sf, aes(el_stdev, ct_soil/n_plots/12, fill=region)) + 
  geom_point(size=3, shape=21, alpha=0.8) + scale_fill_viridis_d() +
  facet_wrap(~region)
ggplot(site.ants.sf, aes(Autre/n_plots, ct_soil/n_plots/12, fill=region)) + 
  geom_point(size=3, shape=21, alpha=0.8) + scale_fill_viridis_d() +
  facet_wrap(~region)









ggplot(ants.soil.site.lc, aes(x=el_mean, 
                              y=ct_soil/n_plots/12, colour=Categorie)) + 
  geom_point() + facet_wrap(~Categorie) + ylim(0, NA) +
  scale_colour_manual(values=lc_cols) + stat_smooth(se=F, span=2) +
  labs(x="Elevation", y="Colonies per L")



site.ants.long <- site.ants.sf %>% st_set_geometry(NULL) %>%
  gather(Habitat, n_plots_habitat, 42:56)




ggplot(site.ants.long, aes(x=el_mean, y=n_plots_habitat/n_plots, 
                           colour=ct_soil/n_plots/12)) +
  geom_point(alpha=0.8, size=2) + 
  scale_colour_viridis_c("Colonies/L") +
  facet_wrap(~Habitat) +
  theme(legend.position=c(0.9,0.1)) +
  labs(x="Elevation (m)", y="Proportion habitat")
ggsave("eda/habProp_by_el.jpg", width=9, height=7, dpi=400)


ggplot(site.ants.sf, aes(x=ForetFeuillus/n_plots, colour=el_mean,
                         y=ct_soil/n_plots/12)) +
  geom_point(alpha=0.9, size=3) + xlim(0,1) + ylim(0,.2) +
  scale_colour_viridis_c("Elevation (m)") +
  labs(x="Proportion broadleaf forest", y="Colonies/L")
ggplot(site.ants.sf, aes(x=ForetConifere/n_plots, colour=el_mean,
                         y=ct_soil/n_plots/12)) +
  geom_point(alpha=0.9, size=3) + xlim(0,1) + ylim(0,.2) +
  scale_colour_viridis_c("Elevation (m)") +
  labs(x="Proportion conifer forest", y="Colonies/L")
ggplot(site.ants.sf, aes(x=ForetMixe/n_plots, colour=el_mean,
                         y=ct_soil/n_plots/12)) +
  geom_point(alpha=0.9, size=3) + xlim(0,1) + ylim(0,.2) +
  scale_colour_viridis_c("Elevation (m)") +
  labs(x="Proportion conifer forest", y="Colonies/L")
ggplot(site.ants.sf, aes(x=(ForetConifere+ForetFeuillus+ForetMixe)/n_plots, 
                         colour=el_mean, y=ct_soil/n_plots/12)) +
  geom_point(alpha=0.9, size=3) + xlim(0,1) + ylim(0,.2) +
  scale_colour_viridis_c("Elevation (m)") +
  labs(x="Proportion forest", y="Colonies/L")
ggplot(site.ants.sf, aes(x=PrairieSeche/n_plots, colour=el_mean,
                         y=ct_soil/n_plots/12)) +
  geom_point(alpha=0.9, size=3) + xlim(0,1) + ylim(0,.2) +
  scale_colour_viridis_c("Elevation (m)") +
  labs(x="Proportion dry meadow", y="Colonies/L")
ggplot(site.ants.sf, aes(x=(Autre+PrairieSeche)/n_plots, colour=el_mean,
                         y=ct_soil/n_plots/12)) +
  geom_point(alpha=0.9, size=3) + xlim(0,1) + ylim(0,.2) +
  scale_colour_viridis_c("Elevation (m)") +
  labs(x="Proportion dry meadow or open", y="Colonies/L")
ggplot(site.ants.sf, aes(x=ZoneConstruite/n_plots, colour=el_mean,
                         y=ct_soil/n_plots/12)) +
  geom_point(alpha=0.9, size=3) + xlim(0,1) + ylim(0,.2) +
  scale_colour_viridis_c("Elevation (m)") +
  labs(x="Proportion built", y="Colonies/L")




ggplot(ants.soil.site, aes(x=MAT_mean/10, y=ct_soil/n_plots/12)) + 
  ylim(0,NA) + stat_smooth(method="loess", span=2) +
  geom_point() + labs(x="Mean Annual Temperature (ºC)", y="Colonies per L")
ggplot(ants.soil.site, aes(x=TAR_mean/10, y=ct_soil/n_plots/12)) + 
  ylim(0,NA) + stat_smooth(method="loess", span=2) +
  geom_point() + labs(x="Temperature Annual Range (ºC)", y="Colonies per L")
ggplot(ants.soil.site, aes(x=Grass_mean, y=ct_soil/n_plots/12)) + 
  ylim(0,NA) + stat_smooth(method="loess", span=2) +
  geom_point() + labs(x="Grass", y="Colonies per L")
ggplot(ants.soil.site, aes(x=Forb_mean, y=ct_soil/n_plots/12)) + 
  ylim(0,NA) + stat_smooth(method="loess", span=2) +
  geom_point() + labs(x="Forb", y="Colonies per L")
ggplot(ants.soil.site, aes(x=Grass_mean + Forb_mean + Shrub_mean, 
                           y=ct_soil/n_plots/12)) + 
  ylim(0,NA) + stat_smooth(method="loess", span=2) +
  geom_point() + labs(x="Grass + Forb + Shrub", y="Colonies per L")

ggplot(ants.tran.site, aes(x=el_mean, y=ct_soil/(2000*n_plots/25))) + 
  ylim(0,NA) + stat_smooth(method="loess", span=2) +
  geom_point() + labs(x="Elevation (m)", y="Colonies per m")





ggplot(ants.soil.site.lc, aes(x=MAT_mean/10, y=ct_soil/n_plots, colour=Categorie)) + 
  geom_point() + facet_wrap(~Categorie) + ylim(0, NA) +
  scale_colour_manual(values=lc_cols) + stat_smooth(se=F, span=2) +
  labs(x="Mean annual temp (site mean)", y="Ant colonies per 12L soil")
ggplot(ants.soil.site.lc, aes(x=TAR_mean/10, y=ct_soil/n_plots, colour=Categorie)) + 
  geom_point() + facet_wrap(~Categorie) + ylim(0, NA) +
  scale_colour_manual(values=lc_cols) + stat_smooth(se=F, span=2) +
  labs(x="Annual temp range (site mean)", y="Ant colonies per 12L soil")
ggplot(ants.soil.site.lc, aes(x=Tmin_mean/10, y=ct_soil/n_plots, colour=Categorie)) + 
  geom_point() + facet_wrap(~Categorie) + ylim(0, NA) +
  scale_colour_manual(values=lc_cols) + stat_smooth(se=F, span=2) +
  labs(x="Annual minimum temp (site mean)", y="Ant colonies per 12L soil")
ggplot(ants.soil.site.lc, aes(x=Tmax_mean/10, y=ct_soil/n_plots, colour=Categorie)) + 
  geom_point() + facet_wrap(~Categorie) + ylim(0, NA) +
  scale_colour_manual(values=lc_cols) + stat_smooth(se=F, span=2) +
  labs(x="Annual maximum temp (site mean)", y="Ant colonies per 12L soil")
ggplot(ants.soil.site.lc, aes(x=Sample_date, y=ct_soil/n_plots, colour=Categorie)) + 
  geom_point() + facet_wrap(~Categorie) + ylim(0, NA) +
  scale_colour_manual(values=lc_cols) + stat_smooth(se=F, span=2) +
  labs(x="Sample date", y="Ant colonies per 12L soil")
ggplot(ants.soil.site.lc, aes(x=SoilTemp_mean, y=ct_soil/n_plots, colour=Categorie)) + 
  geom_point() + facet_wrap(~Categorie) + ylim(0, NA) +
  scale_colour_manual(values=lc_cols) + stat_smooth(se=F, span=2) +
  labs(x="Soil temp (site mean)", y="Ant colonies per 12L soil")
ggplot(ants.soil.site.lc, aes(x=Grass_mean, y=ct_soil/n_plots, colour=Categorie)) + 
  geom_point() + facet_wrap(~Categorie) + ylim(0, NA) +
  scale_colour_manual(values=lc_cols) + stat_smooth(se=F, span=2) +
  labs(x="Grass Coverage (site mean)", y="Ant colonies per 12L soil")
ggplot(ants.soil.site.lc, aes(x=Forb_mean, y=ct_soil/n_plots, colour=Categorie)) + 
  geom_point() + facet_wrap(~Categorie) + ylim(0, NA) +
  scale_colour_manual(values=lc_cols) + stat_smooth(se=F, span=2) +
  labs(x="Forb Coverage (site mean)", y="Ant colonies per 12L soil")
ggplot(ants.soil.site.lc, aes(x=Shrub_mean, y=ct_soil/n_plots, colour=Categorie)) + 
  geom_point() + facet_wrap(~Categorie) + ylim(0, NA) +
  scale_colour_manual(values=lc_cols) + stat_smooth(se=F, span=2) +
  labs(x="Shrub Coverage (site mean)", y="Ant colonies per 12L soil")
ggplot(ants.soil.site.lc, aes(x=Bare_mean, y=ct_soil/n_plots, colour=Categorie)) + 
  geom_point() + facet_wrap(~Categorie) + ylim(0, NA) +
  scale_colour_manual(values=lc_cols) + stat_smooth(se=F, span=2) +
  labs(x="Bare ground Coverage (site mean)", y="Ant colonies per 12L soil")
ggplot(ants.soil.site.lc, aes(x=Litter_mean, y=ct_soil/n_plots, colour=Categorie)) + 
  geom_point() + facet_wrap(~Categorie) + ylim(0, NA) +
  scale_colour_manual(values=lc_cols) + stat_smooth(se=F, span=2) +
  labs(x="Litter Coverage (site mean)", y="Ant colonies per 12L soil")
ggplot(ants.soil.site.lc, aes(x=Moss_mean, y=ct_soil/n_plots, colour=Categorie)) + 
  geom_point() + facet_wrap(~Categorie) + ylim(0, NA) +
  scale_colour_manual(values=lc_cols) + stat_smooth(se=F, span=2) +
  labs(x="Moss Coverage (site mean)", y="Ant colonies per 12L soil")


