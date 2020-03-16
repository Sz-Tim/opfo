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









library(tidyverse); library(sf)
theme_set(theme_bw()); source("code/lc_cols.R")
fonts <- theme(axis.text=element_text(size=14),
               axis.title=element_text(size=16),
               legend.text=element_text(size=12),
               legend.title=element_text(size=14),
               strip.text=element_text(size=14))

# GIS data
gis.dir <- "~/Documents/unil/GIS/data/"
VD <- st_read(paste0(gis.dir, "Vaud_boundaries.shp"))
dem <- raster::raster(paste0(gis.dir, "aster/CH_21781.tif"))
slope <- raster::raster(paste0(gis.dir, "aster/CH_slope_21781.tif"))
sites.sf <- st_read(paste0(gis.dir, "Fourmis_data_shp/MilieuxBDMdiss.shp")) %>%
  st_transform(3857)
site.box.sf <- st_read(paste0(gis.dir, "Fourmis_data_shp/BDMplus.shp")) %>%
  st_transform(st_crs(sites.sf)) %>% mutate(BDM=Id)

# Environmental data
site_summary <- read_csv("data/opfo_site_summary_processed.csv")
lc_summary <- read_csv("data/opfo_lc_summary_processed.csv")
plot_i <- read_csv("plot_env_data_template.csv") %>%
  group_by(BDM_id, BDM, Categorie) %>%
  summarise(n_plots=n()) %>%
  spread(Categorie, n_plots)
plot_i[is.na(plot_i)] <- 0

# Ant data
#-- spatial data
ants.sf <- st_read("data/gis/opfo_sci.shp") %>%
  mutate(BDM_id=site_summary$BDM_id[match(.$Id, site_summary$BDM)])
site.ants.sf <- site.box.sf %>% 
  full_join(., ants.sf %>% st_set_geometry(NULL) %>% 
              filter(TypeEchant=="soil") %>% 
              group_by(BDM_id) %>% summarise(ct_soil=n()) %>%
              full_join(., site_summary, by="BDM_id") %>%
              select(-c(5:16)), by="BDM") %>%
  full_join(., ants.sf %>% st_set_geometry(NULL) %>% 
              filter(TypeEchant=="transect") %>% 
              group_by(BDM_id) %>% summarise(ct_tran=n()) %>%
              full_join(., site_summary, by="BDM_id") %>%
              select(c(ct_tran, BDM)), by="BDM") %>%
  full_join(., plot_i, by="BDM", suffix=c("", "_plots"))
site.ants.sf$ct_tran[is.na(site.ants.sf$ct_tran)] <- 0
site.ants.centroids <- site.ants.sf %>% st_centroid()
#-- identification data



scales::trans_new("log_p1", 
                  function(x) {log(x+1)},
                  function(x) {exp(x)-1})


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



ggplot(ants.soil.lc, aes(x=Categorie, y=ct_soil/n_plots_lc, fill=Categorie)) + 
  geom_bar(stat="identity") + 
  geom_hline(yintercept=sum(ants.soil.lc$ct_soil)/sum(ants.soil.lc$n_plots_lc), 
             colour="black", linetype=2) +
  geom_hline(yintercept=0, colour="gray30") +
  scale_fill_manual(values=lc_cols) +
  labs(x="", y="Ant colonies per 12L soil") +
  theme(axis.text.x=element_text(angle=310, hjust=0, vjust=0.5))
ggplot(ants.soil.lc, aes(x=Categorie, y=ct_soil, fill=Categorie)) + 
  geom_bar(stat="identity") + 
  geom_hline(yintercept=0, colour="gray30") +
  scale_fill_manual(values=lc_cols) +
  labs(x="", y="Ant colonies (total)") +
  theme(axis.text.x=element_text(angle=310, hjust=0, vjust=0.5))




site.props <- sites %>% group_by(Id, Categorie) %>%
  summarise(propArea=sum(Shape_Area)/1e6)
props.df <- site.props %>% st_set_geometry(NULL) %>%
  spread(Categorie, propArea)
props.df[is.na(props.df)] <- 0
prop.totals <- colSums(props.df[,-1])/sum(props.df[,-1])





ggplot() + geom_sf(data=sites.sf) + geom_sf(data=ants)

ggplot(ants, aes(Categorie)) + geom_bar() 

ggplot(filter(ants, is.na(Categorie))) + geom_sf()

ggplot(filter(ants.type, TypeEchant != "tree"), 
       aes(x=TypeEchant, y=(prop-lc_prop)/lc_prop*100)) + 
  geom_hline(yintercept=0) + geom_point() +
  facet_wrap(~Categorie) +
  labs(x="", y="% over-representation: (% ants - % LC) / % LC")
