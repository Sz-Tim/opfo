

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
  mutate(elAnomaly=mnt25-mean(mnt25, na.rm=T)) %>%
  ungroup %>% group_by(GENUSID, SPECIESID) %>%
  mutate(minSpEl=min(mnt25, na.rm=T),
         maxSpEl=max(mnt25, na.rm=T)) 

interp.rng <- ant.df %>% group_by(SPECIESID) %>%
  summarise(minEl=min(mnt25, na.rm=T),
            maxEl=max(mnt25, na.rm=T),
            medEl=median(mnt25, na.rm=T)) %>%
  arrange(minEl, maxEl) %>%
  mutate(SpElOrder=row_number(),
         SPECIESID_el=factor(SPECIESID, levels=SPECIESID[SpElOrder]))
  
ggplot(interp.rng, aes(x=SPECIESID_el, ymin=minEl, ymax=maxEl)) +
  geom_linerange() + coord_flip() + 
  geom_point(aes(y=minEl), size=0.5) + geom_point(aes(y=maxEl), size=0.5)
  
el_seq <- seq(350, 2150, 100)
S <- rep(NA, length(el_seq))
for(i in seq_along(el_seq)) {
  S[i] <- sum(interp.rng$minEl <= el_seq[i] & interp.rng$maxEl >= el_seq[i])
}

plot(el_seq, S, type="b")



ant.df %>% 
  ggplot(aes(x=mnt25, y=SPECIESID)) + 
  facet_grid(GENUSID~., scales="free_y", space="free_y") +
  geom_point(alpha=0.4, size=3, shape=1) + geom_line()


ant.df %>% filter(GENUSID=="Myrmica") %>%
  ggplot(aes(x=mnt25, y=SPECIESID)) + 
  geom_point(alpha=0.4, size=3, shape=1) + geom_line()

ggplot() + geom_sf(data=VD, fill="white") + 
  geom_sf(data=site.box.sf, fill=NA, colour="gray30") +
  geom_point(data=ant.df, aes(X, Y), alpha=0.5, colour="red", size=1) +
  facet_wrap(~GENUSID)
  