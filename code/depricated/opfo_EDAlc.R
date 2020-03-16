

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
lc_i <- readxl::read_xlsx("data/landcover_id.xlsx", 1)



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
site.lc.prop <- readxl::read_xlsx("data/site_summaries_25.xlsx", 1) %>%
  mutate(BDM_id=str_pad(BDM_id, 2, side="left", pad="0"),
         Sample_Date=lubridate::ymd(Sample_Date)) %>%
  mutate(Canopy=lc_i$Canopy[match(Categorie, lc_i$LC)])
  
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
  mutate(Plot_id=as.character(Plot_id),
         SampleDate=lubridate::dmy(SampleDate)) %>%
  left_join(., select(plot.sf, -c(3,21:30)), by="Plot_id") %>%
  mutate(Canopy=lc_i$Canopy[match(Categorie, lc_i$LC)])

ant.lc.prop <- ant.df %>% 
  group_by(BDM, Categorie) %>% summarise(nTubes=n()) %>%
  group_by(BDM) %>% mutate(propTubes=nTubes/sum(nTubes)) %>%
  full_join(select(site.lc.prop, -BDM_id, -Sample_Date), ., 
            by=c("BDM", "Categorie")) %>%
  mutate(propArea=pctArea/100) %>%
  left_join(., site_summary, by="BDM") 

ant.canopy.prop <- ant.df %>% 
  group_by(BDM, Canopy) %>% summarise(nTubes=n()) %>%
  group_by(BDM) %>% mutate(propTubes=nTubes/sum(nTubes)) %>%
  full_join(select(site.lc.prop, -BDM_id, -Sample_Date) %>%
              group_by(BDM, Canopy) %>%
              summarise(pctArea=sum(pctArea), 
                        soil_plots=sum(soil_plots),
                        trans_len=sum(trans_len)), ., 
            by=c("BDM", "Canopy")) %>%
  mutate(propArea=pctArea/100) %>%
  left_join(., site_summary, by="BDM") 

ant.cc <- st_read("data/gis/AntData.shp") %>% st_transform(21781) %>%
  mutate(mnt25=raster::extract(dem, .)) %>%
  rename(TubeNo=CATALOGUEN, SPECIESID=SPECISID)

ant.all <- rbind(ant.df %>% select(TubeNo, mnt25, SPECIESID) %>%
                   mutate(source="structured"),
                 ant.cc %>% select(TubeNo, mnt25, SPECIESID) %>%
                   st_set_geometry(NULL) %>%
                   mutate(source="public"))
maxTubes <- 30
ant.all$SPECIESID[ant.all$SPECIESID=="erraticum"] <- "Tapi_erra"
ant.all$SPECIESID[ant.all$SPECIESID=="Form_lugubris/paralugubris"] <- "Form_lugu/para"
ant.all$SPECIESID[ant.all$SPECIESID=="Form_lugu/para/prat"] <- "Form_lugu/para"
ant.all$SPECIESID[ant.all$SPECIESID=="Lasi_alie gr"] <- "Lasi_alie"
ant.all.sum <- ant.all %>% group_by(SPECIESID) %>%
  summarise(nTubes=n()) %>%
  mutate(nToPin=nTubes)
ant.all.sum$nToPin[ant.all.sum$nToPin > maxTubes] <- maxTubes
sum(ant.all.sum$nToPin)



pointcount = function(r, pts){
  # make a raster of zeroes like the input
  r2 = r
  r2[] = 0
  # get the cell index for each point and make a table:
  counts = table(raster::cellFromXY(r,pts))
  # fill in the raster with the counts from the cell index:
  r2[as.numeric(names(counts))] = counts
  return(r2)
}

rast.cc <- pointcount(dem_agg, st_coordinates(ant.cc))





ggplot(ant.lc.prop, aes(x=propTubes-propArea, colour=Categorie)) + 
  geom_density(size=1) + scale_colour_manual(values=lc_cols)


ggplot(ant.lc.prop, aes(x=el_mean, y=nTubes/soil_plots)) + 
  geom_hline(yintercept=0) + geom_point(alpha=0.5, aes(colour=region)) + 
  stat_smooth(method="loess", se=F, span=2) +
  facet_wrap(~Categorie) + theme(panel.grid=element_blank())

ggplot(ant.lc.prop, aes(x=el_mean, y=propTubes-propArea)) + 
  geom_hline(yintercept=0) + geom_point(alpha=0.5, aes(colour=region)) + 
  stat_smooth(method="loess", se=F, span=2) +
  facet_wrap(~Categorie) + theme(panel.grid=element_blank())

ggplot(ant.canopy.prop, aes(x=el_mean, y=nTubes/soil_plots)) + 
  geom_hline(yintercept=0) + geom_point(alpha=0.5, aes(colour=region)) + 
  stat_smooth(method="loess", se=F, span=2) +
  facet_wrap(~Canopy) + theme(panel.grid=element_blank())

ggplot(ant.canopy.prop, aes(x=SoilTemp_mean, y=nTubes/soil_plots)) + 
  geom_hline(yintercept=0) + geom_point(alpha=0.5, aes(colour=region)) + 
  stat_smooth(method="loess", se=F, span=2) +
  facet_wrap(~Canopy) + theme(panel.grid=element_blank())

ggplot(ant.canopy.prop, aes(x=MAT_mean, y=nTubes/soil_plots)) + 
  geom_hline(yintercept=0) + geom_point(alpha=0.5, aes(colour=region)) + 
  stat_smooth(method="loess", se=F, span=2) +
  facet_wrap(~Canopy) + theme(panel.grid=element_blank())

ggplot(ant.canopy.prop, aes(x=TAR_mean, y=propTubes-propArea)) + 
  geom_hline(yintercept=0) + geom_point(alpha=0.5, aes(colour=region)) + 
  stat_smooth(method="loess", se=F, span=2) +
  facet_wrap(~Canopy) + theme(panel.grid=element_blank())


