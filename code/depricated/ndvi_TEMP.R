library(sf)


gis.dir <- "~/Documents/unil/GIS/data/"
VD <- st_read(paste0(gis.dir, "Vaud_boundaries.shp"))
dem <- raster::raster(paste0(gis.dir, "aster/CH_21781.tif")) %>% 
  raster::crop(., VD)



ndvi.f <- dir(paste0(gis.dir, "ndvi/"), "NDVI")
ndvi.date <- str_sub(ndvi.f, start=15, end=22)
ndvi <- vector("list", length(ndvi.f)) %>% setNames(paste0("NDVI_", ndvi.date))
for(i in 1:length(ndvi.f)) {
  ndvi[[i]] <- raster::raster(paste0(gis.dir, "ndvi/", ndvi.f[i])) %>%
    raster::projectRaster(., dem) %>%
    raster::crop(., dem)
  
}


plot.sf <- st_read("~/Documents/unil/opfo_str_sampling/data/gis/opfo_soil_25per.shp") %>%
  mutate(el_plot=raster::extract(dem, .)) %>% 
  mutate(Plot_id_orig=str_squish(str_replace_all(plot_id, "_", "")))
for(i in 1:length(ndvi)) {
  plot.sf <- plot.sf %>% mutate(NEW=raster::extract(ndvi[[i]], .)) 
  names(plot.sf)[length(plot.sf)] <- paste0(names(ndvi)[i], "_plot")
}

site_summary <- read_csv("~/Documents/unil/opfo_str_sampling/data/opfo_siteSummaryProcessed.csv")
plot_i <- read_csv("~/Documents/unil/opfo_str_sampling/data/opfo_envDataProcessed.csv") %>%
  full_join(., dplyr::select(plot.sf, Plot_id_orig, el_plot,
                      contains("NDVI")) %>% 
              st_set_geometry(NULL), 
            by="Plot_id_orig")



plot.time <- plot_i %>% gather("Date", "NDVI", 31:40)

ggplot(filter(plot.time, !is.na(region)), 
       aes(x=el_plot, y=NDVI, group=region, colour=region)) +
  geom_point(shape=1, alpha=0.3) + stat_smooth(se=F, span=2) + 
  scale_colour_brewer(type="qual") + facet_wrap(~Date)

ggplot(filter(plot.time, !is.na(region)), 
       aes(x=el_plot, y=NDVI, group=Date, colour=Date)) +
  geom_point(shape=1, alpha=0.3) + stat_smooth(se=F, span=2) + 
  scale_colour_viridis_d() + facet_wrap(~region)
