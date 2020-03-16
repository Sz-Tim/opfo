# functions for cleaning and munging data


#' Clean irregularities in species names
#' 
#' @param ant.df dataframe with column SPECIESID
#' 
#' @return dataframe with species names fully standardized

clean_species_names <- function(ant.df) {
  ant.df <- ant.df %>% filter(!is.na(SPECIESID))
  ant.df$SPECIESID[ant.df$SPECIESID=="Form Copt"] <- "Form_copt"
  ant.df$SPECIESID[ant.df$SPECIESID=="Form_lugubris/paralugubris"] <- "Form_lugu/para"
  ant.df$SPECIESID[ant.df$SPECIESID=="Form_lugu/para/prat"] <- "Form_lugu/para"
  ant.df$SPECIESID[ant.df$SPECIESID=="Lasi_alie gr"] <- "Lasi_alie"
  ant.df$SPECIESID[ant.df$SPECIESID=="Temn_nyla gr"] <- "Temn_nyla"
  ant.df$SPECIESID[ant.df$SPECIESID=="Temn_tube gr"] <- "Temn_tube"
  ant.df$SPECIESID[ant.df$SPECIESID=="Temn_unif gr"] <- "Temn_unif"
  ant.df$SPECIESID[ant.df$SPECIESID=="Tetr_caes gr"] <- "Tetr_caes"
  return(ant.df)
}





#' load ant data
#' 
#' @param structured TRUE; load structured samples?
#' @param public TRUE; load public inventory?
#' @param str_type "all"; type of structured samples to include (all, soil,
#'   transect, tree)
#' 
#' @return list with each dataframe and a (simplified) combined dataframe

load_ant_data <- function(structured=TRUE, public=TRUE, 
                          str_type="all", clean_spp=FALSE) {
  library(tidyverse); library(sf); library(googlesheets)
  # load structured samples
  if(structured) {
    df_s <- gs_read(ss=gs_title("Scientific Ant Inventory VD"), ws="Transfer", 
                    col_types=cols(PLOT=col_character())) %>%
      filter(!is.na(lon)) %>% 
      rename(TubeNo=CATALOGUENUMBER, Plot_id=PLOT) %>%
      mutate(Plot_id=as.character(Plot_id),
             SampleDate=lubridate::dmy(SampleDate)) %>%
      st_as_sf(coords=c("lon", "lat")) %>%
      st_set_crs(4326) %>% st_transform(21781)
    if(str_type!="all") df_s <- filter(df_s, TypeOfSample==str_type)
    if(clean_spp) df_s <- clean_species_names(df_s)
  }
  # load public inventory samples
  if(public) {
    df_p <- gs_read(ss=gs_url("https://docs.google.com/spreadsheets/d/19mDCH-A7mmNwelXypK7ixpfcQGzS4x3o7w4MJU8-hJc/edit?usp=sharing"), verbose=F, ws="samples") %>%
      rename(TubeNo=CATALOGUENUMBER, SPECIESID=SPECISID)
    df_p <- rbind(df_p %>% filter(is.na(LATITUDE)) %>%
                    filter(!is.na(SWISSCOORDINATE_X)) %>%
                    st_as_sf(coords=paste0("SWISSCOORDINATE_", c("X", "Y"))) %>%
                    select(TubeNo, SPECIESID) %>%
                    st_set_crs(21781), 
                  df_p %>% filter(!is.na(LONGITUDE)) %>%
                    st_as_sf(coords=c("LONGITUDE", "LATITUDE")) %>% 
                    select(TubeNo, SPECIESID) %>% 
                    st_set_crs(4326) %>% st_transform(21781))
    if(clean_spp) df_p <- clean_species_names(df_p)
  }
  # combine
  if(exists("df_p") && exists("df_s")) {
    df_all <- rbind(df_p %>% mutate(source="p"),
                    df_s %>% select(TubeNo, SPECIESID) %>% mutate(source="s"))
  } else {
    df_all <- NULL
  }
  return(list(pub=df_p, str=df_s, all=df_all))
}




#' aggregate structured site data
#' 
#' @param gis.dir NULL; directory for storing public spreadsheet from dropbox
#' @param clim FALSE; load and summarise climate variables?
#' @param topo FALSE; load and summarise topographic variables?
#' @param opfo.dir "~/Documents/unil/opfo_str_sampling/"
#' @param goal "read"; is the purpose to 'write' a new site.sf or 'read' one?
#' 
#' @return sf object with site squares and all site-level environmental data

agg_str_site_data <- function(gis.dir="~/Documents/unil/opfo_main/2_gis/data/VD_21781/", 
                              clim=FALSE, topo=FALSE, 
                              opfo.dir="~/Documents/unil/opfo_main/1_opfo/data/",
                              goal="read") {
  library(tidyverse); library(sf)
  
  if(goal=="read") {
    site.sf <- st_read(paste0(gis.dir, "site_env_sf.shp"))
    names(site.sf) <- c(read_csv(paste0(gis.dir, "site_env_sf_names.csv"))$full,
                        "geometry")
  } else if(goal=="write") {
    lc_i <- readxl::read_xlsx(paste0(opfo.dir, "landcover_id.xlsx"), 1)
    site.lc.sum <- readxl::read_xlsx(paste0(opfo.dir, "site_summaries_25.xlsx"), 1) %>%
      mutate(BDM_id=str_pad(BDM_id, 2, side="left", pad="0")) %>%
      mutate(Canopy=lc_i$Canopy[match(Categorie, lc_i$LC)]) 
    site.props <- site.lc.sum %>% 
      group_by(BDM, Canopy) %>%
      summarise(area=sum(pctArea), nPlot=sum(soil_plots)) %>%
      pivot_wider(names_from=Canopy, values_from=c(area, nPlot), 
                  values_fill=list(area=0, nPlot=0)) %>%
      ungroup %>%
      mutate(area_Total=rowSums(select(., starts_with("area"))),
             nPlot_Total=rowSums(select(., starts_with("nPlot")))) %>%
      left_join(., site.lc.sum %>% group_by(BDM) %>%
                  summarise(Shann_LC=vegan::diversity(pctArea, index="shannon"),
                            Simps_LC=vegan::diversity(pctArea, index="simpson"),
                            N_LC=n()), 
                by="BDM")
    df_sum <- read_csv(paste0(opfo.dir, "opfo_siteSummaryProcessed.csv")) %>% 
      left_join(., site.props, by="BDM")
    site.sf <- st_read(paste0(gis.dir, "BDMplus_VD_21781.shp")) %>%
      left_join(., df_sum)
    if(clim) {
      clim.f <- dir(gis.dir, "chelsa") %>%
        setNames(., str_remove(., "_chelsa_VD_21781.tif"))
      clim.r <- map(clim.f, ~raster::raster(paste0(gis.dir, .)))
      site.sf <- site.sf %>%
        bind_cols(., map_dfc(clim.r, ~raster::extract(., site.sf, fun=mean)))
    } 
    if(topo) {
      dem <- raster::raster(paste0(gis.dir, "dem_VD_21781.tif"))
      slope <- raster::raster(paste0(gis.dir, "slope_VD_21781.tif"))
      site.sf <- site.sf %>%
        mutate(el=raster::extract(dem, ., fun=mean),
               slope=raster::extract(slope, ., fun=mean))
    }
    st_write(site.sf, paste0(gis.dir, "site_env_sf.shp"), delete_layer=TRUE)
    write_csv(data.frame(full=names(st_set_geometry(site.sf, NULL))), 
              paste0(gis.dir, "site_env_sf_names.csv"))
  }
  return(site.sf)
}
















pointcount <- function(r, pts){
  # make a raster of zeroes like the input
  r2 = r
  r2[] = 0
  # get the cell index for each point and make a table:
  counts = table(raster::cellFromXY(r,pts))
  # fill in the raster with the counts from the cell index:
  r2[as.numeric(names(counts))] = counts
  return(r2)
}