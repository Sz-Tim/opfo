# functions for cleaning and munging data


#' Clean irregularities in species names
#' 
#' @param ant.df dataframe with column SPECIESID
#' 
#' @return dataframe with species names fully standardized

clean_species_names <- function(ant.df) {
  
  # remove rows with SPECIESID as: 
  #   NA
  #   espece exotique
  #   Camp_herc/lign (5 obs; 30/122 for each spp)
  #   Lasi_emar/brun (2 obs; 488/82 for each spp)
  #   Lasi_nige/plat (15 obs; 1567/203 for each spp)
  #   Lasi_sp
  #   Lept_gred/musc (3 obs; 3/1 for each spp)
  ant.df <- ant.df %>% 
    filter(!is.na(SPECIESID)) %>%
    filter(SPECIESID != "Espèce exotique") %>%
    filter(SPECIESID != "Camp_herc/lign") %>%
    filter(SPECIESID != "Lasi_emar/brun") %>%
    filter(SPECIESID != "Lasi_nige/plat") %>%
    filter(SPECIESID != "Lasi_sp") %>%
    filter(SPECIESID != "Lept_gred/musc")
  
  # fix typos
  ant.df$SPECIESID[ant.df$SPECIESID=="LasI_para"] <- "Lasi_para"
  ant.df$SPECIESID[ant.df$SPECIESID=="Form_lugubris/paralugubris"] <- "Form_lugu/para"
  
  # merge morphospecies
  #   Lasius alienus group
  ant.df$SPECIESID[ant.df$SPECIESID=="Lasi_alie gr"] <- "Lasi_alie-GR"
  ant.df$SPECIESID[ant.df$SPECIESID=="Lasi_alie"] <- "Lasi_alie-GR"
  #   Temnothorax nylanderi group
  ant.df$SPECIESID[ant.df$SPECIESID=="Temn_nyla gr"] <- "Temn_nyla-GR"
  ant.df$SPECIESID[ant.df$SPECIESID=="Temn_nyla"] <- "Temn_nyla-GR"
  #   Temnothorax tuberum group
  ant.df$SPECIESID[ant.df$SPECIESID=="Temn_tube gr"] <- "Temn_tube-GR"
  ant.df$SPECIESID[ant.df$SPECIESID=="Temn_tube"] <- "Temn_tube-GR"
  #   Temnothorax unifasciatus group
  ant.df$SPECIESID[ant.df$SPECIESID=="Temn_unif gr"] <- "Temn_unif-GR"
  ant.df$SPECIESID[ant.df$SPECIESID=="Temn_unif"] <- "Temn_unif-GR"
  return(ant.df)
}





#' load ant data
#' 
#' @param structured TRUE; load structured samples?
#' @param public TRUE; load public inventory?
#' @param str_type "all"; type of structured samples to include (all, soil,
#'   transect, tree)
#' @param clean_spp FALSE; standardize species names according to settings in
#'   in clean_spp_names()
#' @param full_pb FALSE; include all original columns in .$pub
#' @param DNA_ID TRUE; update species names with genetic IDs?
#' 
#' @return list with each dataframe and a (simplified) combined dataframe

load_ant_data <- function(structured=TRUE, public=TRUE, 
                          str_type="all", clean_spp=FALSE, full_pub=FALSE,
                          DNA_ID=TRUE) {
  library(tidyverse); library(sf); library(googlesheets)
  # load genetic IDs 
  if(DNA_ID) {
    dna_dir <- "~/Documents/unil/opfo_main/1_opfo/data/DNA_ID_clean/"
    dna_f <- setNames(dir(dna_dir, full.names=T), 
                      str_split_fixed(dir(dna_dir), "_", 2)[,1])
    dna_ids <- map(dna_f, read_csv)
  }
  
  # load structured samples
  if(structured) {
    df_s <- gs_read(ss=gs_title("Scientific Ant Inventory VD"), ws="Transfer", 
                    col_types=cols(PLOT=col_character())) %>%
      filter(!is.na(lon)) %>% 
      rename(TubeNo=CATALOGUENUMBER, Plot_id=PLOT) %>%
      mutate(Plot_id=as.character(Plot_id),
             SampleDate=lubridate::dmy(SampleDate),
             GEN_abbr=str_split_fixed(SPECIESID, "_", 2)[,1]) %>%
      st_as_sf(coords=c("lon", "lat")) %>%
      st_set_crs(4326) %>% st_transform(21781)
    # assign genetic IDs
    for(i in seq_along(dna_ids)) {
      dna_genus_unkSpp <- filter(df_s, GEN_abbr==names(dna_ids)[i] &
                                   !TubeNo %in% dna_ids[[i]]$TubeNo)$TubeNo
      dna_str_tubes <- filter(dna_ids[[i]], 
                              grepl("999[0-9][0-9][0-9][0-9]", TubeNo))
      dna_genus_index <- match(dna_str_tubes$TubeNo, df_s$TubeNo)
      df_s$SPECIESID[dna_genus_index] <- dna_str_tubes$ID
      if(names(dna_ids)[i] %in% c("Camp", "Temn")) {
        # Use morphological IDs for Camp,Temno when genetic ID is unavailable
        # Remove only NA species or Genus spp
        dna_genus_unkSpp <- filter(df_s, TubeNo %in% dna_genus_unkSpp & 
                                     (is.na(SPECIESID) | 
                                        SPECIESID==names(dna_ids)[i]))
      } 
      # Remove tubes without IDs
      df_s <- filter(df_s, !TubeNo %in% dna_genus_unkSpp)
    }
    if(str_type!="all") df_s <- filter(df_s, TypeOfSample==str_type)
    if(clean_spp) df_s <- clean_species_names(df_s)
  }
  # load public inventory samples
  if(public) {
    df_p <- gs_read(ss=gs_url("https://docs.google.com/spreadsheets/d/19mDCH-A7mmNwelXypK7ixpfcQGzS4x3o7w4MJU8-hJc/edit?usp=sharing"), verbose=F, ws="Samples",
                    locale=readr::locale(decimal_mark=",")) %>%
      rename(TubeNo=CATALOGUENUMBER, SPECIESID=SPECISID) %>%
      mutate(SampleDate=lubridate::ymd(DATECOLLECTION))
    # assign genetic IDs
    for(i in seq_along(dna_ids)) {
      dna_genus_unkSpp <- filter(df_p, GENUSID==names(dna_ids)[i] &
                                   !TubeNo %in% dna_ids[[i]]$TubeNo)$TubeNo
      dna_pub_tubes <- filter(dna_ids[[i]], 
                              !grepl("999[0-9][0-9][0-9][0-9]", TubeNo))
      dna_genus_index <- match(dna_pub_tubes$TubeNo, df_p$TubeNo)
      df_p$SPECIESID[dna_genus_index] <- dna_pub_tubes$ID
      if(names(dna_ids)[i] %in% c("Camp", "Temn")) {
        # Use morphological IDs for Camp,Temno when genetic ID is unavailable
        # Remove only NA species or Genus spp
        dna_genus_unkSpp <- filter(df_p, TubeNo %in% dna_genus_unkSpp & 
                                     (is.na(SPECIESID) | 
                                        SPECIESID==names(dna_ids)[i]))
      } 
      # Remove tubes without IDs
      df_p <- filter(df_p, !TubeNo %in% dna_genus_unkSpp)
    }
    if(full_pub) {
      df_p <- rbind(df_p %>% filter(is.na(LATITUDE)) %>%
                      filter(!is.na(SWISSCOORDINATE_X)) %>%
                      st_as_sf(coords=paste0("SWISSCOORDINATE_", c("X", "Y")), 
                               remove=F) %>%
                      st_set_crs(21781) %>%
                      mutate(Long=SWISSCOORDINATE_X, Lat=SWISSCOORDINATE_Y), 
                    df_p %>% filter(!is.na(LONGITUDE)) %>%
                      st_as_sf(coords=c("LONGITUDE", "LATITUDE"),
                               remove=F) %>% 
                      st_set_crs(4326) %>% st_transform(21781) %>%
                      mutate(Long=LONGITUDE, Lat=LATITUDE))
    } else {
      df_p <- rbind(df_p %>% filter(is.na(LATITUDE)) %>%
                    filter(!is.na(SWISSCOORDINATE_X)) %>%
                    st_as_sf(coords=paste0("SWISSCOORDINATE_", c("X", "Y"))) %>%
                    select(TubeNo, SPECIESID, SampleDate) %>%
                    st_set_crs(21781), 
                  df_p %>% filter(!is.na(LONGITUDE)) %>%
                    st_as_sf(coords=c("LONGITUDE", "LATITUDE")) %>% 
                    select(TubeNo, SPECIESID, SampleDate) %>% 
                    st_set_crs(4326) %>% st_transform(21781))
    }
    if(clean_spp) df_p <- clean_species_names(df_p)
  }
  # combine
  if(exists("df_p") && exists("df_s")) {
    df_all <- rbind(df_p %>% select(TubeNo, SPECIESID, SampleDate) %>% 
                      mutate(source="p"),
                    df_s %>% select(TubeNo, SPECIESID, SampleDate) %>% 
                      mutate(source="s"))
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
                              clim=FALSE, topo=FALSE, npp=FALSE,
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
    if(npp) {
      npp <- raster::raster(paste0(gis.dir, "MODIS_2010-2019_VD_21781.tif"))
      site.sf <- site.sf %>%
        mutate(npp=raster::extract(npp, ., fun=mean, na.rm=T))
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