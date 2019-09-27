# Opération Fourmis
# Data Quality Control
# Tim Szewczyk
# Initialized 2019 Sep 27

# This script performs quality control checks on the ant samples and environmental data collected for the Opération Fourmis scientifically structure sampling from 1 July - 5 Sep 2019. 

# The ant samples were (mostly) entered in the field directly into the GIS layer echantillons.shp. This layer should include spatial data for each point, the sample type, and the tube number. The physical tubes were (mostly) bundled into one bag per plot or transect at each site, with that information entered into a spreadsheet for identification.

# First, this script identifies discrepancies between the two datasets, including incorrectly recorded tube numbers, incorrectly recorded BDM square numbers, and incorrectly located points based on habitat type. The recorded locations had, at best, 4m accuracy, and there was variation across squares and collectors in whether the point location was manually edited to reflect the actual sampling location, or if the raw location was used. When there is a disagreement between the habitat type on GIS and the habitat of the plot recorded with the physical tube, I defer to the physical record as more accurate. Assuming the correct plot ID is also the only way to connect the ant samples to the local environmental data.


########
## set up
########

library(tidyverse); library(sf)

# Site & plot data
site.box.sf <- st_read("data/gis/BDMplus.shp") %>%
  st_transform(3857) %>% mutate(BDM=Id)
plot_i <- read_csv("bdm_datasheets/plot_envDataTemplate.csv") %>%
  group_by(BDM_id, BDM, Categorie) %>%
  summarise(n_plots=n()) %>%
  spread(Categorie, n_plots)
plot_i[is.na(plot_i)] <- 0

# Environmental data
site_summary <- read_csv("data/opfo_siteSummaryProcessed.csv")
lc_summary <- read_csv("data/opfo_lcSummaryProcessed.csv")

# Ant data
#-- spatial data
ants.sf <- st_read("data/gis/echantillons.shp") %>% #st_read("data/gis/opfo_sci.shp") %>%
  mutate(BDM_id=site_summary$BDM_id[match(.$Id, site_summary$BDM)]) %>%
  filter(TubeNo >= 9990000)
#-- identification data
ants.df <- readxl::read_xlsx("data/Scientific Inventory IDs.xlsx", 1) %>%
  select(BDMCORRECTED, PLOTCORRECTED, CATALOGUENUMBER, GENUSID, SPECISID) %>%
  mutate(TUBE=str_sub(CATALOGUENUMBER, end=7) %>% as.integer) %>% 
  filter(!is.na(TUBE)) %>% filter(BDMCORRECTED != "??") %>% 
  mutate(BDM=as.numeric(BDMCORRECTED))




########
## identify mismatches in tube numbers between ants.sf and ants.df
########

siteAbund.df <- full_join(ants.df %>% group_by(BDM) %>% 
                            summarise(n_tubes_df=n_distinct(TUBE)),
                          ants.sf %>% st_set_geometry(NULL) %>% 
                            rename(BDM=Id) %>% group_by(BDM) %>%
                            summarise(n_tubes_sf=n_distinct(TubeNo)), by="BDM")
dim(siteAbund.df) # 44 rows = 44 squares
siteAbund.df %>% filter(n_tubes_df != n_tubes_sf) # 13 squares different

site_mismatch <- tibble(BDM=integer(), Tube=integer(), 
                        FoundIn=character(), MissingFrom=character())
for(i in 1:n_distinct(siteAbund.df$BDM)) {
  BDM.i <- siteAbund.df$BDM[i]
  i.df <- filter(ants.df, BDM==BDM.i)
  i.sf <- filter(ants.sf, Id==BDM.i)
  i.df_not_sf <- i.df$TUBE[! i.df$TUBE %in% i.sf$TubeNo]
  i.sf_not_df <- i.sf$TubeNo[! i.sf$TubeNo %in% i.df$TUBE]
  if(length(i.df_not_sf) > 0) {
    site_mismatch <- site_mismatch %>% 
      add_row(BDM=BDM.i, Tube=i.df_not_sf, FoundIn="df", MissingFrom="sf")
  }
  if(length(i.sf_not_df) > 0) {
    site_mismatch <- site_mismatch %>% 
      add_row(BDM=BDM.i, Tube=i.sf_not_df, FoundIn="sf", MissingFrom="df")
  }
}

##--- Identify misaligned tube numbers
nrow(site_mismatch)  # 41
n_distinct(site_mismatch$Tube) # 41
site_mismatch %>% filter(! Tube %in% Tube[duplicated(Tube)]) %>% 
  mutate(FoundIn=as.factor(FoundIn)) %>% summary()
# 57-49=8; 16 (2*8) are assigned the wrong BDM (corrected)
# 20 in df, but not sf
# 21 in sf, but not df

n_distinct(ants.df$TUBE)    # 1431
n_distinct(ants.sf$TubeNo)  # 1440

# tubes with the wrong BDM (corrected -- empty df now)
site_mismatch %>% filter(Tube %in% Tube[duplicated(Tube)]) %>% arrange(Tube)

# tube mismatches by BDM (likely typo check)
site_mismatch %>% arrange(BDM, Tube) %>% print.AsIs()
# output copied to admin/opfo_pointLocationChecks.csv


df_not_sf <- ants.df %>% filter(!TUBE %in% ants.sf$TubeNo)
sf_not_df <- ants.sf %>% filter(!TubeNo %in% ants.df$TUBE) %>%
  select(TypeEchant, TubeNo, Categorie, Id, BDM_id) 

nrow(df_not_sf)  # 12
nrow(sf_not_df)  # 21








