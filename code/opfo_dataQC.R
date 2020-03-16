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
lc_i <- readxl::read_xlsx("data/landcover_id.xlsx", 1)

# Environmental data
site_summary <- read_csv("data/opfo_siteSummaryProcessed.csv")
lc_summary <- read_csv("data/opfo_lcSummaryProcessed.csv")

# Ant data
#-- spatial data
ants.sf <- st_read("data/gis/echantillons.shp") %>%  # issue with crs
  st_join(., site.box.sf, join=st_within) %>%
  filter(!is.na(BDM) & TubeNo >= 9990000)
#-- identification data
ants.df <- readxl::read_xlsx("data/Scientific Inventory IDs_TMS.xlsx", 1) %>%
  select(BDMCORRECTED, PLOTCORRECTED, TUBECORRECTED, GENUSID, SPECISID) %>%
  mutate(TUBE=str_sub(TUBECORRECTED, end=7) %>% as.integer) %>% 
  filter(!is.na(TUBE)) %>% filter(BDMCORRECTED != "??") %>% 
  mutate(BDM=as.numeric(BDMCORRECTED),
         Habitat=str_sub(PLOTCORRECTED, start=3L, end=4L), 
         HabitatName=lc_i$LC[match(Habitat, lc_i$LC_ID)])

st_write(ants.sf, "data/gis/opfo_sci.shp", layer_options="SHPT=POINT", 
         delete_layer=T)
write_csv(ants.df, "data/opfo_antSamplesProcessed.csv")




########
## identify mismatches in tube numbers between ants.sf and ants.df
########

siteAbund.df <- full_join(ants.df %>% group_by(BDM) %>% 
                            summarise(n_tubes_df=n_distinct(TUBE)),
                          ants.sf %>% st_set_geometry(NULL) %>% #rename(BDM=Id) %>%
                            group_by(BDM) %>%
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










# Joined output from echantillon.shp (Olivier added environmental variables)
f_gdr <- "Scientific Ant Inventory VD.xlsx"
f_ssEnv <- "AntSamples_scientificSurvey_withEnvData.xlsx"
ssEnv.df <- readxl::read_xlsx(paste0("~/Desktop/", f_ssEnv), 1, na="<Null>")
id_tim <- readxl::read_xlsx("data/Scientific Inventory IDs_TMS.xlsx", 1)
id_gdr <- readxl::read_xlsx(paste0("~/Desktop/", f_gdr), 1)

ssEnv.df <- ssEnv.df %>% rename(TubeNo=`Tube No`) %>% 
  filter(TubeNo >= 9990000)

ssEnv.df %>% filter(Categorie=="Autre") %>%
  filter(is.na(`Type of open field`)) %>% 
  select(3:4) %>% arrange(TubeNo) %>% 
  left_join(., id_gdr %>% 
              mutate(TubeNo=as.numeric(str_sub(CATALOGUENUMBER, 1, 7))) %>%
              select(BDM, PLOT, TubeNo), by="TubeNo") %>%
  mutate(plotType=str_sub(PLOT, 3, 4)) %>%
  write_csv("~/Desktop/open_mismatch_postCombine.csv")


ssEnv.df %>% select(c(3:6)) %>% arrange(TubeNo) %>% 
  left_join(., id_gdr %>% 
              mutate(TubeNo=as.numeric(str_sub(CATALOGUENUMBER, 1, 7))) %>%
              select(BDM, PLOT, TubeNo), by="TubeNo") %>%
  mutate(plotType=str_sub(PLOT, 3, 4),
         LC_plot=lc_i$LC[match(plotType, lc_i$LC_ID)]) %>%
  filter(Categorie != LC_plot) %>%
  write_csv("~/Desktop/hab_mismatch_postCombine.csv")


ssEnv.df %>% select(c(3:6)) %>% arrange(TubeNo) %>% 
  left_join(., id_gdr %>% 
              mutate(TubeNo=as.numeric(str_sub(CATALOGUENUMBER, 1, 7))) %>%
              select(BDM, PLOT, TubeNo), by="TubeNo") %>%
  mutate(plotType=str_sub(PLOT, 3, 4),
         LC_plot=lc_i$LC[match(plotType, lc_i$LC_ID)]) %>%
  filter(LC_plot=="Autre") %>%
  filter(is.na(`Type of open field`)) %>% 
  write_csv("~/Desktop/plotAutre_noType_postCombine.csv")





# Incorporating Tanja's updates
cat("Rows in id_tim:", nrow(id_tim), 
    "\nRows in id_gdr:", nrow(id_gdr),
    "\nRows in ssEnv:", nrow(ssEnv.df))

cat("Tubes in id_tim:", n_distinct(id_tim$TUBECORRECTED), 
    "\nTubes in id_gdr:", n_distinct(id_gdr$CATALOGUENUMBER),
    "\nTubes in ssEnv:", n_distinct(ssEnv.df$TubeNo))

cat("NA tubes in id_tim:", sum(is.na(id_tim$TUBECORRECTED)), 
    "\nNA tubes in id_gdr:", sum(is.na(id_gdr$CATALOGUENUMBER)),
    "\nNA tubes in ssEnv:", sum(is.na(ssEnv.df$TubeNo)))

cat("Number of tubes from id_tim found in id_gdr:", 
    sum(id_tim$TUBECORRECTED %in% id_gdr$CATALOGUENUMBER),
    "\nNumber of tubes from id_gdr found in id_tim:", 
    sum(id_gdr$CATALOGUENUMBER %in% id_tim$TUBECORRECTED),
    "\nNumber of tubes from ssEnv found in id_tim:", 
    sum(ssEnv.df$TubeNo %in% id_tim$TUBECORRECTED))

cat("Tubes from id_tim NOT found in id_gdr:\n  ", 
    id_tim$TUBECORRECTED[!id_tim$TUBECORRECTED %in% id_gdr$CATALOGUENUMBER],
    "\nTubes from id_gdr NOT found in id_tim:\n  ", 
    id_gdr$CATALOGUENUMBER[!id_gdr$CATALOGUENUMBER %in% id_tim$TUBECORRECTED],
    "Tubes from id_tim NOT found in ssEnv:\n  ", 
    id_tim$TUBECORRECTED[!id_tim$TUBECORRECTED %in% ssEnv.df$TubeNo],
    "\nTubes from ssEnv NOT found in id_tim:\n  ", 
    ssEnv.df$TubeNo[!ssEnv.df$TubeNo %in% id_tim$TUBECORRECTED])


sum(id_tim$TUBECORRECTED == id_gdr$CATALOGUENUMBER)
sum(id_tim$BDMCORRECTED == id_gdr$BDM)
sum(id_tim$PLOTCORRECTED == id_gdr$PLOT, na.rm=T)


cat("Number of plots from id_tim found in id_gdr:", 
    sum(id_tim$PLOTCORRECTED %in% id_gdr$PLOT),
    "\nNumber of plots from id_gdr found in id_tim:", 
    sum(id_gdr$PLOT %in% id_tim$PLOTCORRECTED))

cat("Plots from id_tim NOT found in id_gdr:\n  ", 
    id_tim$PLOTCORRECTED[!id_tim$PLOTCORRECTED %in% id_gdr$PLOT],
    "\nPlots from id_gdr NOT found in id_tim:\n  ", 
    id_gdr$PLOT[!id_gdr$PLOT %in% id_tim$PLOTCORRECTED])

cat("Tubes from id_tim NOT aligned with id_gdr:\n  ", 
    id_tim$TUBECORRECTED[which(id_tim$PLOTCORRECTED != id_gdr$PLOT)],
    "\nTubes from id_gdr NOT aligned with id_tim:\n  ", 
    id_gdr$CATALOGUENUMBER[which(id_gdr$PLOT != id_tim$PLOTCORRECTED)])




# joining plot habitat 
ssEnv.corrected <- ssEnv.df %>% arrange(TubeNo) %>% 
  left_join(., id_gdr %>% 
              mutate(TubeNo=as.numeric(str_sub(CATALOGUENUMBER, 1, 7))) %>%
              select(BDM, PLOT, TubeNo), by="TubeNo") %>%
  mutate(LC_plot=lc_i$LC[match(str_sub(PLOT, 3, 4), lc_i$LC_ID)])
write.csv(ssEnv.corrected, "~/Desktop/ssEnvCorr.csv")



library(googlesheets)
gdr.gsh <- gs_title("Scientific Ant Inventory VD")
gs_ws_ls(gdr.gsh)
id.df <- gs_read(ss=gdr.gsh, ws="Transfer")


hist(id.df$mnt25)



id.df %>% group_by(BDM, PLOT) %>%
  summarise(n_colonies=n(), nGen=n_distinct(GENUSID), el=mean(mnt25)) %>%
  ggplot(aes(el, nGen)) + geom_point(alpha=0.2)

ggplot(id.df, aes(x=GENUSID, y=mnt25)) + 
  geom_point(alpha=0.5, shape=1, size=2) +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

ggplot(id.df, aes(x=SPECIESID, y=mnt25)) + 
  geom_point(alpha=0.5, shape=1, size=2) + 
  facet_wrap(~GENUSID, scales="free_x") +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))
ggplot(id.df, aes(x=SPECIESID, y=mnt25)) + 
  geom_boxplot(outlier.shape=1, outlier.size=0.5) + 
  facet_wrap(~GENUSID, scales="free_x") +
  theme(axis.text.x=element_text(angle=270, hjust=0, vjust=0.5))

par(mfrow=c(3,3))
for(i in c(30,36:61)) {
  plot(id.df[[55]], id.df[[i]], xlab="El", ylab=names(id.df)[i])
}




ggplot(id.df, aes(x=Categorie)) + geom_bar(stat="count") + facet_wrap(~GENUSID)


