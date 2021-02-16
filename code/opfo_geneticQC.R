library(tidyverse); library(sf); library(googlesheets)


################################################################################
####---- Original correction with data/Tetramorium_public_mtDNA_ID.txt
################################################################################

########---- STEP THROUGH LOADING FUNCTION

# copied from 00_fn.R::load_ant_data()
structured=TRUE
public=TRUE
str_type="all"
clean_spp <- TRUE

# load genetic IDs for Tetramorium
tetr.id <- read_tsv("data/Tetramorium_public_mtDNA_ID.txt",
                    col_names=c("TubeNo", "SPECIESID"))
# load structured samples
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

# load public inventory samples
df_p <- gs_read(ss=gs_url("https://docs.google.com/spreadsheets/d/19mDCH-A7mmNwelXypK7ixpfcQGzS4x3o7w4MJU8-hJc/edit?usp=sharing"), verbose=F, ws="Samples",
                locale=readr::locale(decimal_mark=",")) %>%
  rename(TubeNo=CATALOGUENUMBER, SPECIESID=SPECISID) %>%
  mutate(SampleDate=lubridate::ymd(DATECOLLECTION))


  
  
  
  
########---- IDENTIFY MISMATCHES
tetr.clean <- tetr.id
  

# Find TubeNo in tetr.id, but not in database
not_in_df <- filter(tetr.id, !TubeNo %in% df_p$TubeNo)
map_dfr(not_in_df$TubeNo, ~filter(df_p, grepl(.x, TubeNo)) %>% 
          select(RANG, TubeNo, SPECIESID, GENUSID))
# 1407 in tetr.id
filter(tetr.id, TubeNo=="14072")
filter(df_p, GENUSID=="Tetr") %>%
  filter(grepl("[1,4,0,7][1,4,0,7][1,4,0,7][1,4,0,7]", TubeNo)) %>%
  filter(!TubeNo %in% tetr.id$TubeNo) %>%
  select(RANG, TubeNo, SPECIESID, GENUSID)
filter(tetr.id, TubeNo=="4107")
filter(df_p, TubeNo=="4107") %>% select(RANG, TubeNo, SPECIESID, GENUSID)
# 4460 in tetr.id
filter(tetr.id, TubeNo=="14460")
cat("Mis-transcribed:", 
    " 1407 = 4107", 
    " 4460 = 14460", sep="\n")
tetr.clean$TubeNo[tetr.clean$TubeNo=="1407"] <- "4107"
tetr.clean$TubeNo[tetr.clean$TubeNo=="4460"] <- "14460"


# Find TubeNo in both, but that are not Tetr spp in database
df_match_but_not_TetrSPID <- filter(df_p, TubeNo %in% tetr.clean$TubeNo) %>% 
  filter(!grepl("Tetr", SPECIESID)) %>% 
  select(RANG, TubeNo, SPECIESID, GENUSID)
df_match_but_not_TetrSPID %>% print.AsIs
df_match_but_not_TetrSPID %>% filter(GENUSID != "Tetr")

non_tetr <- filter(df_match_but_not_TetrSPID, GENUSID != "Tetr")$TubeNo
map_dfr(non_tetr, ~filter(df_p, grepl(.x, TubeNo)) %>% 
          select(RANG, TubeNo, SPECIESID, GENUSID))
# 9620 in tetr.id
filter(tetr.clean, TubeNo=="9620")
likely_9620s <- filter(df_p, GENUSID=="Tetr") %>% 
  filter(grepl("[0-9][9,6,2,0][9,6,2,0][9,6,2,0]", TubeNo)) %>%
  filter(!TubeNo %in% tetr.clean$TubeNo) %>%
  select(RANG, TubeNo, SPECIESID, GENUSID)
filter(tetr.clean, TubeNo=="8620")
cat("Mis-transcribed:", 
    " 1950 = 11950", 
    " 9620 = 8620",
    "Incompletely transcribed:", 
    " 12937 = 12937B", 
    " 12563 = 12563C", 
    " 12328 = 12328B", sep="\n")
tetr.clean$TubeNo[tetr.clean$TubeNo=="1950"] <- "11950"
tetr.clean$TubeNo[tetr.clean$TubeNo=="9620"] <- "8620"
tetr.clean$TubeNo[tetr.clean$TubeNo=="12937"] <- "12937B"
tetr.clean$TubeNo[tetr.clean$TubeNo=="12563"] <- "12563C"
tetr.clean$TubeNo[tetr.clean$TubeNo=="12328"] <- "12328B"


# Find TubeNo in database, but not in tetr.id
not_in_id <- filter(df_p, GENUSID=="Tetr" & !TubeNo %in% tetr.clean$TubeNo)
table(not_in_id$SPECIESID)
cat(sum(df_p$GENUSID == "Tetr") - nrow(tetr.clean), 
    "more tubes in database than ID file -- just not performed?")
write_csv(not_in_id, "~/Desktop/Tetr_tubes_without_mtDNA.csv")




########---- SAVE CLEAN FILE
write_csv(tetr.clean, "data/Tetramorium_public_mtDNA_ID_CLEAN.csv")













################################################################################
####---- Second pass with data/Tetramorium_ID_GL.txt
################################################################################

########---- STEP THROUGH LOADING FUNCTION

# copied from 00_fn.R::load_ant_data()
structured=TRUE
public=TRUE
str_type="all"
clean_spp <- TRUE

# load genetic IDs for Tetramorium
tetr.id <- read_delim("data/Tetramorium_ID_GL.txt", delim=" ") %>%
  rename(TubeNo=CATALOGUENUMBER)
# load structured samples
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

# load public inventory samples
df_p <- gs_read(ss=gs_url("https://docs.google.com/spreadsheets/d/19mDCH-A7mmNwelXypK7ixpfcQGzS4x3o7w4MJU8-hJc/edit?usp=sharing"), verbose=F, ws="Samples",
                locale=readr::locale(decimal_mark=",")) %>%
  rename(TubeNo=CATALOGUENUMBER, SPECIESID=SPECISID) %>%
  mutate(SampleDate=lubridate::ymd(DATECOLLECTION))






########---- IDENTIFY MISMATCHES
tetr.clean <- tetr.id


# Find TubeNo in tetr.id, but not in database
not_in_df <- filter(tetr.id, !TubeNo %in% df_p$TubeNo)
nrow(not_in_df)
sum(tetr.id$TubeNo %in% df_s$TubeNo)
tetr.id %>% filter(!TubeNo %in% df_p$TubeNo) %>% filter(!TubeNo %in% df_s$TubeNo)
# Mismatches already identified in previous round
cat("Mis-transcribed:", 
    " 1407 = 4107", 
    " 4460 = 14460", sep="\n")
tetr.clean$TubeNo[tetr.clean$TubeNo=="1407"] <- "4107"
tetr.clean$TubeNo[tetr.clean$TubeNo=="4460"] <- "14460"


# Find TubeNo in both, but that are not Tetr spp in database
df_match_but_not_TetrSPID <- filter(df_p, TubeNo %in% tetr.clean$TubeNo) %>% 
  filter(!grepl("Tetr", SPECIESID)) %>% 
  select(RANG, TubeNo, SPECIESID, GENUSID)
df_match_but_not_TetrSPID %>% print.AsIs
df_match_but_not_TetrSPID %>% filter(GENUSID != "Tetr")

non_tetr <- filter(df_match_but_not_TetrSPID, GENUSID != "Tetr")$TubeNo
map_dfr(non_tetr, ~filter(df_p, grepl(.x, TubeNo)) %>% 
          select(RANG, TubeNo, SPECIESID, GENUSID))
# 9620 in tetr.id
filter(tetr.clean, TubeNo=="9620")
likely_9620s <- filter(df_p, GENUSID=="Tetr") %>% 
  filter(grepl("[0-9][9,6,2,0][9,6,2,0][9,6,2,0]", TubeNo)) %>%
  filter(!TubeNo %in% tetr.clean$TubeNo) %>%
  select(RANG, TubeNo, SPECIESID, GENUSID)
filter(tetr.clean, TubeNo=="8620")
cat("Mis-transcribed:", 
    " 1950 = 11950", 
    " 9620 = 8620",
    "Incompletely transcribed:", 
    " 12937 = 12937B", 
    " 12563 = 12563C", 
    " 12328 = 12328B", sep="\n")
tetr.clean$TubeNo[tetr.clean$TubeNo=="1950"] <- "11950"
tetr.clean$TubeNo[tetr.clean$TubeNo=="9620"] <- "8620"
tetr.clean$TubeNo[tetr.clean$TubeNo=="12937"] <- "12937B"
tetr.clean$TubeNo[tetr.clean$TubeNo=="12563"] <- "12563C"
tetr.clean$TubeNo[tetr.clean$TubeNo=="12328"] <- "12328B"


# Repeat for df_s: Find TubeNo in both, but that are not Tetr spp in database
filter(df_s, TubeNo %in% tetr.clean$TubeNo) %>% 
  filter(!grepl("Tetr", SPECIESID)) %>% 
  select(TubeNo, SPECIESID, GENUSID)


# Find TubeNo in database, but not in tetr.id
not_in_id_p <- filter(df_p, GENUSID=="Tetr" & !TubeNo %in% tetr.clean$TubeNo)
table(not_in_id_p$SPECIESID, useNA="ifany")
cat(sum(df_p$GENUSID == "Tetr") - nrow(filter(tetr.clean, nchar(TubeNo) < 7)), 
    "more tubes in database than ID file -- just not performed?")

# Repeat for df_s: Find TubeNo in database, but not in tetr.id
not_in_id_s <- filter(df_s, GENUSID=="Tetramorium" & !TubeNo %in% tetr.clean$TubeNo)
table(not_in_id_s$SPECIESID, useNA="ifany")
cat(sum(df_s$GENUSID == "Tetramorium") - nrow(filter(tetr.clean, nchar(TubeNo) == 7)), 
    "more tubes in database than ID file -- just not performed?")


# Confirm that all tubes in id match either df_p or df_s
sum(tetr.clean$TubeNo %in% df_p$TubeNo) + sum(tetr.clean$TubeNo %in% df_s$TubeNo)
nrow(tetr.clean)

write_csv(bind_rows(not_in_id_p %>% 
                      select(TubeNo, SPECIESID), 
                    st_set_geometry(not_in_id_s, NULL) %>% 
                      select(TubeNo, SPECIESID)), 
          "~/Desktop/Tetr_tubes_without_mtDNA_2.csv")




########---- SAVE CLEAN FILE
write_csv(tetr.clean, "data/Tetramorium_public_mtDNA_ID_CLEAN.csv")


Tetr <- c("tot_p"=sum(df_p$GENUSID=="Tetr"),
          "tot_s"=sum(df_s$GENUSID=="Tetramorium"),
          "tot_id"=nrow(tetr.id),
          "p_in_id"=sum(df_p$TubeNo %in% tetr.clean$TubeNo), 
          "s_in_id"=sum(df_s$TubeNo %in% tetr.clean$TubeNo))
cat("Total Tetramorium in databases: ", 
    Tetr["tot_p"] + Tetr["tot_s"], "\n",
    "Total Tetramorium in new file:  ", 
    Tetr["tot_id"], "\n", 
    "  Mismatch: ", Tetr["tot_p"] + Tetr["tot_s"] - Tetr["tot_id"], 
    "\n\n", 
    Tetr['tot_p']-Tetr['p_in_id'], " public Tetramorium without IDs", "\n",
    Tetr['tot_s']-Tetr['s_in_id'], " structured Tetramorium without IDs", sep="")







################################################################################
####---- Third pass with data/Tetramorium_ID_GL.txt
################################################################################

########---- STEP THROUGH LOADING FUNCTION

# copied from 00_fn.R::load_ant_data()
structured=TRUE
public=TRUE
str_type="all"
clean_spp <- TRUE

# load genetic IDs for Tetramorium
tetr.id <- read_delim("data/Tetramorium_ID_GL_2.txt", delim=" ") %>%
  rename(TubeNo=CATALOGUENUMBER)
# load structured samples
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

# load public inventory samples
df_p <- gs_read(ss=gs_url("https://docs.google.com/spreadsheets/d/19mDCH-A7mmNwelXypK7ixpfcQGzS4x3o7w4MJU8-hJc/edit?usp=sharing"), verbose=F, ws="Samples",
                locale=readr::locale(decimal_mark=",")) %>%
  rename(TubeNo=CATALOGUENUMBER, SPECIESID=SPECISID) %>%
  mutate(SampleDate=lubridate::ymd(DATECOLLECTION))






########---- IDENTIFY MISMATCHES
tetr.clean <- tetr.id


# Find TubeNo in tetr.id, but not in database
not_in_df <- filter(tetr.id, !TubeNo %in% df_p$TubeNo)
nrow(not_in_df)
sum(tetr.id$TubeNo %in% df_s$TubeNo)
tetr.id %>% filter(!TubeNo %in% df_p$TubeNo) %>% filter(!TubeNo %in% df_s$TubeNo)
# Mismatches already identified in previous round
cat("Mis-transcribed:", 
    " 1407 = 4107", 
    " 4460 = 14460", sep="\n")
tetr.clean$TubeNo[tetr.clean$TubeNo=="1407"] <- "4107"
tetr.clean$TubeNo[tetr.clean$TubeNo=="4460"] <- "14460"


# Find TubeNo in both, but that are not Tetr spp in database
df_match_but_not_TetrSPID <- filter(df_p, TubeNo %in% tetr.clean$TubeNo) %>% 
  filter(!grepl("Tetr", SPECIESID)) %>% 
  select(RANG, TubeNo, SPECIESID, GENUSID)
df_match_but_not_TetrSPID %>% print.AsIs
df_match_but_not_TetrSPID %>% filter(GENUSID != "Tetr")

non_tetr <- filter(df_match_but_not_TetrSPID, GENUSID != "Tetr")$TubeNo
map_dfr(non_tetr, ~filter(df_p, grepl(.x, TubeNo)) %>% 
          select(RANG, TubeNo, SPECIESID, GENUSID))
# 9620 in tetr.id
filter(tetr.clean, TubeNo=="9620")
likely_9620s <- filter(df_p, GENUSID=="Tetr") %>% 
  filter(grepl("[0-9][9,6,2,0][9,6,2,0][9,6,2,0]", TubeNo)) %>%
  filter(!TubeNo %in% tetr.clean$TubeNo) %>%
  select(RANG, TubeNo, SPECIESID, GENUSID)
filter(tetr.clean, TubeNo=="8620")
cat("Mis-transcribed:", 
    " 1950 = 11950", 
    " 9620 = 8620",
    "Incompletely transcribed:", 
    " 12937 = 12937B", 
    " 12563 = 12563C", 
    " 12328 = 12328B", sep="\n")
tetr.clean$TubeNo[tetr.clean$TubeNo=="1950"] <- "11950"
tetr.clean$TubeNo[tetr.clean$TubeNo=="9620"] <- "8620"
tetr.clean$TubeNo[tetr.clean$TubeNo=="12937"] <- "12937B"
tetr.clean$TubeNo[tetr.clean$TubeNo=="12563"] <- "12563C"
tetr.clean$TubeNo[tetr.clean$TubeNo=="12328"] <- "12328B"


# Repeat for df_s: Find TubeNo in both, but that are not Tetr spp in database
filter(df_s, TubeNo %in% tetr.clean$TubeNo) %>% 
  filter(!grepl("Tetr", SPECIESID)) %>% 
  select(TubeNo, SPECIESID, GENUSID)


# Find TubeNo in database, but not in tetr.id
not_in_id_p <- filter(df_p, GENUSID=="Tetr" & !TubeNo %in% tetr.clean$TubeNo)
table(not_in_id_p$SPECIESID, useNA="ifany")
cat(sum(df_p$GENUSID == "Tetr") - nrow(filter(tetr.clean, nchar(TubeNo) < 7)), 
    "more tubes in database than ID file -- just not performed?")

# Repeat for df_s: Find TubeNo in database, but not in tetr.id
not_in_id_s <- filter(df_s, GENUSID=="Tetramorium" & !TubeNo %in% tetr.clean$TubeNo)
table(not_in_id_s$SPECIESID, useNA="ifany")
cat(sum(df_s$GENUSID == "Tetramorium") - nrow(filter(tetr.clean, nchar(TubeNo) == 7)), 
    "more tubes in database than ID file -- just not performed?")


# Genus mismatches -- some IDs in tetr.id are not Tetramorium
id_not_tetr <- filter(tetr.clean, !grepl("Tetr", ID, ))


# Confirm that all tubes in id match either df_p or df_s
sum(tetr.clean$TubeNo %in% df_p$TubeNo) + sum(tetr.clean$TubeNo %in% df_s$TubeNo)
nrow(tetr.clean)

write_csv(bind_rows(not_in_id_p %>% 
                      select(TubeNo, SPECIESID), 
                    st_set_geometry(not_in_id_s, NULL) %>% 
                      select(TubeNo, SPECIESID)), 
          "~/Desktop/Tetr_tubes_without_mtDNA_3.csv")




########---- SAVE CLEAN FILE
write_csv(tetr.clean, "data/Tetramorium_public_mtDNA_ID_CLEAN.csv")


Tetr <- c("tot_p"=sum(df_p$GENUSID=="Tetr"),
          "tot_s"=sum(df_s$GENUSID=="Tetramorium"),
          "tot_id"=nrow(tetr.id),
          "p_in_id"=sum(df_p$TubeNo %in% tetr.clean$TubeNo), 
          "s_in_id"=sum(df_s$TubeNo %in% tetr.clean$TubeNo))
cat("Total Tetramorium in databases: ", 
    Tetr["tot_p"] + Tetr["tot_s"], "\n",
    "Total Tetramorium in new file:  ", 
    Tetr["tot_id"], "\n", 
    "  Mismatch: ", Tetr["tot_p"] + Tetr["tot_s"] - Tetr["tot_id"], 
    "\n\n", 
    Tetr['tot_p']-Tetr['p_in_id'], " public Tetramorium without IDs", "\n",
    Tetr['tot_s']-Tetr['s_in_id'], " structured Tetramorium without IDs", sep="")
table(bind_rows(not_in_id_p %>% 
                  select(TubeNo, SPECIESID), 
                st_set_geometry(not_in_id_s, NULL) %>% 
                  select(TubeNo, SPECIESID))$SPECIESID, useNA="ifany")





