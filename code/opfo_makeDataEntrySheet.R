# This script creates a csv for entering data from the data sheets
# Tim Szewczyk
# 2019 Sep 19

library(tidyverse)

# load plot distribution summaries
site_i <- readxl::read_xlsx("data/site_summaries_25.xlsx", 1) %>%
  mutate(BDM_id=str_pad(BDM_id, 2, side="left", pad="0"))
lc_i <- readxl::read_xlsx("landcover_id.xlsx", 1) %>%
  mutate(lc_id=str_pad(LC_ID, 2, side="left", pad="0"))

# data entry sheet needs columns for:
#   BDM_id
#   BDM
#   Categorie
#   Plot_id

data_sheet <- site_i %>% 
  group_by(BDM_id, BDM, Categorie) %>%
  uncount(soil_plots) %>%
  select(BDM_id, BDM, Categorie) %>% 
  mutate(Plot_id=paste0(BDM_id, 
                        lc_i$lc_id[match(Categorie, lc_i$LC)],
                        str_pad(row_number(), 2, side="left", pad="0")))
write_csv(data_sheet, "plot_envDataTemplate.csv")
