# Figures for talks
# Tim Szewczyk

#--- libraries and options
pkgs <- c("tidyverse", "ggspatial", "sf", "viridis", "vegan", "iNEXT", "betapart")
suppressMessages(invisible(lapply(pkgs, library, character.only=T)))
theme_set(theme_bw() + theme(panel.grid.minor=element_blank())); 
walk(paste0("code/", c("lc_cols", "00_opfo_fn"), ".R"), source)
fonts <- theme(axis.text=element_text(size=14),
               axis.title=element_text(size=16),
               legend.text=element_text(size=12),
               legend.title=element_text(size=14),
               strip.text=element_text(size=14),
               panel.grid=element_blank())


#--- load data
# Vaud baundaries
VD_raw <- st_read("../2_gis/data/VD_21781/Vaud_boundaries.shp") %>%
  filter(!grepl("Lac ", NAME)) 
VD <- st_union(VD_raw)
# Topography
dem <- raster::raster("../2_gis/data/VD_21781/dem_VD_21781.tif") %>%
  raster::mask(., st_zm(VD_raw))
dem_bin <- table(dem@data@values %/% 100 * 100)
slope <- raster::raster("../2_gis/data/VD_21781/slope_VD_21781.tif") %>%
  raster::mask(., st_zm(VD_raw))
# Ants
ant <- load_ant_data(clean_spp=T)
# Plots and sites
plot.i <- read_csv("data/opfo_envDataProcessed.csv") %>%
  select(BDM, BDM_id, Categorie, Plot_id, TypeOfOpen, SoilTemp, 
         Grass, Forb, Shrub, Bare, Litter, Moss)
plot.sf <- st_read("../2_gis/data/opfo/opfo_soil_25per.shp") %>%
  mutate(Plot_id=str_replace(str_replace(plot_id, "_", ""), "_", ""),
         mnt25=raster::extract(dem, .)) %>%
  select(Plot_id, mnt25) %>% 
  left_join(., select(plot.i, Plot_id, BDM), by="Plot_id") %>%
  left_join(., ant$str %>% st_set_geometry(NULL) %>%
              select(c(3, 23:25, 29:47, 49:63)) %>% 
              group_by(Plot_id) %>% filter(row_number()==1), by="Plot_id") %>%
  filter(!is.na(mnt25))
site.sf <- full_join(agg_str_site_data(), 
                     plot.i %>% filter(Categorie=="Autre") %>% 
                       group_by(BDM, TypeOfOpen) %>% summarise(n_=n()) %>%
                       group_by(BDM) %>% mutate(pr_=n_/sum(n_)) %>%
                       pivot_wider(names_from="TypeOfOpen", 
                                   values_from=c("n_", "pr_"), 
                                   values_fill=list(n_=0, pr_=0),
                                   names_sep=""),
                     by="BDM") %>%
  mutate(Sample_date=lubridate::yday(lubridate::dmy(Sample_date)))
plot.sf <- plot.sf %>% mutate(region=site.sf$region[match(plot.sf$BDM, site.sf$BDM)])


#--- colors
blue <- "#0571b0"
red <- "#ca0020"
col_region <- c(Alps="#56B4E9", Jura="#E69F00", Low="#666666")
col_mtn <- c(Low="#666666", Mtn="#009E73")







########------------------------------------------------------------------------
## Sampling point maps
##

# Map of structured samples
pdf("eda/map_ants_BDM.pdf", height=7, width=7)
par(mar=c(0.5, 0.5, 0, 0))
raster::plot(dem, legend=F, axes=F, box=F,
             col=colorRampPalette(c("gray60", "white"))(255))
raster::scalebar(d=10000, xy=c(570000, 162000), below="km", 
                 label=c(0, 5, 10), type="bar")
plot(VD, add=TRUE, lwd=0.5)
plot(select(site.sf, region), add=TRUE, col=col_region[site.sf$region])
legend(569000, 180000, legend=c("Alps", "Jura", "Low"), fill=col_region, bty="n")
dev.off()

# Map of public samples
pdf("eda/map_ants_public.pdf", height=7, width=7)
par(mar=c(0.5, 0.5, 0, 0))
raster::plot(dem, legend=F, axes=F, box=F,
             col=colorRampPalette(c("gray60", "white"))(255))
raster::scalebar(d=10000, xy=c(570000, 162000), below="km", 
                 label=c(0, 5, 10), type="bar")
plot(VD, add=TRUE, lwd=0.5)
plot(select(ant$pub, TubeNo), add=TRUE, col=rgb(1,0,0,0.5), cex=0.75)
dev.off()

# Map of all samples
pdf("eda/map_ants_all.pdf", height=7, width=7)
par(mar=c(0.5, 0.5, 0, 0))
raster::plot(dem, legend=F, axes=F, box=F,
             col=colorRampPalette(c("gray60", "white"))(255))
raster::scalebar(d=10000, xy=c(570000, 162000), below="km", 
                 label=c(0, 5, 10), type="bar")
plot(VD, add=TRUE, lwd=0.5)
plot(select(ant$pub, TubeNo), add=TRUE, col=rgb(1,0,0,0.5), cex=0.75)
plot(select(site.sf, region), add=TRUE, col=NA, fill=NA)
dev.off()

# Map of binned elevation
pdf("eda/map_elBins.pdf", height=7, width=7)
par(mar=c(0.5, 0.5, 0, 0))
raster::plot(dem %/% 100 * 100, legend=F, axes=F, box=F, 
             col=terrain.colors(length(dem_bin)))
plot(VD, add=TRUE, lwd=0.5)
dev.off()








########------------------------------------------------------------------------
## Elevational distribution of sites
##

site.sf %>% mutate(region=forcats::lvls_reorder(region, 3:1)) %>%
  ggplot(aes(region, el, colour=region)) + fonts +
  geom_point(size=2, alpha=0.8) + geom_point(shape=1, size=2, colour="black") + 
  scale_colour_manual(values=col_region) + labs(x="", y="Elevation (m)") 
ggsave("eda/site_region_elDist.pdf", width=2.75, height=5)

plot.sf %>% group_by(BDM, region) %>% st_set_geometry(NULL) %>%
  summarise(el_mn=mean(mnt25), el_med=median(mnt25), 
            el_min=min(mnt25), el_max=max(mnt25)) %>%
  mutate(region=forcats::lvls_reorder(region, 3:1)) %>%
  arrange(region, el_mn) %>% ungroup %>%
  mutate(sortOrder=row_number()) %>%
  ggplot(aes(sortOrder, y=el_mn, ymin=el_min, ymax=el_max, colour=region)) + 
  geom_linerange() + fonts +
  geom_point(size=1.5) + geom_point(shape=1, size=1.5, colour="black") + 
  scale_colour_manual(values=col_region) + labs(x="", y="Elevation (m)") + 
  theme(legend.position="none", panel.border=element_blank(),
        axis.line=element_line(size=0.5, colour="gray30"),
        axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave("eda/site_region_elDistRng.pdf", width=3.75, height=5)

plot.sf %>% mutate(elBin=mnt25 %/% 100 * 100) %>% 
  group_by(elBin) %>% summarise(nPlots=n()) %>%
  ggplot(aes(elBin+50, nPlots)) + geom_bar(stat="identity") + fonts +
  labs(x="Elevation (m)", y="Number of plots") + xlim(300, 3100)
ggsave("eda/plot_elDist.pdf", width=4, height=4)








########------------------------------------------------------------------------
## Environmental Covariates
##

# Plots of environmental covariates
full_join(tibble(el=as.numeric(names(dem_bin)), area=c(dem_bin)),
          plot.sf %>% st_set_geometry(NULL) %>% 
            mutate(el=mnt25 %/% 100 * 100) %>% 
            group_by(el) %>% summarise(nPlots=n()), by="el") %>%
  ggplot(aes(area, nPlots)) + geom_point() +
  labs(x="Area", y="Number of plots") + fonts
ggsave("eda/plot_by_area.pdf", width=4, height=4)
ggplot(tibble(el=as.numeric(names(dem_bin)), area=c(dem_bin)), 
       aes(x=el + 50, y=area)) + geom_bar(stat="identity") + fonts + 
  labs(x="Elevation (m)", y="Area")
ggsave("eda/area_elDist.pdf", width=4, height=4)



site.sf.scaled <- site.sf %>% mutate_at(5:48, ~scale(.)[,1]) %>%
  st_set_geometry(NULL) %>%
  pivot_longer(cols=c(5:21, 25, 27:48, 50), names_to="Var", values_to="ScVal") 

ggplot(site.sf.scaled, aes(el, ScVal, colour=region)) + fonts +
  stat_smooth(method="lm", se=F) + geom_point() + 
  scale_colour_manual(values=col_region) + 
  facet_wrap(~Var, scales="free_y")

site.sf %>% st_set_geometry(NULL) %>% select(el, Tmin, MAT, Tmax, region) %>% 
  pivot_longer(2:4, names_to="Var", values_to="Temperature") %>%
  mutate(Var=factor(Var, levels=c("Tmin", "MAT", "Tmax"))) %>%
  ggplot(aes(el, Temperature, colour=Var, group=paste0(Var, region))) + fonts + 
  geom_point() + stat_smooth(method="lm", se=F) +
  scale_colour_manual(values=c(Tmin="blue", MAT="gray30", Tmax="red")) +
  labs(x="Elevation (m)", y="Temperature (ºC)")


site.sf %>% st_set_geometry(NULL) %>% select(el, DTR, Iso, Tseas, region) %>% 
  pivot_longer(2:4, names_to="Var", values_to="Temperature") %>%
  mutate(Var=factor(Var, levels=c("DTR", "Iso", "Tseas"))) %>%
  ggplot(aes(el, Temperature, colour=Var, group=paste0(Var, region))) + fonts + 
  geom_point() + stat_smooth(method="lm", se=F) +
  labs(x="Elevation (m)", y="Temperature (ºC)")

env.base.gg + 
  geom_point(aes(el, DTR), size=3, colour=blue) +
  geom_point(aes(el, Iso), size=3, colour="gray30") +
  geom_point(aes(el, Tseas), size=3, colour=red) +
  labs(x="Elevation (m)", y="Temperature Range (ºC)")


ggplot(site.sf, aes(x=el)) + 
  geom_point(aes(y=pr_pasture), size=3) + 
  geom_point(aes(y=pr_crop), size=3, colour="brown") + 
  geom_point(aes(y=pr_meadow), size=3, colour="green3") +
  geom_point(aes(y=pr_lawn), size=3, colour="blue") +
  stat_smooth(aes(y=pr_pasture), colour="black", method="loess", se=F, span=2) +
  stat_smooth(aes(y=pr_crop), colour="brown", method="loess", se=F, span=2) +
  stat_smooth(aes(y=pr_meadow), colour="green3", method="loess", se=F, span=2) +
  stat_smooth(aes(y=pr_lawn), colour="blue", method="loess", se=F, span=2) 








########------------------------------------------------------------------------
## Elevational Richness and Abundance
##

# Data summaries
antSum_site <- ant$str %>% st_set_geometry(NULL) %>% 
  filter(TypeOfSample=="soil") %>% group_by(BDM) %>% 
  summarise(nTubes=n(), nSp=n_distinct(SPECIESID), nGen=n_distinct(GENUSID),
            nOcc=n_distinct(Plot_id)) %>%
  left_join(., site.sf, by="BDM")
antSum_plot <- ant$str %>% st_set_geometry(NULL) %>% 
  filter(TypeOfSample=="soil") %>% group_by(Plot_id) %>% 
  summarise(nTubes=n(), nSp=n_distinct(SPECIESID), nGen=n_distinct(GENUSID)) %>%
  left_join(plot.sf, ., by="Plot_id") %>%
  replace_na(list(nTubes=0, nSp=0, nGen=0)) %>%
  mutate(BDM=factor(BDM))


# Site-level patterns
ggplot(antSum_site, aes(el, nTubes/nPlot_Total)) + 
  stat_smooth(method="lm", col="gray30", linetype=2, size=0.75, fill="gray85") +
  geom_point(aes(colour=region), size=2) + 
  scale_colour_manual(values=col_region) + fonts + ylim(0, NA) + 
  labs(x="Elevation (m)", y="Site-level colony density\n(mean tubes per plot)")
ggsave("eda/el_AbundObs_site.pdf", height=5, width=6.5)
summary(lm(I(nTubes/nPlot_Total) ~ el, data=antSum_site))

ggplot(antSum_site, aes(el, nOcc/nPlot_Total)) + 
  stat_smooth(method="lm", col="gray30", linetype=2, size=0.75, fill="gray85") +
  geom_point(aes(colour=region), size=2) + 
  scale_colour_manual(values=col_region) + fonts + ylim(0, NA) + 
  labs(x="Elevation (m)", y="Site-level plot occupancy")
ggsave("eda/el_OccObs_site.pdf", height=5, width=6.5)

ggplot(antSum_site, aes(el, nSp/nPlot_Total)) + 
  stat_smooth(method="lm", col="gray30", linetype=2, size=0.75, fill="gray85") +
  geom_point(aes(colour=region), size=2) + 
  scale_colour_manual(values=col_region) + fonts + ylim(0, NA) + 
  labs(x="Elevation (m)", y="Site-level species density\n(mean species per plot)")
ggsave("eda/el_SpRichObs_site.pdf", height=5, width=6.5)
summary(lm(I(nSp/nPlot_Total) ~ el, data=antSum_site))

ggplot(antSum_site, aes(el, nGen/nPlot_Total)) + 
  stat_smooth(method="lm", col="gray30", linetype=2, size=0.75, fill="gray85") +
  geom_point(aes(colour=region), size=2) + 
  scale_colour_manual(values=col_region) + fonts + ylim(0, NA) + 
  labs(x="Elevation (m)", y="Site-level genus density\n(mean genera per plot)")
ggsave("eda/el_GenRichObs_site.pdf", height=5, width=6.5)
summary(lm(I(nGen/nPlot_Total) ~ el, data=antSum_site))
  
# Plot-level elevational patterns
ggplot(antSum_plot, aes(mnt25, nTubes)) + 
  geom_point(aes(colour=region), alpha=0.5) + 
  stat_smooth(method="lm", col="gray30", linetype=2, size=0.75, fill="gray85") +
  scale_colour_manual(values=col_region) + fonts + ylim(0, NA) + 
  labs(x="Elevation (m)", y="Plot-level colony density\n(tubes per plot)")
ggsave("eda/el_AbundObs_plot.pdf", height=5, width=6.5)
summary(lme4::glmer(nTubes ~ mnt25 + (1|BDM), data=antSum_plot, family="poisson"))

ggplot(antSum_plot, aes(mnt25, nSp)) + 
  geom_point(aes(colour=region), alpha=0.5) + 
  stat_smooth(method="lm", col="gray30", linetype=2, size=0.75, fill="gray85") +
  scale_colour_manual(values=col_region) + fonts + ylim(0, NA) + 
  labs(x="Elevation (m)", y="Plot-level species richness\n(species per plot)")
ggsave("eda/el_SpRichObs_plot.pdf", height=5, width=6.5)
summary(lme4::glmer(nSp ~ mnt25 + (1|BDM), data=antSum_plot, family="poisson"))

ggplot(antSum_plot, aes(mnt25, nGen)) + 
  geom_point(aes(colour=region), alpha=0.5) + 
  stat_smooth(method="lm", col="gray30", linetype=2, size=0.75, fill="gray85") +
  scale_colour_manual(values=col_region) + fonts + ylim(0, NA) + 
  labs(x="Elevation (m)", y="Plot-level genus richness\n(genera per plot)")
ggsave("eda/el_GenRichObs_plot.pdf", height=5, width=6.5)
summary(glm(nGen ~ mnt25, data=antSum_plot, family="poisson"))




ant.elBin <- ant$str %>% st_set_geometry(NULL) %>% 
  filter(TypeOfSample=="soil") %>%
  select(Plot_id, TubeNo, SPECIESID, GENUSID) %>%
  full_join(plot.sf, ., by="Plot_id") %>%
  mutate(elBin=mnt25 %/% 100 * 100) %>%
  group_by(Plot_id, elBin) %>%
  summarise(nSp_plot=n_distinct(SPECIESID, na.rm=T),
            nGen_plot=n_distinct(GENUSID, na.rm=T),
            nTube_plot=sum(!is.na(SPECIESID))) 

ggplot(ant.elBin, aes(x=factor(elBin), y=nTube_plot)) + geom_boxplot()
ggplot(ant.elBin, aes(x=factor(elBin), y=nSp_plot)) + geom_boxplot()
ggplot(ant.elBin, aes(x=factor(elBin), y=nGen_plot)) + geom_boxplot()

ant.elBin %>% mutate(presence=nTube_plot>0) %>%
  group_by(elBin) %>% summarise(prPresence=mean(presence)) %>%
  ggplot(aes(elBin, prPresence)) + geom_point()

ggplot(ant.elBin, aes(nSp_plot, fill=elBin)) + 
  geom_hline(yintercept=0, colour="gray30") + geom_bar(colour="gray30") + 
  facet_wrap(~elBin) + scale_fill_viridis()
ggplot(ant.elBin, aes(nGen_plot, fill=elBin)) + 
  geom_hline(yintercept=0, colour="gray30") + geom_bar(colour="gray30") + 
  facet_wrap(~elBin) + scale_fill_viridis()
ggplot(ant.elBin, aes(nTube_plot, fill=elBin)) + 
  geom_hline(yintercept=0, colour="gray30") + geom_bar(colour="gray30") + 
  facet_wrap(~elBin) + scale_fill_viridis()

ggplot(ant.elBin, aes(nTube_plot, nSp_plot, colour=elBin, group=elBin)) + 
  stat_smooth(method="lm", formula=y~sqrt(x), se=F) + geom_point(alpha=0.5) +
  # facet_wrap(~elBin) + 
  scale_colour_viridis()
ggplot(ant.elBin, aes(nTube_plot, nGen_plot, colour=elBin, group=elBin)) + 
  stat_smooth(method="lm", formula=y~sqrt(x), se=F) + geom_point(alpha=0.5) +
  # facet_wrap(~elBin) + 
  scale_colour_viridis()




ant.elBin %>% group_by(elBin, nSp_plot) %>%
  summarise(ct=n()) %>% 
  group_by(elBin) %>% mutate(prop=ct/sum(ct)) %>%
  ggplot(aes(x=elBin, y=prop, colour=nSp_plot, group=nSp_plot)) + 
  stat_smooth(method="loess", se=F) + 
  geom_point(size=2) + scale_colour_viridis()
ant.elBin %>% group_by(elBin, nGen_plot) %>%
  summarise(ct=n()) %>% 
  group_by(elBin) %>% mutate(prop=ct/sum(ct)) %>%
  ggplot(aes(x=elBin, y=prop, colour=nGen_plot, group=nGen_plot)) + 
  stat_smooth(method="loess", se=F) + 
  geom_point(size=2) + scale_colour_viridis()
ant.elBin %>% group_by(elBin, nTube_plot) %>%
  summarise(ct=n()) %>% 
  group_by(elBin) %>% mutate(prop=ct/sum(ct)) %>%
  ggplot(aes(x=elBin, y=prop, colour=nTube_plot, group=nTube_plot)) + 
  stat_smooth(method="loess", se=F) + 
  geom_point(size=2) + scale_colour_viridis()



ant.elBin.sum <- ant.elBin %>% group_by(elBin) %>% 
  summarise(mnSp=mean(nSp_plot), 
            mnGen=mean(nGen_plot), 
            mnTube=mean(nTube_plot))
ggplot(ant.elBin.sum, aes(x=elBin)) + 
  geom_line(aes(y=mnSp)) + 
  geom_line(aes(y=mnGen), size=1) + 
  geom_line(aes(y=mnTube), colour="red")
ggplot(ant.elBin.sum, aes(x=elBin)) + geom_point(aes(y=mnTube/mnSp))
ggplot(ant.elBin.sum, aes(x=elBin)) + geom_point(aes(y=mnTube/mnGen))




ant$str %>% st_set_geometry(NULL) %>% group_by(SPECIESID) %>%
  summarise(minEl=min(mnt25), 
            maxEl=max(mnt25)) %>%
  mutate(medEl=minEl + (maxEl-minEl)/2,
         assemblage=case_when(maxEl < 1000 ~ 'low',
                              maxEl > 1000 & minEl > 1000 ~ 'high',
                              maxEl > 1000 & minEl < 1000 ~ 'mixed')) %>%
  mutate(assemblage=factor(assemblage, levels=c("low", "mixed", "high"))) %>%
  arrange(assemblage, medEl, desc(minEl)) %>% 
  mutate(spOrd=factor(SPECIESID, levels=SPECIESID)) %>%
  ggplot(aes(x=spOrd, y=medEl, ymin=minEl, ymax=maxEl, colour=assemblage)) + 
  geom_pointrange() + coord_flip()








########------------------------------------------------------------------------
## Ordination
##

##-- distances
site.dist.xy <- units::drop_units(st_distance(site.sf))
site.dist.el <- dist(cbind(site.sf$el, site.sf$el))


##--- community matrices
comm.site <- ant$str %>% filter(TypeOfSample=="soil") %>% st_set_geometry(NULL) %>%
  group_by(BDM, SPECIESID) %>% summarise(nObs=n()) %>%
  ungroup %>% spread(SPECIESID, nObs)
comm.site.mx <- comm.site %>% select(-BDM) %>% as.matrix
comm.site.mx[is.na(comm.site.mx)] <- 0
rownames(comm.site.mx) <- comm.site$BDM

comm.plot <- ant$str %>% filter(TypeOfSample=="soil") %>% st_set_geometry(NULL) %>%
  filter(SPECIESID!="Tapi_erra") %>%
  group_by(Plot_id, GENUSID) %>% summarise(nObs=n()) %>%
  ungroup %>% spread(GENUSID, nObs)
comm.plot.mx <- comm.plot %>% select(-Plot_id) %>% as.matrix
comm.plot.mx[is.na(comm.plot.mx)] <- 0
rownames(comm.plot.mx) <- comm.plot$Plot_id 
comm.plot.mx <- comm.plot.mx[rowSums(comm.plot.mx)>1,]

comm.bin <- ant$str %>% filter(TypeOfSample=="soil") %>% st_set_geometry(NULL) %>%
  mutate(elBin=mnt25 %/% 100 * 100) %>%
  group_by(elBin, SPECIESID) %>% summarise(nObs=n()) %>%
  ungroup %>% spread(SPECIESID, nObs)
comm.bin.mx <- comm.bin %>% select(-elBin) %>% as.matrix
comm.bin.mx[is.na(comm.bin.mx)] <- 0
rownames(comm.bin.mx) <- comm.bin$elBin
comm.bin.mx <- comm.bin.mx[rowSums(comm.bin.mx)>1,]




##--- NMDS
nmds.site <- vegan::metaMDS(comm.site.mx, trymax=200, distance="bray", k=2)
nmds.plot <- vegan::metaMDS(comm.plot.mx, trymax=200, distance="jaccard", k=3)
nmds.bin <- vegan::metaMDS(comm.bin.mx, trymax=200, distance="bray", k=2)

nmds.site.df <- as.data.frame(nmds.site$points) %>%
  mutate(BDM=as.numeric(rownames(.))) %>%
  full_join(., st_set_geometry(site.sf, NULL), by="BDM") %>%
  mutate(el_cat=forcats::lvls_revalue(region, c("Mtn", "Mtn", "Low")))
nmds.plot.df <- as.data.frame(nmds.plot$points) %>%
  mutate(Plot_id=rownames(.)) %>%
  full_join(., st_set_geometry(plot.sf, NULL), by="Plot_id") %>%
  mutate(el_cat=forcats::lvls_revalue(region, c("Mtn", "Mtn", "Low")))
nmds.bin.df <- as.data.frame(nmds.bin$points) %>%
  mutate(elBin=rownames(.))



##--- bioenv
mod_v <- paste(names(nmds.site.df)[c(4, 6, 7, 
                                     # 9, 11, 13, 15, 17, 19, 
                                     # 21:23, 27, 29, 
                                     # 32:50,
                                     32, 33, 34, 35, 36, 38, 40, 42, 46, 47, 48,
                                     51:55,
                                     66)], collapse=" + ")
bioenv.site <- bioenv(as.formula(paste("vegdist(comm.site.mx, method='bray') ~", 
                                       mod_v, collapse="")), 
                      upto=10, data=nmds.site.df, parallel=4, metric="gower")
bioenv.plot <- bioenv(as.formula(paste("vegdist(comm.plot.mx, method='bray') ~", 
                                       mod_v, collapse="")), 
                      upto=10, data=nmds.plot.df, parallel=4)



##--- mantel tests
mantel.site <- mantel(vegdist(comm.site.mx, method='bray'), bioenvdist(bioenv.site))
mantel.plot <- mantel(vegdist(comm.plot.mx, method='bray'), bioenvdist(bioenv.plot))



##--- adonis or anosim
adonis(vegdist(comm.site.mx, method='bray') ~ el_cat + Iso + Tmax + Tseas + el, 
       data=nmds.site.df)
anosim(vegdist(comm.site.mx, method='bray'), grouping=nmds.site.df$el_cat)
anosim(vegdist(comm.site.mx, method='bray'), grouping=nmds.site.df$region)



##--- cluster analysis: K-means
# Sites
kmns.site.ls <- map(1:5, ~kmeans(as.matrix(vegdist(comm.site.mx)), ., 1e4, 1e2))
plot(1:5, map_dbl(kmns.site.ls, ~.$tot.withinss), 
     xlab="Clusters", ylab="SS", type="b")
nmds.site.df <- nmds.site.df %>%
  mutate(kmns_2=factor(kmns.site.ls[[2]]$cluster))

# Plots
kmns.plot.ls <- map(1:5, ~kmeans(as.matrix(vegdist(comm.plot.mx)), ., 1e4, 1e2))
plot(1:5, map_dbl(kmns.plot.ls, ~.$tot.withinss), 
     xlab="Clusters", ylab="SS", type="b")
nmds.plot.df <-  left_join(nmds.plot.df, 
                           tibble(kmns_2=factor(kmns.plot.ls[[2]]$cluster), 
                                  Plot_id=names(kmns.plot.ls[[2]]$cluster)),
                           by="Plot_id")

# Bins
kmns.bin.ls <- map(1:5, ~kmeans(as.matrix(vegdist(comm.bin.mx)), ., 1e4, 1e2))
plot(1:5, map_dbl(kmns.bin.ls, ~.$tot.withinss), 
     xlab="Clusters", ylab="SS", type="b")
nmds.bin.df <- nmds.bin.df %>%
  mutate(kmns_2=factor(kmns.bin.ls[[2]]$cluster))



##--- plots
# Sites
ggplot(nmds.site.df, aes(x=MDS1, y=MDS2, colour=region, fill=region)) + fonts +
  ggConvexHull::geom_convexhull(alpha=0.2, size=0.25) + geom_point() + 
  geom_point(data=nmds.site.df %>% group_by(region) %>% 
               summarise(MDS1=mean(MDS1), MDS2=mean(MDS2)), shape=1, size=4) +
  scale_colour_manual(values=col_region) + scale_fill_manual(values=col_region) 
ggsave("eda/NMDS_site_region.pdf", width=5, height=3.5)
full_join(select(site.sf, BDM), nmds.site.df, by="BDM") %>%
  ggplot() + fonts + geom_sf(data=VD) + 
  theme(legend.position="none", axis.text=element_blank()) +
  geom_sf(aes(colour=region, fill=region), size=2) + 
  scale_colour_manual(values=col_region) + scale_fill_manual(values=col_region)
ggsave("eda/NMDS_site_region_MAP.pdf", width=5, height=5)

ggplot(nmds.site.df, aes(x=MDS1, y=MDS2, colour=el_cat, fill=el_cat)) + fonts +
  ggConvexHull::geom_convexhull(alpha=0.2, size=0.25) + geom_point() + 
  geom_point(data=nmds.site.df %>% group_by(el_cat) %>% 
               summarise(MDS1=mean(MDS1), MDS2=mean(MDS2)), shape=1, size=4) +
  scale_colour_manual("", values=col_mtn) + 
  scale_fill_manual("", values=col_mtn)
ggsave("eda/NMDS_site_mtn.pdf", width=5, height=3.5)
full_join(select(site.sf, BDM), nmds.site.df, by="BDM") %>%
  ggplot() + fonts + geom_sf(data=VD) + 
  theme(legend.position="none", axis.text=element_blank()) +
  geom_sf(aes(colour=el_cat, fill=el_cat), size=2) + 
  scale_colour_manual(values=col_mtn) + scale_fill_manual(values=col_mtn)
ggsave("eda/NMDS_site_mtn_MAP.pdf", width=5, height=5)

ggplot(nmds.site.df, aes(x=MDS1, y=MDS2, colour=kmns_2, fill=kmns_2)) + fonts +
  ggConvexHull::geom_convexhull(alpha=0.2, size=0.25) + geom_point() + 
  geom_point(data=nmds.site.df %>% group_by(kmns_2) %>% 
               summarise(MDS1=mean(MDS1), MDS2=mean(MDS2)), shape=1, size=4) +
  scale_colour_manual("Clusters", values=setNames(col_mtn, 1:2)) + 
  scale_fill_manual("Clusters", values=setNames(col_mtn, 1:2))
ggsave("eda/NMDS_site_kmeans2.pdf", width=5, height=3.5)
full_join(select(site.sf, BDM), nmds.site.df, by="BDM") %>%
  ggplot() + fonts + geom_sf(data=VD) + 
  theme(legend.position="none", axis.text=element_blank()) +
  geom_sf(aes(colour=kmns_2, fill=kmns_2), size=2) + 
  scale_colour_manual("Clusters", values=setNames(col_mtn, 1:2)) + 
  scale_fill_manual("Clusters", values=setNames(col_mtn, 1:2)) 
ggsave("eda/NMDS_site_kmeans_MAP.pdf", width=5, height=5)

# Plots
ggplot(nmds.plot.df, aes(x=MDS1, y=MDS2, colour=region, fill=region)) + 
  ggConvexHull::geom_convexhull(alpha=0.2, size=0.25) + geom_point() + 
  geom_point(data=nmds.plot.df %>% group_by(region) %>%
               summarise(MDS1=mean(MDS1, na.rm=T), MDS2=mean(MDS2, na.rm=T)),
             shape=1, size=4) +
  scale_colour_manual(values=col_region) +
  scale_fill_manual(values=col_region) +
  theme(panel.grid=element_blank())
ggsave("eda/NMDS_plot_region.pdf", width=5, height=3.5)
ggplot(nmds.plot.df, aes(x=MDS1, y=MDS2, colour=el_cat, fill=el_cat)) + 
  ggConvexHull::geom_convexhull(alpha=0.2, size=0.25) + geom_point() + 
  geom_point(data=nmds.plot.df %>% group_by(el_cat) %>% 
               summarise(MDS1=mean(MDS1, na.rm=T), MDS2=mean(MDS2, na.rm=T)), 
             shape=1, size=4) +
  scale_colour_manual(values=col_mtn) + scale_fill_manual(values=col_mtn) +
  theme(panel.grid=element_blank())
ggsave("eda/NMDS_plot_mtn.pdf", width=5, height=3.5)
ggplot(nmds.plot.df, aes(x=MDS1, y=MDS2, colour=kmns_2, fill=kmns_2)) + 
  ggConvexHull::geom_convexhull(alpha=0.2, size=0.25) + geom_point() + 
  geom_point(data=nmds.plot.df %>% group_by(kmns_2) %>% 
               summarise(MDS1=mean(MDS1, na.rm=T), MDS2=mean(MDS2, na.rm=T)), 
             shape=1, size=4) +
  scale_colour_manual("Clusters", values=setNames(col_mtn, 2:1)) + 
  scale_fill_manual("Clusters", values=setNames(col_mtn, 2:1)) + 
  theme(panel.grid=element_blank())
ggsave("eda/NMDS_plot_kmeans2.pdf", width=5, height=3.5)

# 3d 
car::scatter3d(nmds.plot.df$MDS1, nmds.plot.df$MDS2, nmds.plot.df$MDS3, 
               point.col=col_region[as.numeric(nmds.plot.df$region)], surface=F)
car::scatter3d(nmds.plot.df$MDS1, nmds.plot.df$MDS2, nmds.plot.df$MDS3, 
               point.col=col_mtn[as.numeric(nmds.plot.df$el_cat)], surface=F)
car::scatter3d(nmds.plot.df$MDS1, nmds.plot.df$MDS2, nmds.plot.df$MDS3, 
               point.col=col_mtn[as.numeric(nmds.plot.df$kmns_2)], surface=F)


ggplot(st_as_sf(nmds.site.df), aes(fill=kmns_2)) + geom_sf()
ggplot(st_as_sf(filter(nmds.plot.df, !is.na(kmns_2))), aes(colour=kmns_2)) + geom_sf()


nmds.plot.full.df <- na.omit(nmds.plot.df %>% select(-TypeOfOpen))
plot.glm.full <- glm(as.formula(paste0("kmns_2 ~ ", 
                                       paste0(names(nmds.plot.full.df)[c(5,
                                                                         7:10,
                                                                         13:15,
                                                                         20,
                                                                         26:30,
                                                                         36:42)],
                                              collapse=" + "))), 
                     family="binomial", na.action="na.fail", 
                     data=nmds.plot.full.df)
plot.dredge <- MuMIn::dredge(plot.glm.full)
plot.dredge %>% mutate_if(is.factor, as.numeric) %>% 
  replace_na(., map_dfc(names(plot.dredge), ~setNames(data.frame(a=0), .x))) %>%
  pivot_longer(., cols=1:(ncol(.)-5), names_to="Variable", values_to="Val") %>%
  group_by(Variable) %>%
  summarise(wt=sum(weight*(Val != 0))) %>% arrange(desc(wt)) %>%
  mutate(wtOrd=factor(Variable, levels=unique(Variable))) %>%
  ggplot(aes(wtOrd, wt)) + geom_bar(stat="identity")



plot.kmns2.glm <- nmds.plot.df %>% filter(!is.na(kmns_2)) %>% 
  mutate(kmns_2=as.numeric(kmns_2==1)) %>%
  glm(kmns_2 ~ mnt25, data=., family="binomial")
summary(plot.kmns2.glm)
nmds.plot.df %>% filter(!is.na(kmns_2)) %>% 
  mutate(kmns_2=as.numeric(kmns_2==1),
         pred=predict(plot.kmns2.glm, type="response")) %>%
  arrange(mnt25) %>%
  ggplot(aes(mnt25, kmns_2)) + geom_point(shape=1) + geom_line(aes(y=pred))

ggplot(filter(nmds.plot.df, !is.na(kmns_2)), aes(x=mnt25, y=as.numeric(kmns_2==1))) + 
  geom_point(alpha=0.75) + stat_smooth() 

ggplot(filter(nmds.site.df, !is.na(kmns_2)), aes(x=el, y=as.numeric(kmns_2==2))) + 
  geom_point(alpha=0.75)






########------------------------------------------------------------------------
## Beta diversity
##

beta.site <- append(beta.pair.abund(comm.site.mx), 
                    beta.pair((comm.site.mx>0)*1, index.family="jaccard"))

beta.df <- tibble(BDM=rep(attr(beta.site$beta.bray, "Label"), each=44),
                  BDM2=rep(attr(beta.site$beta.bray, "Label"), times=44)) %>%
  cbind(map_df(beta.site, ~c(as.matrix(.)))) %>%
  left_join(., site.sf %>% mutate(BDM=as.character(BDM)), by="BDM") %>%
  left_join(., site.sf %>% mutate(BDM=as.character(BDM)) %>% rename(BDM2=BDM), 
            by="BDM2", suffix=c("", "2")) %>%
  mutate(elBin=el %/% 100 * 100, elBin2=el2 %/% 100 * 100,
         el_cat=forcats::lvls_revalue(region, c("Mtn", "Mtn", "Low")),
         el_cat2=forcats::lvls_revalue(region2, c("Mtn", "Mtn", "Low")),
         distance=units::drop_units(c(st_distance(arrange(site.sf, BDM)))),
         bioenv=c(as.matrix(bioenvdist(bioenv.site)))) %>% 
  filter(el != el2)



beta.bin <- append(beta.pair.abund(comm.bin.mx), 
                   beta.pair((comm.bin.mx>0)*1, index.family="jaccard"))

beta.bin.df <- tibble(elBin=rep(attr(beta.bin$beta.bray, "Label"), each=19),
                      elBin2=rep(attr(beta.bin$beta.bray, "Label"), times=19)) %>%
  cbind(map_df(beta.bin, ~c(as.matrix(.)))) %>%
  mutate(elBin=as.numeric(elBin), elBin2=as.numeric(elBin2)) %>%
  filter(elBin != elBin2)
beta.bin.adj <- beta.bin.df %>% filter(elBin2 == elBin + 100)
ggplot(beta.bin.adj, aes(x=elBin, y=beta.bray.bal)) + geom_line() 


ggplot(beta.bin.df, aes(x=elBin, y=beta.bray, colour=elBin2, group=elBin2)) + ylim(0,NA) +
  stat_smooth(method="loess", se=F, size=0.5) + geom_point() + 
  scale_colour_viridis("Elevation (m)\nof site 2") + 
  labs(x="Elevation (m) of site 1", y="Bray-Curtis Dissimilarity")
ggplot(beta.bin.df, aes(x=elBin, y=beta.bray)) + 
  geom_vline(aes(xintercept=elBin2), linetype=2) +
  geom_point(alpha=0.5) + 
  stat_smooth(se=F, method="loess", size=0.5) +
  facet_wrap(~elBin2) + ylim(0,NA) + theme(panel.grid=element_blank()) +
  labs(x="Elevation (m) of site 1", y="Bray-Curtis Dissimilarity") 

ggplot(beta.df, aes(x=el, y=el2, colour=beta.bray)) + geom_point() + 
  scale_colour_viridis()
ggplot(beta.df, aes(x=el, y=el2, colour=beta.jac)) + geom_point() + 
  scale_colour_viridis()

ggplot(beta.df, aes(x=el, y=beta.bray, colour=el2, group=el2)) + ylim(0,NA) +
  stat_smooth(method="loess", span=2, se=F, size=0.5) + geom_point() + 
  scale_colour_viridis("Elevation (m)\nof site 2") + 
  labs(x="Elevation (m) of site 1", y="Bray-Curtis Dissimilarity")
ggsave("eda/beta_bray_el-el.pdf", width=6, height=4.5)
ggplot(beta.df, aes(x=el, y=beta.jac, colour=el2, group=el2)) + ylim(0,NA) +
  stat_smooth(method="loess", span=2, se=F, size=0.75) + geom_point() + 
  scale_colour_viridis("Elevation (m)\nof site 2") + 
  labs(x="Elevation (m) of site 1", y="Jaccard Dissimilarity")


ggplot(beta.df, aes(x=el, y=beta.bray, colour=el_cat2, group=el2)) + ylim(0,NA) +
  stat_smooth(method="loess", span=2, se=F, size=0.5) + geom_point() + 
  scale_colour_manual("Site 2", values=col_mtn) + 
  labs(x="Elevation (m) of site 1", y="Bray-Curtis Dissimilarity")
ggsave("eda/beta_bray_el-mtn.pdf", width=6, height=4.5)
ggplot(beta.df, aes(x=el, y=beta.bray, colour=region2, group=el2)) + ylim(0,NA) +
  stat_smooth(method="loess", span=2, se=F, size=0.5) + geom_point() + 
  scale_colour_manual("Site 2", values=col_region) + 
  labs(x="Elevation (m) of site 1", y="Bray-Curtis Dissimilarity")
ggsave("eda/beta_bray_el-reg.pdf", width=6, height=4.5)
ggplot(beta.df, aes(x=el, y=beta.jac, colour=region2, group=el2)) + ylim(0,NA) +
  stat_smooth(method="loess", span=2, se=F, size=0.75) + geom_point() + 
  scale_colour_manual("Site 2", values=col_region) + 
  labs(x="Elevation (m) of site 1", y="Jaccard Dissimilarity")



# overall
ggplot(beta.df, aes(x=el, y=beta.bray, colour=region)) + 
  geom_vline(aes(xintercept=elBin2), linetype=2) +
  geom_point(alpha=0.5) + 
  stat_smooth(se=F, method="lm", size=0.75) +
  facet_wrap(~round(elBin2, 1)) + ylim(0,NA) + theme(panel.grid=element_blank()) +
  scale_colour_manual("Site 2", values=col_region) + 
  labs(x="Elevation (m) of site 1", y="Bray-Curtis Dissimilarity") 
ggsave("eda/beta_bray_el-region_byEl2.pdf", width=7, height=6)
ggplot(beta.df, aes(x=el, y=beta.jac, colour=region)) + 
  geom_vline(aes(xintercept=elBin2), linetype=2) +
  geom_point(alpha=0.5) + 
  stat_smooth(se=F, method="lm", size=0.75) +
  facet_wrap(~elBin2) + ylim(0,NA) + 
  scale_colour_manual("Site 2", values=col_region) + 
  labs(x="Elevation (m) of site 1", y="Jaccard Dissimilarity")

# turnover
ggplot(beta.df, aes(x=el, y=beta.bray.bal, colour=region)) + 
  geom_vline(aes(xintercept=elBin2), linetype=2) +
  geom_point(alpha=0.5) + 
  stat_smooth(se=F, method="lm", size=0.75) +
  facet_wrap(~elBin2) + ylim(0,NA) + theme(panel.grid=element_blank()) +
  scale_colour_manual("Site 2", values=col_region) + 
  labs(x="Elevation (m) of site 1", y="Bray-Curtis Dissimilarity (species replacement)")
ggsave("eda/beta_brayBal_el-region_byEl2.pdf", width=7, height=6)
ggplot(beta.df, aes(x=elBin, y=beta.jtu, colour=region)) + 
  geom_vline(aes(xintercept=elBin2), linetype=2) +
  geom_point(alpha=0.5) + 
  stat_smooth(se=F, method="lm", size=0.75) +
  facet_wrap(~elBin2) + ylim(0,NA) + 
  scale_colour_manual("Site 2", values=col_region) + 
  labs(x="Elevation (m) of site 1", y="Jaccard Dissimilarity (species replacement)")

# nestedness
ggplot(beta.df, aes(x=el, y=beta.bray.gra, colour=region)) + 
  geom_vline(aes(xintercept=elBin2), linetype=2) +
  geom_point(alpha=0.5) + 
  stat_smooth(se=F, method="lm", size=0.75) +
  facet_wrap(~elBin2) + ylim(0,NA) +
  scale_colour_manual("Site 2", values=col_region) + 
  labs(x="Elevation (m) of site 1", y="Bray-Curtis Dissimilarity (loss of species)")
ggplot(beta.df, aes(x=el, y=beta.jne, colour=region)) + 
  geom_vline(aes(xintercept=elBin2), linetype=2) +
  geom_point(alpha=0.5) + 
  stat_smooth(se=F, method="lm", size=0.75) +
  facet_wrap(~elBin2) + ylim(0,NA) + 
  scale_colour_manual("Site 2", values=col_region) + 
  labs(x="Elevation (m) of site 1", y="Jaccard Dissimilarity (loss of species)")


ggplot(beta.df, aes(x=factor(elBin), y=beta.bray, fill=region2)) + 
  geom_boxplot(width=0.5, size=0.5, outlier.size=0.5) + facet_grid(region~.) +
  scale_fill_manual("Site 2", values=col_region)

ggplot(beta.df, aes(x=factor(elBin), y=beta.bray, fill=region2)) + 
  geom_boxplot(width=0.5, size=0.5, outlier.size=0.5) + facet_wrap(~elBin2) +
  scale_fill_manual("Site 2", values=col_region)

ggplot(beta.df, aes(x=el, y=beta.bray, colour=region2)) + 
  geom_point(shape=1, alpha=0.5) + 
  scale_colour_manual("Site 2", values=col_region) + 
  stat_smooth(method="lm") + facet_grid(region~.) 


ggplot(beta.df, aes(x=abs(el-el2), y=beta.bray, colour=region2)) + 
  geom_point(alpha=0.5) + scale_colour_manual("Site 2", values=col_region) + 
  stat_smooth(method="lm", se=F) + 
  facet_grid(.~region)
ggplot(beta.df, aes(x=distance, y=beta.bray, colour=region2)) + 
  geom_point(alpha=0.5) + scale_colour_manual("Site 2", values=col_region) + 
  stat_smooth(method="lm", se=F) + 
  facet_grid(.~region)
ggplot(beta.df, aes(x=bioenv, y=beta.bray, colour=region2)) + 
  geom_point(alpha=0.5) + scale_colour_manual("Site 2", values=col_region) + 
  stat_smooth(method="lm", se=F) + 
  facet_grid(.~region)
ggplot(beta.df, aes(x=abs(Tmax2-Tmax), y=beta.bray, colour=region2)) + 
  geom_point(alpha=0.5) + scale_colour_manual("Site 2", values=col_region) + 
  stat_smooth(method="lm", se=F) + 
  facet_grid(.~region)


summary(lm(beta.bray ~ bioenv, data=beta.df))
summary(lm(beta.bray ~ distance, data=beta.df))
summary(lm(beta.bray ~ abs(el-el2), data=beta.df))

summary(lme4::lmer(beta.bray ~ bioenv + (1|el_cat2), data=beta.df))


beta.df %>% st_as_sf() %>%
  ggplot(aes(colour=beta.bray, fill=beta.bray)) + geom_sf(size=3) + 
  scale_colour_viridis() + scale_fill_viridis() + facet_wrap(~BDM2, ncol=8)
beta.df %>% st_as_sf() %>%
  ggplot(aes(colour=distance, fill=distance)) + geom_sf(size=3) + 
  scale_colour_viridis() + scale_fill_viridis() + facet_wrap(~BDM2, ncol=8)
beta.df %>% st_as_sf() %>%
  ggplot(aes(colour=abs(el-el2), fill=abs(el-el2))) + geom_sf(size=3) + 
  scale_colour_viridis() + scale_fill_viridis() + facet_wrap(~BDM2, ncol=9)

ggplot(filter(beta.df, el!=el2), aes(x=factor(el), y=beta.bray, fill=region)) + 
  geom_boxplot() + scale_fill_manual("Site 2", values=col_region) + 
  facet_grid(region2~.)


summary(lm(beta.bray ~ distance * region, data=beta.df))


segmod <- segmented(lm(beta.bray ~ elBin * elBin2, 
                       data=beta.df), seg.Z=~elBin)
plot(segmod)

segmod <- segmented(lm(beta.bray ~ el*as.factor(el2), 
                       data=filter(beta.df, el2 > 1000)), seg.Z=~el)
plot(segmod)



plot(site.dist.xy, as.matrix(beta.site$beta.bray.bal))
plot(site.dist.el, beta.site$beta.bray.bal)
plot(bioenvdist(bioenv.site), beta.site$beta.bray.bal)
plot(vegdist(comm.site.mx, method='bray'), beta.site$beta.bray.bal)

plot(site.dist.xy, as.matrix(bioenvdist(bioenv.site)))












rownames(site.dist) <- colnames(site.dist) <- site.sf$BDM
heatmap(site.dist)
heatmap(as.matrix(vegdist(comm.site.mx, method="bray")))
heatmap(as.matrix(bioenvdist(nmds.mod)))

plot(site.dist, as.matrix(vegdist(comm.site.mx, method='bray')))
plot(bioenvdist(bioenv.site), vegdist(comm.site.mx, method='bray'))
plot(site.dist, as.matrix(bioenvdist(bioenv.site)))

mantel(as.matrix(vegdist(comm.site.mx, method='bray')), site.dist)



plot.dist <- units::drop_units(st_distance(plot.sf))
rownames(plot.dist) <- colnames(plot.dist) <- plot.sf$Plot_id
heatmap(plot.dist)











# stressplot(str.nmds)
ordiplot(nmds.site, type="n")
orditorp(nmds.site, display="species", col="red", air=0.01)
orditorp(nmds.site, display="sites", cex=1.25, air=0.01)



region <- nmds.site.df$region
colors <- c("gray", "red", "blue")[as.numeric(as.factor(nmds.site.df$region))]
ordisurf(nmds.site, nmds.site.df$el, col="gray30")
ordiplot(nmds.site, display="sites", type="p", cex=1.5)
for(i in unique(region)) {
  # ordiellipse(nmds.site$point[grep(i,region),], draw="polygon",
           # groups=region[region==i], col=colors[grep(i,region)],label=T) 
  ordihull(nmds.site$point[grep(i,region),], draw="polygon", alpha=0.2,
              groups=region[region==i], col=colors[grep(i,region)],label=T) 
}

ordisurf(nmds.site, nmds.df$MAT/10, col="forestgreen")
ordiellipse(nmds.site, groups=nmds.df$region, draw="polygon", col="grey90", label=T)


cor.nmds <- round(cor(as.matrix(select_if(nmds.df, is.numeric)))[,1:2], 3)
cor.nmds[order(abs(cor.nmds[,1]), decreasing=T),][1:10,]
cor.nmds[order(abs(cor.nmds[,2]), decreasing=T),][1:10,]


ggplot(nmds.site.df, aes(MDS1, MDS2, colour=Tmax/10, shape=region)) + 
  geom_point(size=3) + scale_colour_viridis() 
ggplot(nmds.site.df, aes(MDS1, MDS2, colour=region)) + geom_point(size=3) + 
  scale_colour_manual(values=c("gray30", "red", "blue"))
ggplot(nmds.df, aes(Tmax/10, MDS1, colour=region)) + geom_point(size=3)
ggplot(nmds.df, aes(area_Open, MDS2, colour=region)) + 
  geom_point(size=3)
ggplot(nmds.df, aes(Moss_mean, MDS1, colour=region)) + 
  geom_point(size=3)













########------------------------------------------------------------------------
## Rarefaction and Extrapolation
##

S <- specnumber(comm.site.mx)
rar.S <- rarefy(comm.site.mx, sample=seq(6, 51, by=5))
rar.df <- as.data.frame(cbind(S, rar.S)) %>%
  mutate(BDM=as.numeric(rownames(.))) %>%
  full_join(., site.sf, by="BDM") %>%
  pivot_longer(cols=2:11, names_to="SampleSize", values_to="Richness")
ggplot(rar.df, aes(x=S, y=Richness, colour=SampleSize)) + 
  geom_point() + stat_smooth(se=F, method="loess", span=2) + 
  scale_colour_viridis_d()
ggplot(rar.df, aes(x=el, y=S)) + 
  geom_point() + stat_smooth(se=F, method="loess", span=2) 






chao.out <- map_df(1:nrow(comm.site.mx), 
                   ~ChaoRichness(comm.site.mx[.,], "abundance")) %>%
  mutate(BDM=as.numeric(rownames(comm.site.mx))) %>%
  full_join(., site.sf, by="BDM")
ggplot(chao.out, aes(x=el, y=Estimator)) +
  geom_point() + stat_smooth(method="loess", se=F, span=2)
ggplot(chao.out, aes(x=area_Open, y=Estimator, colour=Observed)) +
  stat_smooth(method="loess", se=F, span=2, colour="gray30", linetype=2) + 
  geom_point() + 
  scale_colour_viridis() + facet_wrap(~region)


sizes <- seq(1, 100)
inext.out <- iNEXT(t(comm.site.mx), size=sizes)
next.df <- do.call('rbind', inext.out$iNextEst) %>% 
  as.data.frame %>%
  mutate(BDM=as.numeric(str_split_fixed(rownames(.), "\\.", 2)[,1])) %>%
  full_join(., site.sf, by="BDM")
divEst.df <- inext.out$AsyEst %>%
  mutate(BDM=as.numeric(rep(rownames(comm.site.mx), each=3))) %>%
  full_join(., site.sf, by="BDM")


ggplot(next.df, aes(m, qD, colour=el, fill=el, group=el)) + 
  geom_line() + scale_colour_viridis() + facet_wrap(~region)

ggplot(filter(next.df, m==100), aes(el, qD)) + geom_point()

ggplot(divEst.df, aes(x=el, y=Estimator, ymin=Estimator-s.e., 
                      ymax=Estimator+s.e., colour=Diversity)) + 
  geom_pointrange(fatten=0.8) + facet_wrap(~region)
ggplot(divEst.df, aes(x=el, y=Observed)) + geom_point() +
  stat_smooth(method="lm", formula=y~x+I(x^2), se=F)
ggplot(divEst.df, aes(x=el, y=Estimator, colour=Diversity)) + geom_point() +
  stat_smooth(method="lm", formula=y~x+I(x^2), se=F)
ggplot(divEst.df, aes(x=el, y=Estimator, colour=Diversity)) + geom_point() +
  stat_smooth(method="loess", se=F)





umap.out <- umap(comm.site.mx)
umap.df <- as.data.frame(umap.out$layout) %>%
  mutate(BDM=as.numeric(rownames(.))) %>%
  full_join(., site.sf, by="BDM")
ggplot(umap.df, aes(V1, V2, colour=Tmax/10, shape=region)) + 
  geom_point(size=3) + scale_colour_viridis() 
cor.umap <- round(cor(as.matrix(select_if(umap.df, is.numeric)))[,1:2], 3)
cor.umap[order(abs(cor.umap[,1]), decreasing=T),][1:10,]
cor.umap[order(abs(cor.umap[,2]), decreasing=T),][1:10,]





site.cov <- st_set_geometry(site.sf, NULL) %>% 
  select(-BDM, -region, -BDM_id, -Sample_date, -SoilTemp_mean, -SoilTemp_sd)
PCA <- rda(site.cov, scale=T)
barplot(as.vector(PCA$CA$eig)/sum(PCA$CA$eig))
sum((as.vector(PCA$CA$eig)/sum(PCA$CA$eig))[1:2])
plot(PCA)
plot(PCA, display = "sites", type = "text")
plot(PCA, display = "species", type = "text")

# In a biplot of a PCA, species' scores are drawn as arrows 
# that point in the direction of increasing values for that variable
biplot(PCA) # biplot of axis 1 vs 2


pca.df <- st_set_geometry(site.sf, NULL) %>% 
  mutate(PCA1=summary(PCA)$sites[,1],
         PCA2=summary(PCA)$sites[,2],
         PCA3=summary(PCA)$sites[,3])
cor.pca <- round(cor(as.matrix(select_if(pca.df, is.numeric)))[,49:51], 3)
cor.pca[order(abs(cor.pca[,1]), decreasing=T),][1:10,]
cor.pca[order(abs(cor.pca[,2]), decreasing=T),][1:10,]
cor.pca[order(abs(cor.pca[,3]), decreasing=T),][1:10,]
  
ggplot(pca.df, aes(PCA1, PCA2, colour=Tmax, shape=region)) + 
  geom_point(size=3) + scale_colour_viridis() 








#--- libraries and options



#--- site, land cover, & plot info
site_summary <- read_csv("data/opfo_siteSummaryProcessed.csv")
lc_summary <- read_csv("data/opfo_lcSummaryProcessed.csv")
plot_i <- read_csv("data/opfo_envDataProcessed.csv")
lc_i <- readxl::read_xlsx("data/landcover_id.xlsx", 1)
site.lc.prop <- readxl::read_xlsx("data/site_summaries_25.xlsx", 1) %>%
  mutate(BDM_id=str_pad(BDM_id, 2, side="left", pad="0"),
         Sample_Date=lubridate::ymd(Sample_Date)) %>%
  mutate(Canopy=lc_i$Canopy[match(Categorie, lc_i$LC)])
site_summary <- site_summary %>%
  left_join(., site.lc.prop %>% group_by(BDM) %>% 
              summarise(tr=sum(trans_len),
                        area=sum(pctArea)), by="BDM")
sites.sf <- st_read(paste0(gis.dir, "Fourmis_data_shp/MilieuxBDMdiss.shp")) %>%
  st_transform(3857)
site.box.sf <- st_read(paste0(gis.dir, "Fourmis_data_shp/BDMplus.shp")) %>%
  st_transform(st_crs(sites.sf)) %>% 
  mutate(BDM=Id, 
         n_plots=site_summary$n_plots[match(BDM, site_summary$BDM)],
         tran_len=site_summary$tr[match(BDM, site_summary$BDM)],
         area=site_summary$area[match(BDM, site_summary$BDM)])


#--- environmental data
VD_raw <- st_read(paste0(gis.dir, "Vaud_boundaries.shp")) 
VD <- st_union(VD_raw)
VD_ext <- raster::extent(matrix(st_bbox(VD), ncol=2))
grd <- raster::raster(ext=VD_ext, crs=st_crs(VD), resolution=1000) %>%
  raster::rasterize(VD_ext, .) %>% 
  raster::rasterToPolygons(., n=4)
grd@data$layer <- 1:raster::ncell(grd)
grd.sf <- st_as_sf(grd) %>% st_set_crs(st_crs(VD)) %>% rename(id=layer) %>%
  mutate(inbd=c(st_intersects(., VD, sparse=F)))



dem <- raster::raster(paste0(gis.dir, "aster/CH_21781.tif")) %>%
  raster::crop(., VD_raw) 
slope <- raster::raster(paste0(gis.dir, "aster/CH_slope_21781.tif")) %>%
  raster::crop(., VD_raw) 



clim_1km.df <- map_dfc(list(MAT="CHELSA_bio10_01.tif",
                            MAP="CHELSA_bio10_12.tif",
                            TAR="CHELSA_bio10_07.tif"), 
                       ~raster::raster(paste0(gis.dir, "chelsa/", .x)) %>%
                         raster::projectRaster(., dem) %>%
                         raster::crop(., VD_raw) %>%
                         raster::extract(., grd.sf, fun=mean, na.rm=T))



#--- ant data
ant.ss <- gs_read(ss=gs_title("Scientific Ant Inventory VD"), ws="Transfer") %>%
  filter(!is.na(lon)) %>% rename(TubeNo=CATALOGUENUMBER, Plot_id=PLOT) %>%
  mutate(Plot_id=as.character(Plot_id),
         SampleDate=lubridate::dmy(SampleDate)) %>%
  # left_join(., select(plot.sf, -c(3,21:30)), by="Plot_id") %>%
  mutate(Canopy=lc_i$Canopy[match(Categorie, lc_i$LC)]) %>%
  st_as_sf(coords=c("lon", "lat")) %>% st_set_crs(4326) %>% st_transform(21781)
ant.ss.boxTot <- ant.ss %>% st_set_geometry(NULL) %>% 
  group_by(BDM, TypeOfSample) %>% 
  summarise(nTubes=n()) %>% ungroup %>% spread(TypeOfSample, nTubes)
ant.ss.boxTot[is.na(ant.ss.boxTot)] <- 0
site.box.sf <- site.box.sf %>%
  full_join(., ant.ss.boxTot, by="BDM")
ant.cs <- readxl::read_xlsx("data/200114_OperationFourmisData.xlsx", 1) %>%
  rename(TubeNo=CATALOGUENUMBER, SPECIESID=SPECISID)
ant.cs.swisscoord <- ant.cs %>% filter(is.na(LATITUDE)) %>%
  filter(!is.na(SWISSCOORDINATE_X)) %>%
  st_as_sf(coords=c("SWISSCOORDINATE_X", "SWISSCOORDINATE_Y")) %>%
  st_set_crs(21781) %>%
  select(TubeNo, SPECIESID)
ant.cs.latlon <- ant.cs %>% filter(!is.na(LONGITUDE)) %>%
  st_as_sf(coords=c("LONGITUDE", "LATITUDE")) %>% st_set_crs(4326) %>% 
  st_transform(21781) %>%
  select(TubeNo, SPECIESID)
ant.cs <- rbind(ant.cs.swisscoord, ant.cs.latlon)
ant.all <- rbind(ant.ss %>% select(TubeNo, SPECIESID) %>% mutate(source="s"),
                 ant.cs %>% mutate(source="c")) %>%
  mutate(el=raster::extract(dem, .))
grid.cs <- pointcount(dem_agg, st_coordinates(ant.cs)) %>%
  raster::as.data.frame(xy=TRUE) %>% 
  st_as_sf(coords=c("x", "y"), remove=F) %>% st_set_crs(21781) %>% 
  rename(nTubes=CH_21781) %>% filter(st_intersects(., VD, sparse=F))






#--- calculate richness
ant.all <- ant$all %>% 
  mutate(el=raster::extract(dem, .))
interp.rng <- ant.all %>% filter(!is.na(SPECIESID)) %>% 
  st_set_geometry(NULL) %>%
  group_by(SPECIESID) %>%
  summarise(minEl=min(el, na.rm=T),
            maxEl=max(el, na.rm=T),
            medEl=median(el, na.rm=T)) %>%
  arrange(minEl, maxEl) %>%
  mutate(SpElOrder=row_number(),
         SPECIESID_el=factor(SPECIESID, levels=SPECIESID[SpElOrder]))
interp.rng.source <- ant.all %>% filter(!is.na(SPECIESID)) %>% 
  st_set_geometry(NULL) %>%
  group_by(SPECIESID, source) %>%
  summarise(minEl=min(el, na.rm=T),
            maxEl=max(el, na.rm=T),
            medEl=median(el, na.rm=T))

ant.S <- data.frame(el_seq=seq(300, 2400, 100),
                    S=NA,
                    S.c=NA,
                    S.s=NA,
                    S.u=NA,
                    S.c.u=NA,
                    S.s.u=NA,
                    area=NA,
                    MAT=NA,
                    MAP=NA,
                    TAR=NA)
ant.N <- data.frame(el_seq=ant.S$el_seq,
                    N=NA,
                    N.c=NA,
                    N.s=NA,
                    area=NA,
                    MAT=NA,
                    MAP=NA,
                    TAR=NA)
for(i in seq_along(ant.S$el_seq)) {
  dem_i <- dem
  dem_i[dem_i < ant.S$el_seq[i]] <- NA
  dem_i[dem_i > (ant.S$el_seq[i]+100)] <- NA
  ant.S$S[i] <- with(interp.rng, 
               sum(minEl < ant.S$el_seq[i]+100 & maxEl >= ant.S$el_seq[i]))
  ant.S$S.c[i] <- with(filter(interp.rng.source, source=="p"),
                 sum(minEl < ant.S$el_seq[i]+100 & maxEl >= ant.S$el_seq[i]))
  ant.S$S.s[i] <- with(filter(interp.rng.source, source=="s"),
                 sum(minEl < ant.S$el_seq[i]+100 & maxEl >= ant.S$el_seq[i]))
  ant.S$S.u[i] <- with(ant.all %>% filter(!is.na(SPECIESID)) %>%
                         filter(el < ant.S$el_seq[i]+100 & el >= ant.S$el_seq[i]),
                       n_distinct(SPECIESID))
  ant.S$S.c.u[i] <- with(ant.all %>% filter(!is.na(SPECIESID)) %>%
                     filter(source=="p") %>%
                     filter(el < ant.S$el_seq[i]+100 & el >= ant.S$el_seq[i]),
                 n_distinct(SPECIESID))
  ant.S$S.s.u[i] <- with(ant.all %>% filter(!is.na(SPECIESID)) %>%
                     filter(source=="s") %>%
                     filter(el < ant.S$el_seq[i]+100 & el >= ant.S$el_seq[i]),
                   n_distinct(SPECIESID))
  ant.S$area[i] <- ant.N$area[i] <- sum(dem@data@values %/% 100 * 100 == ant.S$el_seq[i], na.rm=T)
  # ant.S$MAT[i] <- ant.N$MAT[i] <- raster::cellStats(raster::mask(MAT, dem_i), 'mean')
  # ant.S$MAP[i] <- ant.N$MAP[i] <- raster::cellStats(raster::mask(MAP, dem_i), 'mean')
  # ant.S$TAR[i] <- ant.N$TAR[i] <- raster::cellStats(raster::mask(TAR, dem_i), 'mean')
  ant.N$N[i] <- with(ant.all %>% filter(!is.na(SPECIESID)) %>%
                       filter(el < ant.S$el_seq[i]+100 & el >= ant.S$el_seq[i]),
                     length(el))
  ant.N$N.c[i] <- with(ant.all %>% filter(!is.na(SPECIESID)) %>%
                         filter(source=="p") %>%
                         filter(el < ant.S$el_seq[i]+100 & el >= ant.S$el_seq[i]),
                     length(el))
  ant.N$N.s[i] <- with(ant.all %>% filter(!is.na(SPECIESID)) %>%
                         filter(source=="s") %>%
                         filter(el < ant.S$el_seq[i]+100 & el >= ant.S$el_seq[i]),
                     length(el))
}

names(ant.S) <- c("Elevation", 
                  "S_all_int", "S_public_int", "S_structured_int",
                  "S_all_obs", "S_public_obs", "S_structured_obs", 
                  "area", "MAT", "MAP", "TAR")
ant.S.long <- gather(ant.S, source, S, 2:7) %>% 
  mutate(dataset=factor(str_split_fixed(source, "_", 3)[,2]),
         interp=str_split_fixed(source, "_", 3)[,3],
         dataset=forcats::lvls_reorder(dataset, c(2,3,1)))

names(ant.N) <- c("Elevation", 
                  "N_all", "N_public", "N_structured",
                  "area", "MAT", "MAP", "TAR")
ant.N.long <- gather(ant.N, source, N, 2:4) %>% 
  mutate(dataset=factor(str_split_fixed(source, "_", 2)[,2]),
         dataset=forcats::lvls_reorder(dataset, c(2,3,1)))



# Richness curves
ggplot(ant.S.long, aes(x=Elevation, y=S, colour=dataset)) + fonts +
  geom_point() + geom_line() + facet_grid(.~interp) +
  scale_colour_brewer(type="qual", palette=2) + 
  theme(legend.position=c(0.9, 0.8))
ggsave("eda/S_all_el.pdf", width=9, height=4.5, units="in")

ggplot(ant.S.long, aes(x=log(area), y=log(S), colour=dataset)) + fonts +
  stat_smooth(method="loess", se=F, span=2, size=0.5) + 
  geom_point() + facet_grid(.~interp) +
  scale_colour_brewer(type="qual", palette=2) + 
  theme(legend.position=c(0.9, 0.2))
ggsave("eda/S_all_area.pdf", width=9, height=4.5, units="in")

ggplot(ant.S.long, aes(x=MAT, y=S, colour=dataset)) + fonts +
  stat_smooth(method="loess", se=F, span=2, size=0.5) + 
  geom_point() + facet_grid(.~interp) +
  scale_colour_brewer(type="qual", palette=2) + 
  theme(legend.position=c(0.1, 0.8))
ggsave("eda/S_all_MAT.pdf", width=9, height=4.5, units="in")

ggplot(ant.S.long, aes(x=MAP, y=S, colour=dataset)) + fonts +
  stat_smooth(method="loess", se=F, span=2, size=0.5) + 
  geom_point() + facet_grid(.~interp) +
  scale_colour_brewer(type="qual", palette=2) + 
  theme(legend.position=c(0.9, 0.8))
ggsave("eda/S_all_MAP.pdf", width=9, height=4.5, units="in")

ggplot(ant.S.long, aes(x=TAR, y=S, colour=dataset)) + fonts +
  stat_smooth(method="loess", se=F, span=2, size=0.5) + 
  geom_point() + facet_grid(.~interp) +
  scale_colour_brewer(type="qual", palette=2) + 
  theme(legend.position=c(0.9, 0.8))
ggsave("eda/S_all_TAR.pdf", width=9, height=4.5, units="in")




# Abundance curves
ggplot(ant.N.long, aes(x=Elevation, y=N, colour=dataset)) + fonts +
  geom_point() + geom_line() + 
  scale_colour_brewer(type="qual", palette=2) + 
  theme(legend.position=c(0.85, 0.8))
ggsave("eda/N_all_el.pdf", width=5, height=4.5, units="in")

ggplot(ant.N.long, aes(x=log(area), y=log(N+1), colour=dataset)) + fonts +
  stat_smooth(method="loess", se=F, span=2, size=0.5) + 
  geom_point() + 
  scale_colour_brewer(type="qual", palette=2) + 
  theme(legend.position=c(0.85, 0.2))
ggsave("eda/N_all_area.pdf", width=5, height=4.5, units="in")

ggplot(ant.N.long, aes(x=MAT, y=log(N+1), colour=dataset)) + fonts +
  stat_smooth(method="loess", se=F, span=2, size=0.5) + 
  geom_point() + 
  scale_colour_brewer(type="qual", palette=2) + 
  theme(legend.position=c(0.15, 0.8))
ggsave("eda/N_all_MAT.pdf", width=5, height=4.5, units="in")

ggplot(ant.N.long, aes(x=MAP, y=log(N+1), colour=dataset)) + fonts +
  stat_smooth(method="loess", se=F, span=2, size=0.5) + 
  geom_point() + 
  scale_colour_brewer(type="qual", palette=2) + 
  theme(legend.position=c(0.15, 0.2))
ggsave("eda/N_all_MAP.pdf", width=5, height=4.5, units="in")

ggplot(ant.N.long, aes(x=TAR, y=log(N+1), colour=dataset)) + fonts +
  stat_smooth(method="loess", se=F, span=2, size=0.5) + 
  geom_point() + 
  scale_colour_brewer(type="qual", palette=2) + 
  theme(legend.position=c(0.15, 0.2))
ggsave("eda/N_all_TAR.pdf", width=5, height=4.5, units="in")










# Map: point locations
ggplot() + fonts + labs(x="", y="") +
  geom_sf(data=ant.cs, alpha=0.25) + 
  geom_sf(data=VD, colour="gray70", fill=NA) +
  geom_sf(data=site.box.sf, colour="red", size=1)
ggsave("eda/map_sampDesign.pdf", width=7, height=6.5, units="in")

ggplot() + fonts + labs(x="", y="") +
  geom_sf(data=ant.all, aes(colour=source), alpha=0.25) + 
  geom_sf(data=VD, colour="gray70", fill=NA) + 
  scale_colour_brewer(type="qual", palette=2, labels=c("Public", "Structured")) +
  theme(legend.position=c(0.12, 0.88))
ggsave("eda/map_allPoints.pdf", width=7, height=6.5, units="in")

ggplot() + fonts + labs(x="", y="") + 
  geom_tile(data=grid.cs, aes(x, y, fill=nTubes>0)) + 
  geom_sf(data=VD, colour="gray70", fill=NA) + 
  scale_fill_manual("nTubes", values=c("gray90", "red"), labels=c("0", "1+")) +
  theme(legend.position=c(0.12, 0.88))
ggsave("eda/map_csPres.pdf", width=7, height=6.5, units="in")

ggplot() + fonts + labs(x="", y="") +
  geom_tile(data=grid.cs %>% mutate(nTubes=na_if(nTubes, 0)), aes(x, y, fill=nTubes)) +
  geom_sf(data=VD, colour="gray70", fill=NA) +
  scale_fill_viridis() +
  theme(legend.position=c(0.12, 0.83))
ggsave("eda/map_csNTubes.pdf", width=7, height=6.5, units="in")

ggplot() + fonts + labs(x="", y="") +
  geom_sf(data=VD, colour="gray70", fill=NA) + 
  geom_sf(data=site.box.sf, 
          aes(fill=(soil+transect+tree)/(area/100),
              colour=(soil+transect+tree)/(area/100)), size=1) + 
  scale_fill_viridis("Tubes/km2\n(all)") +
  scale_colour_viridis("Tubes/km2\n(all)") +
  theme(legend.position=c(0.12, 0.83))
ggsave("eda/map_ssNTubes.pdf", width=7, height=6.5, units="in")

ggplot() + fonts + labs(x="", y="") +
  geom_sf(data=VD, colour="gray70", fill=NA) + 
  geom_sf(data=site.box.sf, aes(fill=soil/n_plots, colour=soil/n_plots), size=1) + 
  scale_fill_viridis("nTubes/.75m2\n(soil)") +
  scale_colour_viridis("nTubes/.75m2\n(soil)") +
  theme(legend.position=c(0.12, 0.83))
ggsave("eda/map_ssNTubesSoil.pdf", width=7, height=6.5, units="in")

ggplot() + fonts + labs(x="", y="") +
  geom_sf(data=VD, colour="gray70", fill=NA) + 
  geom_sf(data=site.box.sf, aes(fill=transect/tran_len, 
                                colour=transect/tran_len), size=1) + 
  scale_fill_viridis("nTubes/m\n(transect)") +
  scale_colour_viridis("nTubes/m\n(transect)") +
  theme(legend.position=c(0.12, 0.83))
ggsave("eda/map_ssNTubesTransect.pdf", width=7, height=6.5, units="in")

ggplot() + fonts + labs(x="", y="") +
  geom_tile(data=grid.cs %>% mutate(nTubes=na_if(nTubes, 0)), 
            aes(x, y, fill=nTubes/max(grid.cs$nTubes))) +
  geom_sf(data=site.box.sf, fill=NA,
          aes(colour=(soil+transect+tree)/83), size=1) + 
  geom_sf(data=VD, colour="gray70", fill=NA) +
  scale_fill_viridis("Relative\nDensity", limits=c(0,1)) +
  scale_colour_viridis("Relative\nDensity", limits=c(0,1)) +
  theme(legend.position=c(0.12, 0.83))
ggsave("eda/map_allRelDens1km2.pdf", width=7, height=6.5, units="in")






# plot-level
ant.canopy.prop <- ant.ss %>% st_set_geometry(NULL) %>%
  group_by(BDM, Canopy) %>% 
  summarise(nTubes=n(),
            bio1_tmean_8110=mean(bio1_tmean_8110, na.rm=T),
            bio7_tar_8110=mean(bio7_tar_8110, na.rm=T),
            bio12_p_8110=mean(bio12_p_8110, na.rm=T),
            mnt25=mean(mnt25, na.rm=T)) %>%
  group_by(BDM) %>% mutate(propTubes=nTubes/sum(nTubes)) %>%
  full_join(select(site.lc.prop, -BDM_id, -Sample_Date) %>%
              group_by(BDM, Canopy) %>%
              summarise(pctArea=sum(pctArea), 
                        soil_plots=sum(soil_plots),
                        trans_len=sum(trans_len)), ., 
            by=c("BDM", "Canopy")) %>%
  mutate(propArea=pctArea/100) %>%
  left_join(., site_summary, by="BDM") 

ggplot(ant.canopy.prop, aes(mnt25, propTubes)) + 
  geom_hline(yintercept=0, colour="gray80") + geom_point(aes(colour=region)) + 
  stat_smooth(se=F, method="loess", span=2, colour="gray30", size=0.5) + 
  fonts + facet_wrap(~Canopy) + 
  labs(x="Elevation (m)", y=expression(frac(tubes[i],sum(tubes[site]))))
ggsave("eda/canopy_El_propTube.pdf", width=10, height=4, units="in")

ggplot(ant.canopy.prop, aes(propArea, propTubes)) + xlim(0,1) + ylim(0,1) +
  geom_abline(colour="gray80") + geom_point(aes(colour=region)) + 
  stat_smooth(se=F, method="loess", span=2, colour="gray30", size=0.5) + 
  fonts + facet_wrap(~Canopy) + 
  labs(x=expression(frac(area[i],sum(area[site]))), 
       y=expression(frac(tubes[i],sum(tubes[site]))))
ggsave("eda/canopy_propArea_propTube.pdf", width=10, height=4, units="in")

ggplot(ant.canopy.prop, aes(mnt25, propTubes-propArea)) + 
  geom_hline(yintercept=0, colour="gray80") + geom_point(aes(colour=region)) + 
  stat_smooth(se=F, method="loess", span=2, colour="gray30", size=0.5) + 
  fonts + facet_wrap(~Canopy) + 
  labs(x="Elevation (m)", 
       y=expression(frac(tubes[i],sum(tubes[site]))-frac(area[i],sum(area[site]))))
ggsave("eda/canopy_El_propDiff.pdf", width=10, height=4, units="in")

ggplot(ant.canopy.prop, aes(bio1_tmean_8110/100, propTubes-propArea)) + 
  geom_hline(yintercept=0, colour="gray80") + geom_point(aes(colour=region)) + 
  stat_smooth(se=F, method="loess", span=2, colour="gray30", size=0.5) + 
  fonts + facet_wrap(~Canopy) + 
  labs(x="Mean Annual Temperature (ºC)", 
       y=expression(frac(tubes,sum(site))-frac(area,sum(site))))

ggplot(ant.canopy.prop, aes(bio7_tar_8110/100, propTubes-propArea)) + 
  geom_hline(yintercept=0, colour="gray80") + geom_point(aes(colour=region)) + 
  stat_smooth(se=F, method="loess", span=2, colour="gray30", size=0.5) + 
  fonts + facet_wrap(~Canopy) + 
  labs(x="Annual Temperature Range (ºC)", 
       y=expression(frac(tubes,sum(site))-frac(area,sum(site))))

ggplot(ant.canopy.prop, aes(bio12_p_8110, propTubes-propArea)) + 
  geom_hline(yintercept=0, colour="gray80") + geom_point(aes(colour=region)) + 
  stat_smooth(se=F, method="loess", span=2, colour="gray30", size=0.5) + 
  fonts + facet_wrap(~Canopy) + 
  labs(x="Annual Precipitation (mm)", 
       y=expression(frac(tubes,sum(site))-frac(area,sum(site))))



ant.plot.prop <- ant.ss %>% st_set_geometry(NULL) %>%
  group_by(BDM, Canopy, Plot_id) %>% 
  summarise(nTubes=n(),
            SoilTemp=mean(SoilTemp, na.rm=T),
            SoilTempAnomaly=mean(SoilTempAnomaly, na.rm=T),
            bio1_tmean_8110=mean(bio1_tmean_8110, na.rm=T),
            bio7_tar_8110=mean(bio7_tar_8110, na.rm=T),
            bio12_p_8110=mean(bio12_p_8110, na.rm=T),
            mnt25=mean(mnt25, na.rm=T),
            Grass=mean(Grass),
            Shrub=mean(Shrub),
            Bare=mean(Bare),
            Litter=mean(Litter)) %>%
  group_by(BDM) %>% mutate(propTubes=nTubes/sum(nTubes)) %>%
  full_join(select(site.lc.prop, -BDM_id, -Sample_Date) %>%
              group_by(BDM, Canopy) %>%
              summarise(pctArea=sum(pctArea), 
                        soil_plots=sum(soil_plots),
                        trans_len=sum(trans_len)), ., 
            by=c("BDM", "Canopy")) %>%
  mutate(propArea=pctArea/100) %>%
  left_join(., site_summary, by="BDM") 

ggplot(ant.plot.prop, aes(SoilTempAnomaly, propTubes)) + 
  geom_hline(yintercept=0, colour="gray80") + geom_point(aes(colour=region)) + 
  stat_smooth(se=F, method="loess", span=2, colour="gray30", size=0.5) + 
  fonts + facet_wrap(~Canopy) + 
  labs(x="Soil Temp Anomaly (ºC)", y=expression(frac(tubes[i],sum(tubes[site]))))

ggplot(ant.plot.prop, aes(Litter, propTubes)) + 
  geom_hline(yintercept=0, colour="gray80") + geom_point(aes(colour=region)) + 
  stat_smooth(se=F, method="loess", span=2, colour="gray30", size=0.5) + 
  fonts + #facet_wrap(~Canopy) + 
  labs(x="Litter", y=expression(frac(tubes[i],sum(tubes[site]))))

ggplot(ant.plot.prop, aes(propArea, propTubes)) + 
  geom_hline(yintercept=0, colour="gray80") + geom_point(aes(colour=region)) + 
  stat_smooth(se=F, method="loess", span=2, colour="gray30", size=0.5) + 
  fonts + facet_wrap(~Canopy) + 
  labs(x="Elevation (m)", y=expression(frac(tubes[i],sum(tubes[site]))))
