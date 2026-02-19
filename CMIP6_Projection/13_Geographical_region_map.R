# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2025.6.23

# This code is to make maps of geographic regions and divisions

#########
# Global
#########

dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
library(raster)
library(ncdf4)
library(sf)
library(ggplot2)
library(cowplot)
library(wesanderson)
library(dplyr)
library(tidyr)
library(ggpattern)
# This is for the strokes of texts
library(shadowtext)
library(ggrepel)

# CONUS
CONUS <- st_read("/fs/ess/PAS2204/Data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
CONUS <- CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",]

# Region
CONUS_region <- st_read("/fs/ess/PAS2204/Data/cb_2018_us_region_20m/cb_2018_us_region_20m.shp") %>%
  select(NAME)
# Division
CONUS_division <- st_read("/fs/ess/PAS2204/Data/cb_2018_us_division_20m/cb_2018_us_division_20m.shp") %>%
  select(NAME)
# Crop Region and Division
CONUS_region_shp <- st_intersection(CONUS_region,st_union(CONUS))
CONUS_division_shp <- st_intersection(CONUS_division,st_union(CONUS))

# Input path of combined 0.25D rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")
# names of the 13 sites
Contour_features_df <- read.csv("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/Validation_Final/Contour_feature_df.csv")
# Their coordinates
Site_info_ls <- read.csv("/fs/ess/PAS2204/SharedData/AmeriFlux_All_Sites/AMF_All_Sites_info.csv")
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Figures_V3/"

############
# Functions
############
my_theme <- theme(
  panel.background = element_blank(),
  legend.key = element_blank(),
  legend.text = element_text(size=14),
  plot.title = element_text(size=14),
  #plot.margin = margin(80,80,80,80),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size=14)
  #legend.margin = margin(0,80,0,80)
  #axis.ticks.length = unit(-0.8,"cm")
)

########
# Main
########

# Rasterize CONUS shapefile to match raster_all
# For Region
CONUS_region <- CONUS_region_shp %>%
  filter(st_geometry_type(.) %in% c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION"))
CONUS_region <- st_collection_extract(CONUS_region, "POLYGON")
CONUS_region <- st_cast(CONUS_region, "MULTIPOLYGON")
CONUS_region$region_id <- as.numeric(factor(CONUS_region$NAME))
region_raster <- rasterize(CONUS_region,raster_all$Theta,field = "region_id")
raster_all$region <- region_raster
rm(region_raster)

# For Division
CONUS_division <- CONUS_division_shp %>%
  filter(st_geometry_type(.) %in% c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION"))
CONUS_division <- st_collection_extract(CONUS_division, "POLYGON")
CONUS_division <- st_cast(CONUS_division, "MULTIPOLYGON")
CONUS_division$division_id <- as.numeric(factor(CONUS_division$NAME))
division_raster <- rasterize(CONUS_division,raster_all$Theta,field = "division_id")
raster_all$division <- division_raster
rm(division_raster)

# Get region names for each region_id
region_names <- CONUS_region$NAME
names(region_names) <- CONUS_region$region_id

# Get division names for each division_id
#division_names <- CONUS_division$NAME
#names(division_names) <- CONUS_division$division_id

# Find centroids for each division
CONUS_division_shp <- st_make_valid(CONUS_division_shp)
division_centro <- CONUS_division_shp %>%
  st_centroid(of_largest_polygon = TRUE) %>%
  mutate(label = gsub(" ","\n",NAME),
         x = st_coordinates(geometry)[,1],
         y = st_coordinates(geometry)[,2])

# Make map of geographical regions
raster_df <- as.data.frame(raster_all$region,xy=TRUE)

g_GEO <- ggplot()+
  geom_tile(data=raster_df,aes(x=x,y=y,fill=factor(region)))+
  scale_fill_brewer(palette = "Set2",na.value=NA,labels = region_names)+
  geom_sf(data=CONUS,fill=NA,color="grey",alpha=0.7,size = 0.1)+
  geom_sf(data=CONUS_division_shp,fill=NA,color="black",alpha=0.7,size=5)+
  geom_shadowtext(data=division_centro,
                  aes(x = x, y=y, label=label),
                  size = 3.5,
                  color="black",
                  fontface = "bold",
                  bg.color = "white")+
  my_theme+
  labs(x="",y="",fill="Region")

# Make map for PFT
# Get the coordinates of the 13 sites
Site_ls <- data.frame(Site_ID = Contour_features_df$Site_ID)
Site_info_ls <- Site_info_ls[c("Site.ID","Lat","Long")]
names(Site_info_ls)[1] <- "Site_ID"
Site_ls <- merge(Site_ls,Site_info_ls,by="Site_ID")

# Get the PFT df
#raster_all <- mask(raster_all,mask=raster_all$Theta)
raster_df <- as.data.frame(raster_all$LC,xy=TRUE)
names(raster_df)[3] <- "PFT"
raster_df <- raster_df %>%
  mutate(PFT = replace(PFT,PFT==41,"DF"),
         PFT = replace(PFT,PFT==42,"EF"),
         PFT = replace(PFT,PFT==43,"MF"),
         PFT = replace(PFT,PFT==51,"Dwarf scrub"),
         PFT = replace(PFT,PFT==52,"SHB"),
         PFT = replace(PFT,PFT==71,"GRA"),
         PFT = replace(PFT,PFT==72,"Sedge/Herbaceous"),
         PFT = replace(PFT,PFT==81,"PAS"),
         PFT = replace(PFT,PFT==82,"CRO"))
raster_df$PFT <- factor(raster_df$PFT,levels = c("CRO","SHB","GRA","PAS","EF","DF","MF"))

# Make a map of these 13 sites

g_PFT <- ggplot()+
  geom_tile(data=raster_df,aes(x=x,y=y,fill=PFT))+
  scale_fill_brewer(palette = "Set2",na.value=NA)+
  geom_sf(data=CONUS,fill=NA,color="black",alpha=0.7)+
  geom_point(data=Site_ls,aes(x=Long,y=Lat),shape = 21, color="black",fill="yellow",size=3)+
  geom_label_repel(data=Site_ls,aes(x=Long,y=Lat,label=Site_ID),
                   point.padding = 0.5,
                   box.padding = 0.5,
                   segment.color="grey")+
  my_theme+
  labs(x="",y="")

# Put the two maps together
g <- plot_grid(g_PFT,g_GEO,nrow = 1,
               align = "hv",
               axis = "lbrt",
               labels = "auto")
pdf(paste0(Output_path,"PFT+GEO_maps.pdf"),
    height=5,width=15)
print(g)
dev.off()

saveRDS(raster_all,paste0(Output_path,"raster_all_with_GEO.rds"))



