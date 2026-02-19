# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.8.22

# This code is to make maps for the two boundaries of projected frequency
# By rotating alpha around a fixed point
# Results are mean of 15 models

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

# Input path of projected frequency rasters
Freq_rasters <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_kb/Model_mean_End_ssp585_freq_kb.rds")
# Make CONUS boundary
# Whole US map
CONUS <- st_read("/fs/ess/PAS2204/Data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
# CONUS outer boundary map
#CONUS <- st_union(CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",])
CONUS <- CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",]

# This is the color for frequency
color_freq <- RColorBrewer::brewer.pal(9,"YlOrRd")
############
# Functions
############

my_theme <- theme(
  #axis.line=element_line(color="black"),
  panel.background = element_blank(),
  #text = element_text(size=100),
  #panel.border = element_rect(colour="black",fill=NA),
  legend.key = element_blank(),
  #legend.key.size = unit(6,"cm"),
  #aspect.ratio = 1/1,
  #legend.key.size = unit(0.3,'cm'),
  legend.text = element_text(size=16),
  plot.title = element_text(size=16),
  legend.title = element_text(size=16),
  #plot.margin = margin(80,80,80,80),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank()
  #legend.margin = margin(0,80,0,80)
  #axis.ticks.length = unit(-0.8,"cm")
)

# This function is to make frequency maps
# Input raster is the raster to plot
# my_color is the color to use 
# my_legend is the legend title
# my_title is the title of the plot
# min and max are the min and max for color scales
g_map <- function(raster,my_color,my_legend,my_title,min,max){
  raster_df <- as.data.frame(raster,xy=TRUE)
  g <- ggplot()+
    geom_tile(data=raster_df,aes(x=x,y=y,fill=raster_df[,3]))+
    geom_sf(data=CONUS,fill=NA,color="black",alpha=0.5)+
    my_theme+
    scale_fill_gradientn(na.value="white",
                         colours = my_color,
                         limits = c(min,max),
                         guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))+
    labs(x="",y="",fill=my_legend,title=my_title)+
    theme(legend.position = "right")  
  return(g)
}

########
# Main
########
g1 <- g_map(Freq_rasters[[1]],color_freq,
            my_legend = "Frequency",
            my_title = "Boundary 1",
            min=0,max=1)
g2 <- g_map(Freq_rasters[[2]],color_freq,
            my_legend = "Frequency",
            my_title = "Boundary 2",
            min=0,max=1)



