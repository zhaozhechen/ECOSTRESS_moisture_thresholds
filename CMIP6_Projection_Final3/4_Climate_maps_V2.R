# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2025.6.23

# This code is to make projected change in VPD and SM

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

# Input path for historical observation
#hist_obs_daily_CONUS_mean <- read.csv("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Hist_obs/Hist_obs_daily_CONUS_mean.csv") 
hist_obs_gs_mean <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Hist_obs/Hist_obs_gs_mean.rds")
# Input path for adjusted CMIP projection
CMIP_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/CMIP6_adjusted/"
# Output path for plots
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Figures_V3/"
# Make CONUS boundary
# Whole US map
CONUS <- st_read("/fs/ess/PAS2204/Data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
# CONUS outer boundary map
#CONUS <- st_union(CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",])
CONUS <- CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",]

# Color palettes for VPD and SM
color_VPD <- rev(heat.colors(100))
color_SM <- RColorBrewer::brewer.pal(11, "BrBG")
color_SM <- c(color_SM[1:5],color_SM[9:10])

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

# This function is to read CMIP raster
# time_name is "Mid" or "End"
# ssp is "ssp245" or "ssp585"
# varname is "daily_SM" or "daily_maxVPD"
CMIP_raster <- function(time_name,ssp,varname){
  raster <- readRDS(paste0(CMIP_path,"Mean_models_gs_mean_",time_name,"_",ssp,"_",varname,".rds"))
  return(raster)
}

# This function is to make maps
# Input raster is the raster to plot
# my_color is either color_VPD or color_SM
# min and max are the min and max values for color bars
# my_title is the title for the plot
# my_legend is the title for the legend
p_maps <- function(raster,my_color,min,max,my_title,my_legend){
  raster_df <- as.data.frame(raster,xy=TRUE)
  g <- ggplot()+
    geom_tile(data=raster_df,aes(x=x,y=y,fill=raster_df[,3]))+
    geom_sf(data=CONUS,fill=NA,color="black",alpha=0.5)+
    my_theme+
    scale_fill_gradientn(na.value="white",
                         colours = my_color,
                         limits = c(min,max),
                         # values out of bounds are limited to the bounds
                         oob = scales::squish,
                         guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))+
    labs(x="",y="",title = my_title,fill = my_legend)
  return(g)
}

##########
# Main
##########
# Make maps of delta_VPD
# Get delta_VPD raster
delta_VPD <- CMIP_raster("End","ssp585","daily_maxVPD") - hist_obs_gs_mean$daily_maxVPD
# Make a map for it
g_delta_VPD <- p_maps(delta_VPD,color_VPD,0,1.6,"",bquote(Delta~VPD~(kPa)))
# Get delta_SM raster
delta_SM <- CMIP_raster("End","ssp585","daily_SM") - hist_obs_gs_mean$daily_SM
# Make a map for it
g_delta_SM <- p_maps(delta_SM,color_SM,-0.03,0.01,"",bquote(Delta~SM~"(m"^3~"/m"^3~")"))

# Combine these 2 maps
g <- plot_grid(g_delta_VPD,g_delta_SM,
               nrow = 1,
               labels = "auto")
pdf(paste0(Output_path,"delta_Climate_maps_End_ssp585_SM_VPD.pdf"),
    height=4,
    width=15)
print(g)
dev.off()






