# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Update Date: 2024.8.12

# This code is to make maps of climates and pdf of daily CONUS mean

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
hist_obs_daily_CONUS_mean <- read.csv("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Hist_obs/Hist_obs_daily_CONUS_mean.csv") 
hist_obs_gs_mean <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Hist_obs/Hist_obs_gs_mean.rds")
# Input path for adjusted CMIP projection
CMIP_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/CMIP6_adjusted/"
# Output path for plots
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Figures_V2/"
# Make CONUS boundary
# Whole US map
CONUS <- st_read("/fs/ess/PAS2204/Data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
# CONUS outer boundary map
#CONUS <- st_union(CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",])
CONUS <- CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",]

# Color palettes for VPD and SM
color_VPD <- wes_palette("Zissou1",100,type="continuous")
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

pdf_theme <- theme(panel.background = element_blank(),
                   panel.border = element_rect(colour="black",fill=NA),
                   legend.key = element_blank(),
                   #legend.key.size = unit(6,"cm"),
                   legend.text = element_text(size=16),
                   plot.title = element_text(size=16),
                   legend.title = element_text(size=16),
                   aspect.ratio = 1/1)

# This function is to read CMIP raster
# time_name is "Mid" or "End"
# ssp is "ssp245" or "ssp585"
# varname is "daily_SM" or "daily_maxVPD"
CMIP_raster <- function(time_name,ssp,varname){
  raster <- readRDS(paste0(CMIP_path,"Mean_models_gs_mean_",time_name,"_",ssp,"_",varname,".rds"))
  return(raster)
}

# This function is to read CMIP daily CONUS mean
# time_name is "Mid" or "End"
# ssp is "ssp245" or "ssp585"
# varname is "daily_SM" or "daily_maxVPD"
CMIP_var <- function(time_name,ssp,varname){
  var <- read.csv(paste0(CMIP_path,"Mean_models_daily_CONUS_mean_",time_name,"_",ssp,"_",varname,".csv"))
  return(var)
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
                         guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))+
    labs(x="",y="",title = my_title,fill = my_legend)
  return(g)
}

# This function is to make all 6 climate maps for the current scenario
# Input time is "Mid" or "End"
# ssp is "ssp245" or ""ssp585
All_maps <- function(time_name,ssp){
  # Stack all rasters to be plot
  raster_stack <- stack(hist_obs_gs_mean$daily_maxVPD,
                        CMIP_raster(time_name,ssp,"daily_maxVPD"),
                        hist_obs_gs_mean$daily_SM,
                        CMIP_raster(time_name,ssp,"daily_SM"),
                        CMIP_raster(time_name,ssp,"daily_maxVPD") - hist_obs_gs_mean$daily_maxVPD,
                        CMIP_raster(time_name,ssp,"daily_SM") - hist_obs_gs_mean$daily_SM)
  # Titles of the maps
  titles <- c("Historical",
              paste0(time_name,"-Century ",ssp),
              "Historical",
              paste0(time_name,"-Century ",ssp),
              "Change in daily maximum VPD",
              "Change in daily mean SM")
  # Use the same color scale for the first two maps, and the second two maps
  VPD_max <- max(max(values(raster_stack[[1]]),na.rm=T),
                 max(values(raster_stack[[2]]),na.rm=T))
  VPD_min <- min(min(values(raster_stack[[1]]),na.rm=T),
                 min(values(raster_stack[[2]]),na.rm=T))
  SM_max <- max(max(values(raster_stack[[3]]),na.rm=T),
                max(values(raster_stack[[4]]),na.rm=T))
  SM_min <- min(min(values(raster_stack[[3]]),na.rm=T),
                min(values(raster_stack[[4]]),na.rm=T))
  # Make maps
  #g1 <- p_maps(raster_stack[[1]],color_VPD,VPD_min,VPD_max,titles[1],"VPD (kPa)")
  #g2 <- p_maps(raster_stack[[2]],color_VPD,VPD_min,VPD_max,titles[2],"VPD (kPa)")
  #g3 <- p_maps(raster_stack[[3]],color_SM,SM_min,SM_max,titles[3],bquote("SM (m"^3~"/m"^3~")"))
  #g4 <- p_maps(raster_stack[[4]],color_SM,SM_min,SM_max,titles[4],bquote("SM (m"^3~"/m"^3~")"))
  g5 <- p_maps(raster_stack[[5]],color_VPD,NA,NA,titles[5],bquote(Delta~VPD~(kPa)))
  g6 <- p_maps(raster_stack[[6]],color_SM,-0.03,0.01,titles[6],bquote(Delta~SM~"(m"^3~"/m"^3~")"))
  # Combine these 6 maps
  g <- plot_grid(g5,g6,
                 nrow=1,
                 labels = "auto",
                 align = "hv",
                 axis="brlt")
  # Output this map
  pdf(paste0(Output_path,"Climate_maps_",time_name,"_",ssp,"_V2.pdf"),
      height=4,
      width=15)
  print(g)
  dev.off()
}

# This function outputs pdf for VPD and SM
p_pdf <- function(time_name){
  # Initialize a list for output
  out <- list()
  # Make df of SM and VPD
  df <- data.frame(
    Type = rep(c("Historical",
                 paste0(time_name," (SSP245)"),
                 paste0(time_name," (SSP585)")),
               each = length(hist_obs_daily_CONUS_mean$daily_SM)),
    SM = c(hist_obs_daily_CONUS_mean$daily_SM,
           CMIP_var(time_name,"ssp245","daily_SM")$x,
           CMIP_var(time_name,"ssp585","daily_SM")$x),
    VPD = c(hist_obs_daily_CONUS_mean$daily_maxVPD,
            CMIP_var(time_name,"ssp245","daily_maxVPD")$x,
            CMIP_var(time_name,"ssp585","daily_maxVPD")$x))
  # Get mean of VPD and SM at different scenario
  df_mean <- df %>%
    group_by(Type) %>%
    summarise_at(.vars = c("VPD","SM"),
                 .funs = "mean",
                 na.rm=T)
  
  g_VPD <- ggplot(df,aes(x=VPD,after_stat(density),color=Type))+
    stat_density(position="identity",geom="line")+
    geom_vline(data=df_mean,aes(xintercept=VPD,color=Type),
               linetype = "dashed",
               show.legend = FALSE)+
    pdf_theme+
    theme(legend.position = "none")+
    labs(x="Daily maximum VPD (kPa)",color="")
  
  g_SM <- ggplot(df,aes(x=SM,after_stat(density),color=Type))+
    stat_density(position="identity",geom="line")+
    geom_vline(data=df_mean,aes(xintercept=SM,color=Type),
               linetype = "dashed",
               show.legend = FALSE)+
    pdf_theme+
    labs(x=bquote("Daily mean SM (m"^3~"/m"^3~")"),color="")
  out[[1]] <- g_VPD
  out[[2]] <- g_SM
  return(out)
}
#######
# Main
#######
# Make maps of climates ==========================================================
# For each scenario (One period, One ssp)
#time_name <- "Mid"
#ssp <- "ssp245"
#All_maps("Mid","ssp245")
#All_maps("Mid","ssp585")
#All_maps("End","ssp245")
All_maps("End","ssp585")
 
# Make pdf of daily max VPD and daily SM =========================================
pdf_Mid <- p_pdf("Mid")
pdf_End <- p_pdf("End")

# Combine plots
g <- plot_grid(pdf_Mid[[1]],pdf_Mid[[2]],pdf_End[[1]],pdf_End[[2]],
               align="hv",
               axis="brlt")
pdf(paste0(Output_path,'Climate_pdf.pdf'),
    height=4,
    width=8)
print(g)
dev.off()

