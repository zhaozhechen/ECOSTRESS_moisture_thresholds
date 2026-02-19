# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.6.20

# This code is to make pdf plots and maps of frequency and delta_VPD

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

# Input path of frequency and delta_VPD rasters
Input_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks/"
# Input path of combined 0.25D rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")
# Make CONUS boundary
# Whole US map
CONUS <- st_read("/fs/ess/PAS2204/Data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
# CONUS outer boundary map
#CONUS <- st_union(CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",])
CONUS <- CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",]
# Output path
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Figures/"

# This is the color for frequency
color_freq <- RColorBrewer::brewer.pal(9,"YlOrRd")
# This is the color for delta frequency
color_delta_freq <- rev(RColorBrewer::brewer.pal(11,"PiYG"))[2:10]
# This is the color for delta_VPD
color_delta_VPD <- wes_palette("Zissou1",100,type="continuous")

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
                   legend.title = element_text(size=16))

# This function is to make maps
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

#######
# Main
#######
for(time_name in c("Mid","End")){
  for(ssp in c("ssp245","ssp585")){
    
    # Frequency maps ==============================================================================
    # Historical frequency and delta_VPD
    Hist_freq <- readRDS(paste0(Input_path,"Hist_obs_freq.rds"))
    # Mask by theta
    Hist_freq <- mask(Hist_freq,raster_all$Theta)
    # Map of historical frequency
    g1 <- g_map(Hist_freq$Freq,color_freq,"Frequency","Historical",0,1)
    
    # Model mean frequency and delta_VPD
    Proj_freq <- readRDS(paste0(Input_path,"Model_mean_freq_",time_name,"_",ssp,".rds"))
    # Mask by theta
    Proj_freq <- mask(Proj_freq,raster_all$Theta)
    # Map of projected frequency mean
    g2 <- g_map(Proj_freq$Freq_mean,color_freq,"Frequency",paste0(time_name,"-Century ",ssp),0,1)
    # Map of projected frequency sd
    g3 <- g_map(Proj_freq$Freq_sd,color_freq,"sd(Frequency)",paste0(time_name,"-Century ",ssp),NA,NA)
    
    # Calculate frequency cv as a mask
    freq_mask <- Proj_freq$Freq_sd/Proj_freq$Freq_mean
    # pdf of cv
    df <- as.data.frame(freq_mask)
    g6 <- ggplot(df,aes(x=layer,after_stat(density)))+
      stat_density(position = "identity",geom="line")+
      pdf_theme+
      labs(x="cv(frequency)")
    
    # cv > 1 are masked
    freq_mask[freq_mask>1] <- NA
    # Map of projected frequency mean, masked by cv > 1
    g4 <- g_map(mask(Proj_freq$Freq_mean,freq_mask),
                color_freq,"Frequency",
                paste0(time_name,"-Century ",ssp," (masked by cv>1)"),0,1)
    
    # Calculate delta_Frequency
    delta_freq <- Proj_freq$Freq_mean - Hist_freq$Freq
    # Mask by frequency cv > 1
    delta_freq <- mask(delta_freq,freq_mask)
    #delta_freq[delta_freq>0.4] <- 0.4
    #delta_freq[delta_freq< (-0.4)] <- (-0.4)
    g5 <- g_map(delta_freq,color_delta_freq,bquote(Delta~Frequency),
                paste0(time_name,"-Century ",ssp," (masked by cv>1)"),NA,NA)
    
    # Combine plots for frequency
    g <- plot_grid(g2,g3,
                   nrow=1)
    pdf(paste0(Output_path,"Frequency_maps_",time_name,"_",ssp,"_V2.pdf"),
        height=4,
        width=16)
    print(g)
    dev.off()
    
    # Delta_VPD maps ==========================================================================
    # Historical delta_VPD
    Hist_delta_VPD <- Hist_freq$delta_VPD
    Hist_delta_VPD[Hist_delta_VPD < 0] <- 0
    Hist_delta_VPD[Hist_delta_VPD > 3] <- 3
    # Map of historical delta_VPD
    g1 <- g_map(Hist_delta_VPD,color_delta_VPD,bquote(Delta~Thresholds[VPD]),"Historical",NA,NA)
    
    # Projected delta_VPD
    Proj_delta_VPD <- Proj_freq$Delta_VPD_mean
    Proj_delta_VPD[Proj_delta_VPD < 0] <- 0
    Proj_delta_VPD[Proj_delta_VPD > 3] <- 3
    # Map of projected delta_VPD
    g2 <- g_map(Proj_delta_VPD,color_delta_VPD,bquote(Delta~Thresholds[VPD]),
                paste0(time_name,"-Century ",ssp),NA,NA)
    
    # Mask projected delta_VPD by cv(frequency)
    Proj_delta_VPD <- mask(Proj_delta_VPD,freq_mask)
    # Map of projected delta_VPD (masked by cv)
    g3 <- g_map(Proj_delta_VPD,color_delta_VPD,bquote(Delta~Thresholds[VPD]),
                paste0(time_name,"-Century ",ssp, " (masked by cv>1)"),NA,NA)
    
    # Change in delta_VPD
    delta_delta_VPD <- Proj_delta_VPD - Hist_delta_VPD
    g4 <- g_map(delta_delta_VPD,color_delta_VPD,bquote(Delta~Thresholds[VPD]),
                paste0(time_name,"-Century ",ssp, " Change from historical (masked by cv>1)"),NA,NA)
    # Combine plots for frequency
    g <- plot_grid(g1,g2,g3,g4,
                   nrow=2)
    pdf(paste0(Output_path,"Delta_VPD_maps_",time_name,"_",ssp,".pdf"),
        height=8,
        width=15)
    print(g)
    dev.off()
  }
}

