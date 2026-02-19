# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.8.7
# Update: 2025.6.23

# This code is to plot delta_freqency
# between frequency from historical climate
# vs frequency from 15 models (changing only SM or VPD, 15 models)
# For 3 scenarios:
# (1) Only Change SM to CMIP6 SM in 2100 (ssp585)
# (2) Only Change VPD to CMIP6 VPD in 2100 (ssp585)
# While keeping the other variable at the original historical values
# (3) Change both VPD and SM to CMIP6 projection in 2100 (ssp585)
# (4) 3 - (1+2)

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

# Input path of projected frequency changing only SM or VPD
Input_path1 <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_SM_VPD/"
# Input path of historical frequency and changing both SM and VPD
Input_path2 <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_V2/"
# Input path of combined 0.25D rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")
# Make CONUS boundary
# Whole US map
CONUS <- st_read("/fs/ess/PAS2204/Data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
# CONUS outer boundary map
#CONUS <- st_union(CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",])
CONUS <- CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",]
# Output path of figures
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Figures_V3/"

# This is the color for frequency
color_freq <- RColorBrewer::brewer.pal(9,"YlOrRd")
# This is the color for delta frequency
color_delta_freq <- rev(RColorBrewer::brewer.pal(11,"PiYG"))[c(2,6,7,8,9,10)]

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

# This function is to make delta_frequency maps
# Input raster is the raster to plot
# p_value is the corresponding p-value to mask the map
# my_color is the color to use
# my_title is the title of the plot
g_delta_map <- function(raster,p_value,my_color,my_title){
  df <- as.data.frame(raster,xy=TRUE)
  df$p_value <- values(p_value)
  # Only show significant pixels
  df$layer[df$p_value > 0.05] <- NA
  df$layer[df$layer >= 0.4] <- 0.4
  df$layer[df$layer <= -0.1] <- -0.1
  g <- ggplot()+
    geom_tile(data=df,aes(x=x,y=y,fill=layer))+
    geom_sf(data=CONUS,fill=NA,color="black",alpha=0.5)+
    my_theme+
    scale_fill_gradientn(na.value="white",
                         colours = my_color,
                         guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"),
                         breaks = c(-0.1,0,0.1,0.2,0.3,0.4),
                         labels = c(expression(""<= -0.1),0,0.1,0.2,0.3,expression("">=0.4)))+
    labs(x="",y="",fill=bquote(Delta~Frequency),
         title = paste0(time_name,"-Century ",ssp," ",my_title))+
    theme(legend.position = "right") 
  
  return(g)
}

########
# Main
########
# Read in data and pre-processing ============================================
time_name <- "End"
ssp <- "ssp585"

# Historical frequency
Hist_freq <- readRDS(paste0(Input_path2,"Hist_obs_freq.rds"))$Freq
Hist_freq[Hist_freq == 1] <- NA
# Mask by theta
Hist_freq <- mask(Hist_freq,raster_all$Theta)

# Model mean frequency, changing both VPD and SM
Proj_freq <- readRDS(paste0(Input_path2,"Model_mean_freq_",time_name,"_",ssp,".rds"))$Freq_mean
# Mask by theta
Proj_freq_SM_VPD <- mask(Proj_freq,raster_all$Theta)
# Calculate delta_frequency,using both projected SM and VPD
delta_freq_SM_VPD <- Proj_freq_SM_VPD - Hist_freq
# The corresponding p-value raster
p_value_SM_VPD <- readRDS(paste0(Input_path2,"p_value_rasters.rds"))
p_value_SM_VPD <- p_value_SM_VPD[[which(names(p_value_SM_VPD) == paste0("p_value_",time_name,"_",ssp))]]

# Model mean frequency, changing only SM
Proj_freq <- readRDS(paste0(Input_path1,"Model_mean_",time_name,"_",ssp,"_SM_chg_freq.rds"))$Freq_mean
# Mask by theta
Proj_freq_SM <- mask(Proj_freq,raster_all$Theta)
# Calculate delta_frequency,using only projected SM
delta_freq_SM <- Proj_freq_SM - Hist_freq
# The corresponding p-value raster
p_value_SM <- readRDS(paste0(Input_path1,"p_value_rasters.rds"))$p_value_End_ssp585_SM_chg

# Model mean frequency, changing only VPD
Proj_freq <- readRDS(paste0(Input_path1,"Model_mean_",time_name,"_",ssp,"_VPD_chg_freq.rds"))$Freq_mean
# Mask by theta
Proj_freq_VPD <- mask(Proj_freq,raster_all$Theta)
# Calculate delta_frequency,using only projected VPD
delta_freq_VPD <- Proj_freq_VPD - Hist_freq
# The corresponding p-value raster
p_value_VPD <- readRDS(paste0(Input_path1,"p_value_rasters.rds"))$p_value_End_ssp585_VPD_chg

# Make maps of frequencies =====================================================
# Historical frequency
Hist_map <- g_map(Hist_freq,color_freq,
                  my_legend = "Frequency",my_title = "Historical",
                  min = 0,max=1)
# Projected frequency,changing both SM and VPD
Proj_SM_VPD_map <- g_map(Proj_freq_SM_VPD,color_freq,
                         my_legend = "Frequency",my_title = paste(time_name,ssp,"Projected SM and VPD"),
                         min = 0,max=1)
# Projected frequency,changing only SM
Proj_SM_map <- g_map(Proj_freq_SM,color_freq,
                     my_legend = "Frequency",my_title = paste(time_name,ssp,"Projected SM and Historical VPD"),
                     min = 0,max=1)
# Projected frequency,changing only VPD
Proj_VPD_map <- g_map(Proj_freq_VPD,color_freq,
                     my_legend = "Frequency",my_title = paste(time_name,ssp,"Projected VPD and Historical SM"),
                     min = 0,max=1)
# Combine these 4 maps
g <- plot_grid(Hist_map,Proj_SM_map,Proj_VPD_map,Proj_SM_VPD_map,
               nrow=2)
pdf(paste0(Output_path,"Frequency_maps_",time_name,"_",ssp,"_SM_VPD.pdf"),
    height=8,
    width=15)
print(g)
dev.off()

# Make maps of delta_frequency =================================================
delta_freq_SM_VPD_map <- g_delta_map(delta_freq_SM_VPD,p_value_SM_VPD,color_delta_freq,"Projected SM and VPD")
delta_freq_SM_map <- g_delta_map(delta_freq_SM,p_value_SM,color_delta_freq,"Projected SM and Historical VPD")
delta_freq_VPD_map <- g_delta_map(delta_freq_VPD,p_value_VPD,color_delta_freq,"Projected VPD and Historical SM")
# Make a fourth one, which is delta_freq_SM_VPD_map - (delta_freq_SM_map + delta_freq_VPD_map)

# Combine these 3 maps
g <- plot_grid(delta_freq_SM_map,delta_freq_VPD_map,delta_freq_SM_VPD_map,
               nrow=2)
#pdf(paste0(Output_path,"delta_Frequency_maps_",time_name,"_",ssp,"_SM_VPD.pdf"),
#    height=8,
#    width=15)
#print(g)
#dev.off()




