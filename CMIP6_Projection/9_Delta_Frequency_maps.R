# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.6.27

# This code is to make maps of delta_frequency, masked by p-value

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

# This is the color for delta frequency
color_delta_freq <- rev(RColorBrewer::brewer.pal(11,"PiYG"))

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

########
# Main
########
time_name <- "End"
ssp <- "ssp585"

# Historical frequency
Hist_freq <- readRDS(paste0(Input_path,"Hist_obs_freq.rds"))$Freq
# Mask by theta
Hist_freq <- mask(Hist_freq,raster_all$Theta)

# Model mean frequency
Proj_freq <- readRDS(paste0(Input_path,"Model_mean_freq_",time_name,"_",ssp,".rds"))$Freq_mean
# Mask by theta
Proj_freq <- mask(Proj_freq,raster_all$Theta)

# Calculate delta_Frequency
delta_freq <- Proj_freq - Hist_freq

# The corresponding p-value raster
p_value <- readRDS(paste0(Input_path,"p_value_rasters.rds"))
p_value <- p_value[[which(names(p_value) == paste0("p_value_",time_name,"_",ssp))]]

# Make map of delta_frequency, masked by p-value > 0.05
df <- as.data.frame(delta_freq,xy=TRUE)
df$p_value <- values(p_value)
df$p_value[df$p_value <= 0.05] <- NA

df$layer[df$layer >= 0.4] <- 0.4
df$layer[df$layer <= (-0.2)] <- (-0.2)

g <- ggplot()+
  geom_tile(data=df,aes(x=x,y=y,fill=layer))+
  geom_sf(data=CONUS,fill=NA,color="black",alpha=0.5)+
  geom_point(data=df[!is.na(df$p_value),],
             aes(x=x,y=y),
             shape = 16,
             fill = "black",
             size = 0.5,
             alpha = 0.3)+
  my_theme+
  scale_fill_gradient2(na.value="white",
                       midpoint = 0,
                       #colours = color_delta_freq,
                       high = color_delta_freq[10],
                       low = color_delta_freq[2],
                       guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))+
  labs(x="",y="",fill=bquote(Delta~Frequency),
       title = paste0(time_name,"-Century ",ssp))
pdf(paste0(Output_path,"Delta_frequency_maps_",time_name,"_",ssp,"_V2.pdf"),
    height = 4,
    width=8)
print(g)
dev.off()


  

