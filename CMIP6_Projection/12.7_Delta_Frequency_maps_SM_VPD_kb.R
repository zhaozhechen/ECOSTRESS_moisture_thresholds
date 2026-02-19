# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.9.17

# This code is to plot delta_freqeucny
# Future projections - historical frequency
# For (1) Only change SM
# (2) Only change VPD
# For two kb combinations: kb1 and kb2
# Masked by p-value rasters (only keep p-value <= 0.05)

# Note: k1 was calculated using mean slope + sd
# k2 was calculated using mean slope - sd

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

# Projected frequency changing only SM or VPD, calculated using two kb combinations
Freq_proj_SM_VPD <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_SM_VPD_kb/Model_mean_End_ssp585_freq_SM_VPD_kb.rds")
# Projected frequency changing both SM and VPD, calculated using two kb combinations
Freq_proj_both <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_kb/Model_mean_End_ssp585_freq_kb.rds")
# Input path of Historical frequency calculated using two kb combinations
Freq_hist_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_kb/"
# Input path of p-value rasters
p_value_stack <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_SM_VPD_kb/p_value_rasters.rds")
# Input path of combined 0.25D rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")
# Make CONUS boundary
# Whole US map
CONUS <- st_read("/fs/ess/PAS2204/Data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
# CONUS outer boundary map
#CONUS <- st_union(CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",])
CONUS <- CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",]
# Output path for the figure
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Figures_V2/"

# This is the color for delta frequency
color_delta_freq <- rev(RColorBrewer::brewer.pal(11,"PiYG"))[c(2,4,6,8,10)]

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

# This function is to make delta_frequency maps
# Input raster is the raster to plot
# my_color is the color to use
# my_title is the title of the plot
g_delta_map <- function(raster,my_color,my_title){
  df <- as.data.frame(raster,xy=TRUE)
  df$layer[df$layer >=0.4] <- 0.4
  df$layer[df$layer <= -0.4] <- -0.4
  g <- ggplot()+
    geom_tile(data=df,aes(x=x,y=y,fill=layer))+
    geom_sf(data=CONUS,fill=NA,color="black",alpha=0.5)+
    my_theme+
    scale_fill_gradientn(na.value="white",
                         colours = my_color,
                         guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"),
                         breaks = c(-0.4,-0.2,0,0.2,0.4),
                         labels = c(expression(""<= -0.4),-0.2,0,0.2,expression("">=0.4)))+
    labs(x="",y="",fill=bquote(Delta~Frequency),
         title = paste0(my_title))+
    theme(legend.position = "right") 
  
  return(g)
}

#######
# Main
#######
# Read in data and pre-processing =======================
time_name <- "End"
ssp <- "ssp585"

# Rename Projected frequency for changing both SM and VPD
names(Freq_proj_both) <- c("kb1","kb2")

# Initialize a list to store all maps
p_ls <- list()
for(var_chg in c("SM_","VPD_","")){
  for(kb in c(1,2)){
    # Historical frequency with this kb combination
    Freq_hist <- readRDS(paste0(Freq_hist_path,"Historical_freq_kb",kb,".rds"))
    
    # Projected frequency with this kb combination, for this var_chg scenario
    if(var_chg == ""){
      Freq_proj_tmp <- Freq_proj_both
    }else{
      Freq_proj_tmp <- Freq_proj_SM_VPD
    }
    # Get the corresponding raster layer for the projected frequency
    Freq_proj_tmp <- Freq_proj_tmp[[paste0(var_chg,"kb",kb)]]
    
    # Calculate delta_frequency
    delta_freq <- Freq_proj_tmp - Freq_hist
    # Mask by theta
    delta_freq <- mask(delta_freq,raster_all$Theta)
    
    # Get the corresponding p-value raster
    if(var_chg == ""){
      p_name <- paste0("p_value_",time_name,"_",ssp,"_kb",kb)
    }else{
      p_name <- paste0("p_value_",time_name,"_",ssp,"_",var_chg,"chg_kb",kb)
    }
    p_value <- p_value_stack[[p_name]]
    # Only keep p_value <= 0.05
    p_value[p_value>0.05] <- NA
    # Mask delta_frequency by p-value
    delta_freq <- mask(delta_freq,p_value)

    # Add one pixel of value to make the color scales consistent
    delta_freq[1] <- 0.4
    delta_freq[2] <- -0.4
    # Make map
    if(kb == 1){title <- expression(alpha + sd)}
    delta_freq_map <- g_delta_map(delta_freq,color_delta_freq,paste("Changing",var_chg,kb))  
    # Store this map
    p_ls[[length(p_ls)+1]] <- delta_freq_map
  }
}

# Combine these plots
g <- plot_grid(plotlist = p_ls,nrow=3)
# Output the figure
pdf(paste0(Output_path,"delta_Freuquency_maps_",time_name,"_",ssp,"_SM_VPD_kb_pvalue_masked.pdf"),
    height = 8,
    width = 12)
print(g)
dev.off()