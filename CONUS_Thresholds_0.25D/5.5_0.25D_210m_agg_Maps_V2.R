# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.4.15

# This code is to output 0.25D thresholds maps (aggregated from 210m)
# Revise panel b and c to VPD50 and SM50 quantiles instead of absolute values

#########
# Global
#########
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)

library(raster)
library(sf)
library(ggplot2)
library(dplyr)
library(cowplot)
library(wesanderson)

# Input path of combined 0.25D rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Figures/"

# Make CONUS boundary
# Whole US map
CONUS <- st_read("/fs/ess/PAS2204/Data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
# CONUS outer boundary map
#CONUS <- st_union(CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",])
CONUS <- CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",]

# Variable list to plot
var_ls <- c("Theta","SM_50th","VPD_50th")

###########
# Functions
###########
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

#######
# Main
#######

# Mask all layers with Theta
raster_all <- mask(raster_all,mask = raster_all$Theta)

# Initialize a list to store all maps
g_all <- list()
for(i in 1:length(var_ls)){
  var_name <- var_ls[i]
  raster_df <- as.data.frame(raster_all[[var_name]],xy=TRUE)
  g <- ggplot()+
    geom_tile(data=raster_df,aes(x=x,y=y,fill=.data[[var_name]]))+
    geom_sf(data=CONUS,fill=NA,color="black",alpha=0.5)+
    my_theme+
    labs(x="",y="",title = var_name,fill="")+
    scale_x_continuous(labels = scales::label_number(accuracy = 1)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 1))
  # Use different color scales for different variables
  if(var_name == "Theta"){
    my_color <- RColorBrewer::brewer.pal(11,"PRGn")
    my_color <- c(my_color[2:5],my_color[7:10])
    g <- g + 
      scale_fill_gradientn(na.value="white",
                          colors = my_color,
                          breaks = c(-10,20,50,80),
                          guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))+
      labs(title = "",fill=bquote(alpha~ (degree)))+
      theme(legend.position = "bottom")
  }else if(var_name == "VPD_50th"){
    BrBG_color <- RColorBrewer::brewer.pal(11, "BrBG")
    BrBG_color <- c(BrBG_color[2:5],BrBG_color[7:10])
    g <- g + 
      scale_fill_gradientn(na.value="white",
                           colours = BrBG_color,
                           breaks = c(25,50,75),
                           guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))+
      labs(title = "",fill = bquote(SM[50]))+
      theme(legend.position = "bottom")
  }else if(var_name == "SM_50th"){
    #my_color <- RColorBrewer::brewer.pal(11,"RdYlGn")
    my_color <- wes_palette("Zissou1",100,type="continuous")
    #my_color <- rev(c(my_color[2:5],my_color[7:10]))
    g <- g + 
      scale_fill_gradientn(na.value="white",
                          colours = my_color,
                          breaks = c(25,50,75),
                          guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))+
      labs(title = "",fill = bquote(VPD[50]))+
      theme(legend.position = "bottom")
  }
  
  # Add this figure to the figure list
  g_all[[i]] <- g
  print(paste0("Complete Map",i))
}

# Combine three maps
g <- plot_grid(plotlist = g_all,nrow=2,
               labels="auto",
               align = "hv",
               axis="lrbt")

pdf(paste0(Output_path,"0.25D_210m_agg_Maps_V3.pdf"),
    height=10,
    width=16)
print(g)
dev.off()


