# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2025.7.26

# This code is to output 0.25D theta +- sd(Theta) maps

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

# Input path of combined 0.25D rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Figures/"

# Make CONUS boundary
# Whole US map
CONUS <- st_read("/fs/ess/PAS2204/Data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
# CONUS outer boundary map
#CONUS <- st_union(CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",])
CONUS <- CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",]

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


# ------ Main --------
# Convert theta back to slope
raster_all$Slope <- tan(raster_all$Theta/180*pi)
# Get Slope+sd and Slope-sd
raster_all$Slope_plussd <- raster_all$Slope + raster_all$Slope_sd
raster_all$Slope_minsd <- raster_all$Slope - raster_all$Slope_sd
# Convert them to theta
raster_all$Theta_plussd <- atan(raster_all$Slope_plussd)/pi*180
raster_all$Theta_minsd <- atan(raster_all$Slope_minsd)/pi*180

var_ls <- c("Theta_plussd","Theta_minsd")
my_color <- RColorBrewer::brewer.pal(11,"PRGn")
my_color <- c(my_color[2:5],my_color[7:10])

# Initialize a list to store both figures
g_all <- list()
for(i in 1:length(var_ls)){
  var_name <- var_ls[i]
  raster_df <- as.data.frame(raster_all[[var_name]],xy=TRUE)
  # Add two values to make sure the plot covers full color range
  raster_df[1,3] <- 90
  raster_df[2,3] <- -90
  g <- ggplot()+
    geom_tile(data=raster_df,aes(x=x,y=y,fill=.data[[var_name]]))+
    geom_sf(data=CONUS,fill=NA,color="black",alpha=0.5)+
    my_theme+
    labs(x="",y="",title = var_name,fill="")+
    scale_x_continuous(labels = scales::label_number(accuracy = 1)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 1))
  g <- g + 
    scale_fill_gradientn(na.value="white",
                         colors = my_color,
                         breaks = c(-90,-45,0,45,90),
                         guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))+
    labs(title = "",fill=bquote(alpha~ (degree)))+
    theme(legend.position = "bottom")
  g_all[[i]] <- g
}

# Combine the two plots
g <- plot_grid(plotlist = g_all,nrow=1,labels="auto")
# Output this figure
pdf(paste0(Output_path,"0.25D_210m_agg_Maps_theta_sd.pdf"),
    height=5,
    width=16)
print(g)
dev.off()


