# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2025.6.22

# This code is to output 0.25D threshold maps (aggregated from 210m)
# Make absolute values of SM50 - SM25, SM75 - SM50, and VPD50 - VPD25, and VPD75 - VPD50
# Note: the variable names are the opposite

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
  legend.position = "bottom",
  plot.title = element_text(size=16),
  legend.title = element_text(size=16),
  #plot.margin = margin(80,80,80,80),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank()
  #legend.margin = margin(0,80,0,80)
  #axis.ticks.length = unit(-0.8,"cm")
)

# This function is to plot change in moisture thresholds
# The input are the two variable names
# This function will plot var1 - var2
plot_delta_thresholds <- function(var1_name, var2_name){
  raster_df <- as.data.frame(raster_all[[c(var1_name,var2_name)]],xy=TRUE)
  # Get the target variable
  raster_df <- raster_df %>%
    mutate(var = .data[[var1_name]] - .data[[var2_name]])
  
  g <- ggplot()+
    geom_tile(data=raster_df,aes(x=x,y=y,fill=var))+
    geom_sf(data=CONUS,fill=NA,color="black",alpha=0.5)+
    my_theme+
    labs(x="",y="",title = "",fill="")+
    scale_x_continuous(labels = scales::label_number(accuracy = 1)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 1))
  
  # Use different color scales for different variables
  if(grepl("VPD",var1_name)){
    my_color <- RColorBrewer::brewer.pal(11, "BrBG")
    my_color <- c(my_color[2:5],my_color[7:10])
    g <- g +
      scale_fill_gradientn(na.value = "white",
                           colours = my_color,
                           limits = c(-0.2,0.2),
                           # values out of bounds are limited to the bounds
                           oob = scales::squish,
                           guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))
  }else{
    my_color <- wes_palette("Zissou1",100,type="continuous")
    g <- g +
      scale_fill_gradientn(na.value = "white",
                           colours = my_color,
                           limits = c(-2,2),
                           # values out of bounds are limited to the bounds
                           oob = scales::squish,
                           guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))
  }
  
  return(g)  
}

########
# Main
########

# Mask all layers with Theta
raster_all <- mask(raster_all,mask = raster_all$Theta)

# -----------
# SM25 - SM50 in absolute values
g1 <- plot_delta_thresholds("VPD_25_abs","VPD_50_abs")
# SM75 - SM50 in absolute values
g2 <- plot_delta_thresholds("VPD_75_abs","VPD_50_abs")

# VPD25 - VPD50 in absolute values
g3 <- plot_delta_thresholds("SM_25_abs","SM_50_abs")
# VPD75 - VPD50 in absolute values
g4 <- plot_delta_thresholds("SM_75_abs","SM_50_abs")

# Combine four maps
g <- plot_grid(g1,g2,g3,g4,
               nrow=2,
               labels = "auto",
               align = "hv",
               axis="lrbt")

pdf(paste0(Output_path,"0.25D_210m_agg_delta_thresholds.pdf"),
    height=10,
    width=16)
print(g)
dev.off()

# Calculate statistics -------------
# Calculate absolute SM threshold changes when VPD changes from 50 to 75 quantile
df <- as.data.frame(raster_all[[c("VPD_50_abs","VPD_75_abs","AI")]])
df <- df %>%
  mutate(dif = VPD_75_abs - VPD_50_abs,
         percent = dif/VPD_50_abs)
# Group by aridity category
# Get dryland type classification
df$Type[df$AI < 0.2] <- "1"
df$Type[df$AI >= 0.2 & df$AI < 0.5] <- "2"
df$Type[df$AI >= 0.5] <- "3"
df$Type <- as.factor(df$Type)
levels(df$Type) <- c("Arid","Semiarid","Semihumid_and_humid")

df_summary <- df %>%
  group_by(Type) %>%
  summarise(mean_percent = mean(percent,na.rm=TRUE))

