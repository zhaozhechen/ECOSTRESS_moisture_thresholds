# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.7.10

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
Input_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_V2/"
# Input path of combined 0.25D rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")
# Make CONUS boundary
# Whole US map
CONUS <- st_read("/fs/ess/PAS2204/Data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
# CONUS outer boundary map
#CONUS <- st_union(CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",])
CONUS <- CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",]
# Output path
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Figures_V2/"

# This is the color for frequency
color_freq <- RColorBrewer::brewer.pal(9,"YlOrRd")
# This is the color for delta frequency
color_delta_freq <- rev(RColorBrewer::brewer.pal(11,"PiYG"))
# This is the color for pdf of frequency
pdf_color <- c("#8DD3C7","#FB8072")
# This is the color for delta_VPD
#color_delta_VPD <- wes_palette("Zissou1",100,type="continuous")
color_delta_VPD <- RColorBrewer::brewer.pal(9,"OrRd")
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
# Read in data and pre-processing ============================================
time_name <- "End"
ssp <- "ssp585"
# Historical frequency
Hist_freq <- readRDS(paste0(Input_path,"Hist_obs_freq.rds"))$Freq
Hist_freq[Hist_freq == 1] <- NA
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

# Make maps =====================================================
# Historical frequency
Hist_map <- g_map(Hist_freq,color_freq,
                  my_legend = "Frequency",my_title = "Historical",
                  min = 0,max=1)
# Projected frequency
Proj_map <- g_map(Proj_freq,color_freq,
                  my_legend = "Frequency",my_title = paste(time_name,ssp),
                  min = 0,max=1)
# Make map of delta_frequency, masked by p-value > 0.05
df <- as.data.frame(delta_freq,xy=TRUE)
df$p_value <- values(p_value)
df$p_value[df$p_value <= 0.05] <- NA

df$layer[df$layer >= 0.4] <- 0.4
#df$layer[df$layer <= (-0.2)] <- (-0.2)

Delta_map <- ggplot()+
  geom_tile(data=df,aes(x=x,y=y,fill=layer))+
  geom_sf(data=CONUS,fill=NA,color="black",alpha=0.5)+
  geom_point(data=df[!is.na(df$p_value),],
             aes(x=x,y=y),
             shape = 16,
             fill = "black",
             size = 0.5,
             alpha = 0.3)+
  my_theme+
  scale_fill_gradientn(na.value="white",
                       colours = color_delta_freq[c(2,6,7,8,9,10)],
                       guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"),
                       breaks = c(-0.1,0,0.1,0.2,0.3,0.4),
                       labels = c(-0.1,0,0.1,0.2,0.3,expression("">=0.4)))+
  labs(x="",y="",fill=bquote(Delta~Frequency),
       title = paste0(time_name,"-Century ",ssp))

# pdf of frequency ==============================================
freq_df <- data.frame(Frequency = c(values(Hist_freq),values(Proj_freq)),
                      Type = rep(c("Historical","Projected"),each = 24000))
hist_mean <- mean(values(Hist_freq),na.rm=T)
proj_mean <- mean(values(Proj_freq),na.rm=T)

g_pdf <- ggplot(freq_df,aes(x=Frequency,color=Type))+
  geom_density()+
  geom_vline(xintercept = hist_mean,linetype = 2,color=pdf_color[1])+
  geom_vline(xintercept = proj_mean,linetype = 2,color=pdf_color[2])+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA),
        legend.key = element_blank(),
        #legend.key.size = unit(6,"cm"),
        legend.text = element_text(size=14),
        plot.title = element_text(size=16),
        axis.text = element_text(size=16),
        legend.position = c(0.8,0.8),
        axis.title = element_text(size=16),
        aspect.ratio = 1/1.5)+
  scale_color_manual(values=c(Historical = pdf_color[1],Projected = pdf_color[2]))+
  labs(x= "Frequency",color="")
  
# Combine these four plots
g <- plot_grid(Hist_map,Proj_map,Delta_map,g_pdf,
               nrow=2,labels = "auto")
pdf(paste0(Output_path,"Frequency_maps_",time_name,"_",ssp,".pdf"),
    height=8,
    width=15)
print(g)
dev.off()

# Maps of delta_VPD thresholds ===================================================
# Historical delta_VPD
Hist_delta_VPD <- readRDS(paste0(Input_path,"Hist_obs_freq.rds"))$delta_VPD
Hist_delta_VPD[Hist_delta_VPD>6] <- 6
# Map of historical delta_VPD
g1 <- g_map(Hist_delta_VPD,color_delta_VPD,bquote(Delta~Thresholds[VPD](kPa)),
            "Historical",NA,NA)

# Projected delta_VPD
Proj_delta_VPD <- readRDS(paste0(Input_path,"Model_mean_freq_",time_name,"_",ssp,".rds"))$Delta_VPD_mean
Proj_delta_VPD[Proj_delta_VPD>6] <- 6
# Map of projected delta_VPD
g2 <- g_map(Proj_delta_VPD,color_delta_VPD,bquote(Delta~Thresholds[VPD](kPa)),
            paste0(time_name,"-Century ",ssp),NA,NA)
# Combine the two plots
g <- plot_grid(g1,g2,nrow=1)
pdf(paste0(Output_path,"delta_VPD_maps_",time_name,"_",ssp,".pdf"),
    height=4,
    width=15)
print(g)
dev.off()





