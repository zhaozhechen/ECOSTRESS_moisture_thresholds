# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.6.26

# This code is to test significant difference between frequency from historical climate
# And frequency from 15 models (future projection, 15 values) for each pixel

#########
# Global
#########

dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
library(raster)
library(ncdf4)
library(sf)
library(dplyr)
library(ggplot2)
library(cowplot)

# Path to the frequency stacks
Freq_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_V2/"
# Make CONUS boundary
# Whole US map
CONUS <- st_read("/fs/ess/PAS2204/Data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
# CONUS outer boundary map
#CONUS <- st_union(CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",])
CONUS <- CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",]
# Output path for plots
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Figures_V2/"

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

#######
# Main
#######

# Get p-value rasters =============================================================
# Historical frequency
hist_freq <- readRDS(paste0(Freq_path,"Hist_obs_freq.rds"))$Freq
hist_freq[hist_freq == 1] <- NA

# Initialize a stack to store p-value rasters for all four scenarios
p_value_stack <- stack()
for(time_name in c("Mid","End")){
  for(ssp in c("ssp245","ssp585")){
    # Get all file names for this scenario (they are results from the 15 models)
    file_names <- dir(Freq_path)[grepl(paste0(time_name,"_",ssp,"_freq"),dir(Freq_path))]
    # Stack these frequencies together
    CMIP_freq <- stack()
    for(i in 1:length(file_names)){
      CMIP_freq_tmp <- readRDS(paste0(Freq_path,file_names[i]))$Freq
      CMIP_freq <- stack(CMIP_freq,CMIP_freq_tmp)
    }
    print("Complete stacking CMIP frequencies")
    
    # Compare the historical frequency vs 15 CMIP frequencies for each pixel
    # Using Wilcox Signed-rank test
    # Initialize a vector to store all p-values
    p_all <- c()
    for(i in 1:ncell(hist_freq)){
      # Get 15 models' frequency for this pixel
      CMIP_freq_px <- values(CMIP_freq)[i,]
      if(is.na(hist_freq[i])){
        # If no historical frequency,p-value = NA
        p_tmp <- NA
      }else{
        # Do Wilcox test
        p_tmp <- wilcox.test(CMIP_freq_px,mu= values(hist_freq)[i])$p.value
      }
      p_all <- c(p_all,p_tmp)
    }
    
    # Construct a raster for p-value
    p_raster <- hist_freq
    values(p_raster) <- p_all
    names(p_raster) <- paste0("p_value_",time_name,"_",ssp)
    # Store this p-value raster
    p_value_stack <- stack(p_value_stack,p_raster)
    print(paste("Complete",time_name,ssp))
  }
}

# Output p-value rasters
saveRDS(p_value_stack,paste0(Freq_path,"p_value_rasters.rds"))

# Make pdf and map of p-values ===============================
# Initialize a list to store figures
p_ls <- list()
for(i in 1:nlayers(p_value_stack)){
  # Make pdf
  p_df <- as.data.frame(p_value_stack[[i]],xy=TRUE)
  p_pdf <- ggplot(p_df,aes(x=p_df[,3],after_stat(density)))+
    stat_density(position="identity",geom="line")+
    labs(x="p-value")+
    pdf_theme
  # Make map
  # Assign colors for p-values
  fill <- rep(NA,nrow(p_df))
  fill[p_df[,3] <= 0.05] <- "a"
  fill[p_df[,3] > 0.05 & p_df[,3] <= 0.1] <- "b"
  fill[p_df[,3] > 0.1] <- "c"
  p_df$fill <- fill
  
  p_map <- ggplot(p_df)+
    geom_tile(data=p_df,aes(x=x,y=y,fill=fill))+
    geom_sf(data=CONUS,fill=NA,color="black",alpha=0.5)+
    my_theme+
    scale_fill_discrete(breaks = c("a","b","c"),
                        labels = c(bquote(p-value <= 0.05),
                                   bquote(p-value <= 0.1),
                                   bquote(p-value > 0.1)),
                        na.value = "white")+
    labs(x="",y="",fill="",
         title = names(p_df)[3])
  # Combine p_map and p_pdf
  p <- plot_grid(p_map,p_pdf,nrow=1,
                 rel_widths = c(1.5,1))
  # Store this p
  p_ls[[i]] <- p
  print(i)
}

# Combine figure for all four scenarios
p <- plot_grid(plotlist = p_ls,
               nrow=4)
pdf(paste0(Output_path,"p_value_maps.pdf"),
    height =16,width=12)
print(p)
dev.off()
  




