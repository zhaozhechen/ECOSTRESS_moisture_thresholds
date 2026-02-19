# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.8.11

# The ultimate goal is to quantify frequency uncertainty caused by alpha uncertainty
# In this code, we first get combinations of slope and intercept
# At each pixel, change slope to slope +- sd, and fix one point: ((50+VPD_50th)/2,(50+SM_50th)/2)
# And calculate the corresponding intercept, as b1 and b2

# Output a map of intercept and slope combinations

###########
# Global
###########
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)

library(raster)
library(ncdf4)
library(sf)
library(dplyr)

# Input path of combined 0.25D rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")
# Output path
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/"
########
# Main
######## 

# Get Slope from theta, which is the mean slope
slope_mean <- tan(raster_all$Theta/180*pi)
# Fixed point on each GS filtration WLS line
# Should be ((50+VPD_50th)/2,(50+SM_50th)/2)
fix_x <- (50+raster_all$VPD_50th)/2
fix_y <- (50+raster_all$SM_50th)/2
# Change slope to slope+slope_sd and slope-slope_sd
k1 <- slope_mean + raster_all$Slope_sd
k2 <- slope_mean - raster_all$Slope_sd

# Calculate corresponding intercept
b1 <- fix_y - k1*fix_x
b2 <- fix_y - k2*fix_x

# Output k1,b1,k2,b2
out <- stack(k1,b1,k2,b2)  
names(out) <- c("k1","b1","k2","b2")
saveRDS(out,paste0(Output_path,"Combined_k_b.rds"))
  

  

