# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.7.10

# This code is to get frequency of Historical obs that passes
# JOINT thresholds
# And average delta_VPD over the five years (VPD absolute value - VPD_star)
# For each of the 15 models

# When calculating frequency for each pixel, get the SM quantile within its own SM range,
# Then calculate VPD absolute values at this SM quantile

# Note: this code uses the absolute value of these climate conditions and
# Threshold values, rather than their quantiles
# This is to avoid the large amount of quantile 1 for future VPD

# Delta_VPD was calculation used average of only positive delta_VPD

#########
# Global
#########

dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
library(raster)
library(ncdf4)
library(sf)
library(dplyr)

# Path to daily adjusted CMIP6 projection
CMIP6_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/CMIP6_adjusted/"
# Input path of daily historical observations
hist_obs_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final2/Hist_obs/"
# Historical SM range list (for each pixel)
Hist_SM_range_ls <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Hist_obs/Hist_SM_range_ls.rds")
# 0.25D raster stacks, including slope, intercept, VPD_gs_mean, SM_gs_mean, VPD_max, and VPD_min
raster_stack <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final2/0.25D_intercepts/Combined_0.25D_Intercepts.rds")
# Output path of frequency and delta_VPD stacks
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_V2/"

############
# Functions
############
# This function is to get CMIP6 SM quantiles as in historical range
# Input is CMIP6 SM absolute values
Cal_SM_qt <- function(CMIP6_SM_abs){
  # Loop over each pixel
  for(i in 1:ncell(CMIP6_SM_abs)){
    if(all(!is.na(values(CMIP6_SM_abs)[i,]))){
      # Convert all daily SM for this pixel from absolute values to quantiles as in historical range
      values(CMIP6_SM_abs)[i,] <- Hist_SM_range_ls[[i]](values(CMIP6_SM_abs)[i,])*100
    }
    if(i %% 1000 == 0){
      print(paste("Complete",i,"out of",ncell(CMIP6_SM_abs),"cells")) 
    }
  }
  return(CMIP6_SM_abs)  
}

# This function is to calculate total frequency and delta_VPD over the five years
# For historical observations
# Input are the daily climate conditions to process
Cal_freq_all <- function(VPD_abs,SM_abs){
  # Initialize a stack to store total number of days that pass thresholds
  freq_all <- stack()
  # Initialize a stack to store all delta_VPD
  delta_VPD_all <- stack()
  # Store total number of days for calculation
  n_day <- 0
  
  # Get SM quantiles (as in historical range)
  SM_qt <- Cal_SM_qt(SM_abs)
  # Record number of days
  n_day <- n_day + nlayers(VPD_abs)
  # Use the SM quantiles to calculate the corresponding VPD_star in a GS plot as quantiles
  VPD_star <- SM_qt * raster_stack$Slope + raster_stack$Intercept
  # Convert VPD star to its original scale
  VPD_star <- VPD_star/100 * (raster_stack$VPD_max - raster_stack$VPD_min) + raster_stack$VPD_min
  # Compare if VPD absolute value is greater than VPD_star
  freq <- VPD_abs > VPD_star
  # Sum them up, as the number of days that VPD absolute values that pass VPD_star
  freq <- calc(freq,fun = sum,na.rm=T)
  freq_all <- stack(freq_all,freq)
  # Calculate the difference between CMIP6 VPD absolute value and VPD_star
  delta_VPD <- VPD_abs - VPD_star
  # Only keep positive delta_VPD
  delta_VPD[delta_VPD <0] <- NA
  # Get the average of delta_VPD
  delta_VPD_all <- stack(delta_VPD_all,mean(delta_VPD,na.rm=T))
  # Get the average of delta_VPD over the five years
  delta_VPD <- mean(delta_VPD_all,na.rm=T)
  delta_VPD <- mask(delta_VPD,raster_stack$Slope)
  if(nlayers(freq_all)>1){
    freq <- calc(freq_all,fun=sum,na.rm=T)/n_day
  }else{
    freq <- freq_all/n_day
  }
  freq <- mask(freq,raster_stack$Slope)
  out <- stack(freq,delta_VPD)
  names(out) <- c("Freq","delta_VPD")
  return(out)
}

######
# Main
######
# Combine all historical obs files into one raster stack, for daily_max_VPD and daily_SM, respectively
# Files for historical observations
hist_obs_files <- paste0(hist_obs_path,dir(hist_obs_path))
# Historical obs does not have 2018.5 and 2018.6 SM, so skip
hist_obs_files <- hist_obs_files[3:25]
VPD_abs <- stack()
SM_abs <- stack()
for(i in 1:length(hist_obs_files)){
  data <- readRDS(hist_obs_files[i])
  # Stack all daily_max_VPD
  VPD_abs <- stack(VPD_abs,data$daily_maxVPD)
  # Stack all daily_SM
  SM_abs <- stack(SM_abs,data$daily_SM)
}

VPD_abs <- mask(VPD_abs,raster_stack$Slope)
SM_abs <- mask(SM_abs,raster_stack$Slope)

# Get frequency and delta_VPD for historical obs =====================================
out_hist <- Cal_freq_all(VPD_abs,SM_abs)
# Output this historical freq stack
saveRDS(out_hist,paste0(Output_path,"Hist_obs_freq.rds"))
print("Complete Historical obs")
