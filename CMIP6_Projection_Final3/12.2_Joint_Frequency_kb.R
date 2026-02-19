# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.8.11

# This code is to get frequency of daily climate conditions that passes JOINT THRESHOLDS
# Using new combinations of k and b 
# For historical obs and CMIP6 End ssp585

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
# New k and b combinations
kb_raster <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Combined_k_b.rds")
# Output path
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_kb/"

# List of 15 models
models_ls <- c("CMCC-ESM2","CanESM5","CanESM5-1","EC-Earth3","INM-CM4-8",
               "INM-CM5-0","IPSL-CM6A-LR","KACE-1-0-G","MIROC6","MPI-ESM1-2-HR",
               "MPI-ESM1-2-LR","MRI-ESM2-0","NorESM2-LM","NorESM2-MM","TaiESM1")

time_name <- "End"
ssp <- "ssp585"

arrayid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# The first digit is for var_ID
var_ID <- as.numeric(substr(arrayid,1,1))
if(var_ID == 1){
  k <- kb_raster$k1
  b <- kb_raster$b1
}else if(var_ID == 2){
  k <- kb_raster$k2
  b <- kb_raster$b2
}
# The second and third digit is for model
job_ID <- as.numeric(substr(arrayid,2,3))
if(job_ID < 16){
  model_name <- models_ls[job_ID]
}

############
# Functions
############
# This function is to get CMIP6 SM quantiles as in historical range
# Input is CMIP6 SM absolute values (could also be historical SM absolute values)
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

# This function is to calculate total frequency
# Under the given VPD and SM combinations
# Input are the daily climate conditions to process
# k and b are slope and intercept to use
Cal_freq_all <- function(VPD_abs,SM_abs,k,b){
  # Record number of days
  n_day <- nlayers(VPD_abs)
  # Get SM quantiles (as in historical range)
  SM_qt <- Cal_SM_qt(SM_abs)
  # Use the SM quantiles to calculate the corresponding VPD_star in a GS plot as quantiles
  VPD_star <- SM_qt * k + b
  # Convert VPD star to its original scale
  VPD_star <- VPD_star/100 * (raster_stack$VPD_max - raster_stack$VPD_min) + raster_stack$VPD_min
  # Compare if VPD absolute value is greater than VPD_star
  freq <- VPD_abs > VPD_star
  # Sum them up, as the number of days that VPD absolute values that pass VPD_star
  freq <- calc(freq,fun = sum,na.rm=T)
  freq <- freq/n_day
  freq <- mask(freq,raster_stack$Slope)
  return(freq)
}

########
# Main
########
# If Job_ID is 16, use historical VPD and SM
if(job_ID == 16){
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
  # Calculate frequency
  freq <- Cal_freq_all(VPD_abs,SM_abs,k,b)
  saveRDS(freq,paste0(Output_path,"Historical_freq_kb",var_ID,".rds"))
# If Job_ID is 1-15, get CMIP6 projected VPD and SM for 15 models
}else{
  CMIP_files <- paste0(CMIP6_path,model_name,"_",time_name,"_",ssp,".rds")
  # Use CMIP6 VPD
  VPD_abs <- readRDS(CMIP_files)$daily_maxVPD
  # Use CMIP6 SM
  SM_abs <- readRDS(CMIP_files)$daily_SM
  VPD_abs <- mask(VPD_abs,raster_stack$Slope)
  SM_abs <- mask(SM_abs,raster_stack$Slope)
  # Calculate frequency
  freq <- Cal_freq_all(VPD_abs,SM_abs,k,b)
  saveRDS(freq,paste0(Output_path,model_name,"_",time_name,"_",ssp,"_freq_kb",var_ID,".rds"))
}

print("All done!!")


