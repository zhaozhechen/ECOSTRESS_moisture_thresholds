# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.6.19

# This code is to get frequency of CMIP6 projection that passes
# JOINT thresholds
# And average delta_VPD over the five years (VPD absolute value - VPD_star)
# For each of the 15 models

# Note: this code uses the absolute value of these climate conditions and
# Threshold values, rather than their quantiles
# This is to avoid the large amount of quantile 1 for future VPD

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
# 0.25D raster stacks, including slope, intercept, VPD_gs_mean, SM_gs_mean, VPD_max, and VPD_min
raster_stack <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final2/0.25D_intercepts/Combined_0.25D_Intercepts.rds")
# Output path of frequency and delta_VPD stacks
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks/"

# List of 15 models
models_ls <- c("CMCC-ESM2","CanESM5","CanESM5-1","EC-Earth3","INM-CM4-8",
               "INM-CM5-0","IPSL-CM6A-LR","KACE-1-0-G","MIROC6","MPI-ESM1-2-HR",
               "MPI-ESM1-2-LR","MRI-ESM2-0","NorESM2-LM","NorESM2-MM","TaiESM1")
arrayid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
model_name <- models_ls[arrayid]

############
# Functions
############
# This function is to get CMIP6 SM quantiles as in historical range
# Input is CMIP6 SM absolute values
Cal_SM_qt <- function(CMIP6_SM_abs){
  # Loop over each layer in CMIP6_SM_abs
  for(i in 1:nlayers(CMIP6_SM_abs)){
    values(CMIP6_SM_abs[[i]]) <- hist_range(values(CMIP6_SM_abs[[i]]))*100
  }
  return(CMIP6_SM_abs)
}

# This function is to calculate total frequency and delta_VPD over the five years
# Input file_name are the daily climate conditions to process
Cal_freq_all <- function(file_name){
  # Initialize a stack to store total number of days that pass thresholds
  freq_all <- stack()
  # Initialize a stack to store all delta_VPD
  delta_VPD_all <- stack()
  # Store total number of days for calculation
  n_day <- 0
  # Loop over files for this scenario
  for(i in 1:length(file_name)){
    # Read in climate data
    data <- readRDS(file_name[i])
    # Get VPD absolute values
    VPD_abs <- data$daily_maxVPD
    # Get SM quantiles (as in historical range)
    SM_qt <- Cal_SM_qt(data$daily_SM)
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
    # Get the average of delta_VPD
    delta_VPD_all <- stack(delta_VPD_all,mean(delta_VPD,na.rm=T))
  }
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
# SM historical range
hist_range <- ecdf(values(raster_stack$SM_gs_mean))

# Get frequency and delta_VPD for historical obs =====================================
# Files for historical observations
hist_obs_files <- paste0(hist_obs_path,dir(hist_obs_path))
# Historical obs does not have 2018.5 and 2018.6 SM, so skip
hist_obs_files <- hist_obs_files[3:25]
out_hist <- Cal_freq_all(hist_obs_files)
# Output this historical freq stack
saveRDS(out_hist,paste0(Output_path,"Hist_obs_freq.rds"))
print("Complete Historical obs")

# Get frequency and delta_VPD for the current CMIP6 model =========================
for(ssp in c("ssp245","ssp585")){
  for(time_name in c("Mid","End")){
    CMIP_files <- paste0(CMIP6_path,model_name,"_",time_name,"_",ssp,".rds")
    out_CMIP <- Cal_freq_all(CMIP_files)
    saveRDS(out_CMIP,paste0(Output_path,model_name,"_",time_name,"_",ssp,"_freq.rds"))
    print(paste("Complete",ssp,time_name))
  }
}
print("All done")


