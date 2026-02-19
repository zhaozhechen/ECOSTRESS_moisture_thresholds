# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.6.18

# This code is to get adjusted CMIP6 projection for the two periods:
# Mid-century: 2046-2050
# End-century: 2096-2100
# For the two ssp: ssp245 and ssp585

# Methods to adjust:
# Get one map of 5-year growing season average for historical observations 
# Get one map of 5-year growing season average for CMIP6 historical
# Get the delta of these two, then used to adjust CMIP projection
# CMIP6_proj_adj (daily) = CMIP6_proj (daily) - (CMIP6_hist (5-year gs mean) - hist_obs (5-year gs mean))

# The output are filtered by T>0
# Also, VPD and SM <0 are adjusted to 0
# Also output a df of CONUS daily average

#########
# Global
#########
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
library(raster)
library(ncdf4)
library(sf)

# Input path of historical observations
hist_obs_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final2/Hist_obs/"
# Input path of CMIP6 15 model raster stacks
CMIP_raw_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/CMIP6_raw/"
# Input path of combined 0.25D rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")
# List of 15 models
models_ls <- c("CMCC-ESM2","CanESM5","CanESM5-1","EC-Earth3","INM-CM4-8",
               "INM-CM5-0","IPSL-CM6A-LR","KACE-1-0-G","MIROC6","MPI-ESM1-2-HR",
               "MPI-ESM1-2-LR","MRI-ESM2-0","NorESM2-LM","NorESM2-MM","TaiESM1")
# Output path for historical obs
Hist_Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Hist_obs/"
# Output path for adjusted CMIP
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/CMIP6_adjusted/"
# Determines which model to process
arrayid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
model_name <- models_ls[arrayid]

###########
# Functions
###########

# This function is to get the 5-year gs mean and daily CONUS mean
# The input file names determines which files to process
# varname determines which variable to process
gs_mean <- function(file_name,varname){
  # Initialize a list to store both gs_mean and CONUS daily
  out <- list()
  # Initialize a vector to store daily CONUS mean
  daily_CONUS <- c()
  # Initialize a stack to store all monthly(yearly) mean for this variable
  long_term_all <- stack()
  
  # Loop over all files within this file_name
  for(i in 1:length(file_name)){
    # Read in the raster stack list
    stack_all <- readRDS(file_name[i])
    # Only keep the target variable
    var_stack <- stack_all[[varname]]
    
    # Convert unit from K to degree C for daily T (only historical obs used K for unit)
    if(varname == "daily_T"){
      var_stack <- var_stack - 273.15
    }
    
    # Record the daily average of CONUS, this is for visualization of climate pdf,
    # Not directly related to CMIP6 projection adjustment
    for(layer in 1:nlayers(stack_all[[2]])){
      # If there is no daily_SM or only one layer, give NA
      if(length(var_stack)<=24001){
        daily_CONUS_tmp <- NA
      }else{
        daily_CONUS_tmp <- mean(values(var_stack[[layer]]),na.rm=T)
      }
      daily_CONUS <- c(daily_CONUS,daily_CONUS_tmp)
    }
    
    # Get the mean of this file
    var_stack <- mean(var_stack)
    # If this raster is NA, ignore, else store it to the 5-year gs stack
    if(length(var_stack)==24000){
      long_term_all <- stack(long_term_all,var_stack)
    }
  }
  # Get the 5-year gs mean
  long_term_all <- mean(long_term_all)
  # Mask it with theta
  long_term_all <- mask(long_term_all,raster_all$Theta)
  out[[1]] <- long_term_all
  out[[2]] <- daily_CONUS
  return(out)
}

#######
# Main
#######
# Get gs mean of historical observations, and delta between gs mean of historical obs and historical CMIP ==================
# List of variables to process
varname_ls <- c("daily_SM","daily_maxVPD","daily_T")
# Files for historical observations
hist_obs_files <- paste0(hist_obs_path,dir(hist_obs_path))
# Files for historical CMIP6
hist_245_files <- paste0(CMIP_raw_path,model_name,"_",2018:2022,"_ssp245.rds")
hist_585_files <- paste0(CMIP_raw_path,model_name,"_",2018:2022,"_ssp585.rds")

# Initialize a stack to store 5-year gs mean for all target variables for historical obs
hist_obs_gs_mean <- stack()
# Initialize a df to store daily CONUS mean for all target variables
hist_obs_daily_CONUS_mean <- c()
# Initialize a stack to store difference between historical CMIP gs mean and historical observation gs mean
hist_CMIP_245_gs_mean_dif <- stack()
hist_CMIP_585_gs_mean_dif <- stack()

# Loop over all target variables
for(i in 1:length(varname_ls)){
  varname <- varname_ls[i]
  # Get gs mean for historical observations
  hist_obs_gs_mean_tmp <- gs_mean(hist_obs_files,varname)[[1]]
  # Get the difference between gs mean for historical CMIP6 and gs mean for historical observations
  hist_CMIP_245_gs_mean_dif <- stack(hist_CMIP_245_gs_mean_dif,gs_mean(hist_245_files,paste0(varname,"_ssp245"))[[1]] -
                                       hist_obs_gs_mean_tmp)
  hist_CMIP_585_gs_mean_dif <- stack(hist_CMIP_585_gs_mean_dif,gs_mean(hist_585_files,paste0(varname,"_ssp585"))[[1]] -
                                       hist_obs_gs_mean_tmp)
  
  hist_obs_gs_mean <- stack(hist_obs_gs_mean,hist_obs_gs_mean_tmp)
  # Get daily CONUS mean for historical observations
  hist_obs_daily_CONUS_mean <- cbind(hist_obs_daily_CONUS_mean,gs_mean(hist_obs_files,varname)[[2]])
  print(paste("Complete",varname))
}

names(hist_obs_gs_mean) <- varname_ls
# Output this raster stack for plots
saveRDS(hist_obs_gs_mean,paste0(Hist_Output_path,"Hist_obs_gs_mean.rds"))

colnames(hist_obs_daily_CONUS_mean) <- varname_ls
# Output this df for plots
write.csv(hist_obs_daily_CONUS_mean,paste0(Hist_Output_path,"Hist_obs_daily_CONUS_mean.csv"))

hist_CMIP_gs_mean_dif <- stack(hist_CMIP_245_gs_mean_dif,hist_CMIP_585_gs_mean_dif)
names(hist_CMIP_gs_mean_dif) <- c(paste0(varname_ls,"_ssp245"),paste0(varname_ls,"_ssp585"))
rm(hist_CMIP_245_gs_mean_dif,hist_CMIP_585_gs_mean_dif)

# Adjust CMIP projection ================================================================
for(time_name in c("Mid","End")){
  if(time_name == "Mid"){
    # Mid-century
    time_ls <- 2046:2050
  }else{
    # End-of-century
    time_ls <- 2096:2100
  }
  for(ssp in c("ssp245","ssp585")){
    # Adjust CMIP6 projection for all target variables, for this ssp, during this time, and for this model
    # Initialize a list to store output
    out <- list()
    for(i in 1:length(varname_ls )){
      varname <- varname_ls[i]
      # Initialize a raster stack for this variable for all five years
      var_stack <- stack()
      # Stack five years together
      for(time in time_ls){
        CMIP6_proj <- readRDS(paste0(CMIP_raw_path,model_name,"_",time,"_",ssp,".rds"))
        # Get the variable
        var_stack <- stack(var_stack,CMIP6_proj[[paste0(varname,"_",ssp)]])
      }
      var_time <- names(var_stack)
      # Adjust CMIP6 projection
      var_stack <- var_stack - hist_CMIP_gs_mean_dif[[paste0(varname,"_",ssp)]]
      names(var_stack) <- var_time
      out[[i]] <- var_stack
    }
    names(out) <- varname_ls
    # Apply filters
    out$daily_SM[out$daily_SM > 1] <- 1
    out$daily_SM[out$daily_SM < 0] <- 0
    out$daily_maxVPD[out$daily_maxVPD < 0] <- 0
    out$daily_SM[out$daily_T < 0] <- NA
    out$daily_maxVPD[out$daily_T < 0] <- NA
    # Only output daily SM and daily max VPD
    out <- list(out$daily_SM,out$daily_maxVPD)
    names(out) <- c("daily_SM","daily_maxVPD")
    # Output for this model, this ssp, and this time
    saveRDS(out,paste0(Output_path,model_name,"_",time_name,"_",ssp,".rds"))
    print(paste("Complete",time_name,ssp))
  }
}


