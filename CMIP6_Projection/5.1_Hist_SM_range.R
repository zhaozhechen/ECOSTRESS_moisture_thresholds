# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.7.9

# This code is to calculate historical SM range for each pixel
# Store each SM range cdf for each pixel as one element in the output list

#########
# Global
#########

dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
library(raster)
library(ncdf4)
library(sf)
library(dplyr)

# Input path of daily historical observations
hist_obs_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final2/Hist_obs/"
# This is the output to store the list of historical SM range
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Hist_obs/"

#######
# Main
#######

# Get full historical SM data ================================
# Files for historical observations
hist_obs_files <- paste0(hist_obs_path,dir(hist_obs_path))
# Historical obs does not have 2018.5 and 2018.6 SM, so skip
hist_obs_files <- hist_obs_files[3:25]
# Initialize a stack to store all historical SM
hist_SM <- stack()
for(i in 1:length(hist_obs_files)){
  hist_SM_tmp <- readRDS(hist_obs_files[i])$daily_SM
  # Stack all these historical SM rasters
  hist_SM <- stack(hist_SM,hist_SM_tmp)
}
print("Get full range SM data")

# Get the historical SM range for each pixel =============================
hist_SM_range_ls <- list()
for(i in 1:ncell(hist_SM)){
  # Check if all SM for this pixel is NA
  if(all(is.na(values(hist_SM)[i,]))){
    hist_SM_range_ls[[i]] <- NA
  }else{
    hist_SM_range_ls[[i]] <- ecdf(values(hist_SM)[i,]) 
  }
  
  if(i%%100 == 0){
    print(i)
  }
}

# Output this list
saveRDS(hist_SM_range_ls,paste0(Output_path,"Hist_SM_range_ls.rds"))

