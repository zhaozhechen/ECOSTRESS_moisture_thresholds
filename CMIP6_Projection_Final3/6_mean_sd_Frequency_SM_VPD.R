# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.8.7

# This code is to get mean and sd of frequency of projections that pass
# JOINT thresholds
# For two scenarios:
# (1) Only Change SM to CMIP6 SM in 2100 (ssp585)
# (2) Only Change VPD to CMIP6 VPD in 2100 (ssp585)
# While keeping the other variable at the original historical values

#########
# Global
#########

dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
library(raster)
library(ncdf4)
library(sf)
library(dplyr)

# Path to the frequency stacks
Freq_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_SM_VPD/"

#######
# Main
#######
time_name <- "End"
ssp <- "ssp585"
var_chg <- "SM_chg"
for(var_chg in c("SM_chg","VPD_chg")){
  # Get all file names for this scenarios (they are results from 15 models)
  file_names <- dir(Freq_path)[grepl(paste0(time_name,"_",ssp,"_",var_chg),dir(Freq_path))]
  # Loop over the files, and stack freq for 15 models
  # Initialize a raster stack to store all freq
  freq <- stack()
  for(i in 1:length(file_names)){
    freq <- stack(freq,readRDS(paste0(Freq_path,file_names[i])))
  }
  # Get the average and sd of 15 models
  freq_mean <- stackApply(freq,rep(1,nlayers(freq)),fun=mean,na.rm=T)
  freq_sd <- stackApply(freq,rep(1,nlayers(freq)),fun=sd,na.rm=T)
  out <- stack(freq_mean,freq_sd)
  names(out) <- c("Freq_mean","Freq_sd")
  saveRDS(out,paste0(Freq_path,"Model_mean_",time_name,"_",ssp,"_",var_chg,"_freq.rds"))  
}



