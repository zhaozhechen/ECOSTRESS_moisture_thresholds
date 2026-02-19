# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.6.20

# This code is to get mean and sd of frequency and delta_VPD
# From results of the 15 models
# For the two time periods, and two ssp

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
Freq_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_V2/"

#######
# Main
#######

#time_name <- "Mid"
#ssp <- "ssp585"
for(time_name in c("Mid","End")){
  for(ssp in c("ssp245","ssp585")){
    
    # Get all file names for this scenario (they are results from the 15 models)
    file_names <- dir(Freq_path)[grepl(paste0(time_name,"_",ssp,"_freq"),dir(Freq_path))]
    # Loop over the files, and stack freq and delta_VPD for 15 models
    # Initialize a raster to store all freq
    freq <- stack()
    # Initialize a raster to store all delta_VPD
    delta_VPD <- stack()
    for(i in 1:length(file_names)){
      freq <- stack(freq,readRDS(paste0(Freq_path,file_names[i]))$Freq)
      delta_VPD <- stack(delta_VPD,readRDS(paste0(Freq_path,file_names[i]))$delta_VPD)
    }
    # Get the average and sd of 15 models
    freq_mean <- stackApply(freq,rep(1,nlayers(freq)),fun=mean,na.rm=T)
    freq_sd <- stackApply(freq,rep(1,nlayers(freq)),fun=sd,na.rm=T)
    delta_VPD_mean <- stackApply(delta_VPD,rep(1,nlayers(freq)),fun=mean,na.rm=T)
    delta_VPD_sd <- stackApply(delta_VPD,rep(1,nlayers(freq)),fun=sd,na.rm=T)
    out <- stack(freq_mean,freq_sd,delta_VPD_mean,delta_VPD_sd)
    names(out) <- c("Freq_mean","Freq_sd","Delta_VPD_mean","Delta_VPD_sd")
    saveRDS(out,paste0(Freq_path,"Model_mean_freq_",time_name,"_",ssp,".rds"))    
  }
}







