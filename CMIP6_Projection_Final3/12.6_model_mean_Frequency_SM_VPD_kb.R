# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.8.29

# This code is to get mean of frequency upper and lower boundaries
# From 15 models
# For the two k and b combinations
# For only SM or VPD change

# Note: k1 was calculated using mean slope + sd
# k2 was calculated using mean slope - sd

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
Freq_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_SM_VPD_kb/"

# List of 15 models
models_ls <- c("CMCC-ESM2","CanESM5","CanESM5-1","EC-Earth3","INM-CM4-8",
               "INM-CM5-0","IPSL-CM6A-LR","KACE-1-0-G","MIROC6","MPI-ESM1-2-HR",
               "MPI-ESM1-2-LR","MRI-ESM2-0","NorESM2-LM","NorESM2-MM","TaiESM1")
time_name <- "End"
ssp <- "ssp245"

########
# Main
########
# Initialize a stack to store mean frequency for the two kb combinations
# For the two changing variables (SM or VPD)
out <- stack()
for(kb in 1:2){
  for(var_chg in c("SM","VPD")){
    # Initialize a stack to store all rasters from 15 models, for the current kb combination, and the current changing variable
    freq_stack <- stack()
    for(model_name in models_ls){
      freq <- readRDS(paste0(Freq_path,model_name,"_",time_name,"_",ssp,"_",var_chg,"_chg_freq_kb",kb,".rds"))
      freq_stack <- stack(freq_stack,freq)
    }
    # Get mean of this stack
    freq <- stackApply(freq_stack,rep(1,nlayers(freq_stack)),fun=mean,na.rm=T)  
    out <- stack(out,freq)
    print(c(var_chg,kb))
  }
}
names(out) <- c("SM_kb1","VPD_kb1","SM_kb2","VPD_kb2")
saveRDS(out,paste0(Freq_path,"Model_mean_",time_name,"_",ssp,"_freq_SM_VPD_kb.rds"))
