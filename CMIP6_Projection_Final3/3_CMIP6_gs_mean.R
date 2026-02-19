# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.6.19

# This code is to get average of 15 CMIP6 models for gs mean
# And get daily CONUS mean (for making pdf)

#########
# Global
#########
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
library(raster)
library(ncdf4)
library(sf)
# Input path of adjusted CMIP6
Input_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/CMIP6_adjusted/"
models_ls <- c("CMCC-ESM2","CanESM5","CanESM5-1","EC-Earth3","INM-CM4-8",
               "INM-CM5-0","IPSL-CM6A-LR","KACE-1-0-G","MIROC6","MPI-ESM1-2-HR",
               "MPI-ESM1-2-LR","MRI-ESM2-0","NorESM2-LM","NorESM2-MM","TaiESM1")
arrayid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#######
# Main
#######
ssp <- c("ssp245","ssp585")[as.numeric(substr(arrayid,1,1))]
time_name <- c("Mid","End")[as.numeric(substr(arrayid,2,2))]
varname <- c("daily_SM","daily_maxVPD")[as.numeric(substr(arrayid,3,3))]

# Stack 15 models for daily variables =============================================
# Initialize a stack to store 15 models
model_stack <- stack()
for(i in 1:length(models_ls)){
  model_name <- models_ls[i]
  CMIP_adj <- readRDS(paste0(Input_path,model_name,"_",time_name,"_",ssp,".rds"))
  # Get the target variable
  CMIP_adj <- CMIP_adj[[varname]]
  # Record layers and names
  layers <- nlayers(CMIP_adj)
  CMIP_time <- names(CMIP_adj)
  model_stack <- stack(model_stack,CMIP_adj)
  print(paste0("Complete",round(i/15,2)*100,"%"))
}
# Get daily mean of the 15 models
index <- rep(seq(1:layers),length(models_ls))
model_stack <- stackApply(model_stack,index,fun=mean,na.rm=T)
names(model_stack) <- CMIP_time

# Get daily CONUS mean ===========================================================
# Initialize a vector to store daily CONUS mean
daily_CONUS <- c()
for(layer in 1:nlayers(model_stack)){
  daily_CONUS_tmp <- mean(values(model_stack[[layer]]),na.rm=T)
  daily_CONUS <- c(daily_CONUS,daily_CONUS_tmp)
}
write.csv(daily_CONUS,paste0(Input_path,"Mean_models_daily_CONUS_mean_",time_name,"_",ssp,"_",varname,".csv"))

# Get gs mean =====================================================================
model_stack <- mean(model_stack,na.rm=T)
saveRDS(model_stack,paste0(Input_path,"Mean_models_gs_mean_",time_name,"_",ssp,"_",varname,".rds"))

print("All done!!")



