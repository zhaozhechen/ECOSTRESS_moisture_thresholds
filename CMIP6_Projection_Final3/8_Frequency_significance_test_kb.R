# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.9.9

# This code is to test significant difference between frequency from historical climate
# vs frequency from 15 models (changing only SM or VPD, and changing both)
# For two kb combinations (alpha+sd, alpha-sd)
# For each pixel

#########
# Global
#########

dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
library(raster)
library(ncdf4)
library(sf)
library(dplyr)
library(ggplot2)
library(cowplot)

# For changing both SM and VPD, projected (of 15 models) and historical frequency path
freq_path_both <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_kb/"
# For changing only SM or VPD, projected frequency (of 15 models) path
# Note, for the same kb combination, it uses the same as in above historical frequency
freq_path_SM_VPD <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_SM_VPD_kb/"

# List of 15 models
models_ls <- c("CMCC-ESM2","CanESM5","CanESM5-1","EC-Earth3","INM-CM4-8",
               "INM-CM5-0","IPSL-CM6A-LR","KACE-1-0-G","MIROC6","MPI-ESM1-2-HR",
               "MPI-ESM1-2-LR","MRI-ESM2-0","NorESM2-LM","NorESM2-MM","TaiESM1")

##########
# Function
##########
# This function is to output p-value raster
# Input includes kb (1 or 2)
# ssp (ssp585 or ssp245)
# time_name ("Mid" or "End")
# var_chg ("" or "SM_chg_" or "VPD_chg_")
Cal_pvalue <- function(kb,ssp,time_name,var_chg){
  if(var_chg == ""){
    # Based on var_chg, determine the path of frequency
    freq_path <- freq_path_both  
  }else{
    freq_path <- freq_path_SM_VPD
  }
  
  # Historical frequency for the current kb combination, this is the same for
  # Either changing only one variable or changing both variables
  Hist_freq <- readRDS(paste0(freq_path_both,"Historical_freq_kb",kb,".rds"))
  # Stack frequency of the 15 models for this scenario
  Proj_freq <- stack()
  # Loop over each model
  for(i in 1:length(models_ls)){
    model_name <- models_ls[i]
    # Get projected frequency 
    Proj_freq_tmp <- readRDS(paste0(freq_path,model_name,"_",time_name,"_",ssp,"_",var_chg,"freq_kb",kb,".rds"))
    Proj_freq <- stack(Proj_freq,Proj_freq_tmp)
  }
  # Compare the historical frequency vs 15 CMIP frequencies for each pixel
  # Using Wilcox Signed-rank test
  # Initialize a vector to store all p-values
  p_all <- c()
  for(i in 1:ncell(Hist_freq)){
    # Get 15 models' frequency for this pixel
    Proj_freq_px <- values(Proj_freq)[i,]
    if(is.na(Hist_freq[i])){
      # If no historical frequency,p-value = NA
      p_tmp <- NA
    }else{
      # Do Wilcox test
      p_tmp <- wilcox.test(Proj_freq_px,mu= values(Hist_freq)[i])$p.value
    }
    p_all <- c(p_all,p_tmp)
  }
  # Construct a raster for p-value
  p_raster <- Hist_freq
  values(p_raster) <- p_all
  names(p_raster) <- paste0("p_value_",time_name,"_",ssp,"_",var_chg,"kb",kb)
  print(paste("Complete",time_name,ssp,var_chg,kb))  
  return(p_raster)
}

########
# Main
########
# Get p-value rasters ======================================================
# Initialize a stack to store p-value rasters for all scenarios
p_value_stack <- stack()
ssp <- "ssp585"
time_name <- "End"
for(kb in c(1,2)){
  for(var_chg in c("","SM_chg_","VPD_chg_")){
    p_value <- Cal_pvalue(kb,ssp,time_name,var_chg)
    p_value_stack <- stack(p_value_stack,p_value)
  }
}

# Output p-value rasters
saveRDS(p_value_stack,paste0(freq_path_SM_VPD,"p_value_rasters.rds"))



