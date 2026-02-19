# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.3.10

# This code is to get RF variables at 210m
# Add soil hydraulic traits, diversity variables, and aridity index
# Output full rasters (no theta or uncertainty mask)
# Reference code: /fs/ess/PAS2204/Code/CONUS_Threshold/CONUS_RF_Final_4/RF_0.25D_Data_Processing.R

##########
# Global
##########
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)

library(readxl)
library(raster)
library(ncdf4)
library(vegan)

# Determines which block to process
#Block_ID <- 1
Block_ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID')) 

# Input for Thresholds along with other already resampled variables
feature_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/RF_Samples/210m_all/"
# HWSD soil texture sand fraction
T_Sand <- raster("/fs/ess/PAS2204/Data/HWSD_1247/data/T_SAND.nc4")
# Data path to canopy height (https://glad.umd.edu/dataset/gedi)
CH <- raster("/fs/ess/PAS2204/SharedData/GLAD2019/Forest_height_2019_NAM.tif")
# This is the path to DEM data
DEM_path <- "/fs/ess/PAS2376/terrain_DEM/"
# NLCD Land cover
LC <-raster("/fs/ess/PAS2204/SharedData/nlcd_2019/nlcd_2019_land_cover_l48_20210604.img")
# Soil hydraulic traits path (use surface parameter)
SHT_path <- "/fs/ess/PAS2204/SharedData/Soil_Hydraulics/Hydraul_Param_SoilGrids_Schaap_sl1.nc"
# Aridity index path
AI_path <- "/fs/ess/PAS2204/SharedData/Aridity_Index/"

# Output path
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/210m_rasters_all/"

########
# Main
########

############################
# DEM
# List the extent of the DEM rasters
DEM_file <- dir(DEM_path)
# extents for these layers
DEM_xmin <- -as.numeric(substr(DEM_file,6,8))
DEM_xmax <- DEM_xmin + 1
DEM_ymin <- as.numeric(substr(DEM_file,2,3))
DEM_ymax <- DEM_ymin + 1
DEM_extent <- list()
for(i in 1:length(DEM_file)){
  D_extent <- extent(DEM_xmin[i],DEM_xmax[i],DEM_ymin[i],DEM_ymax[i])
  DEM_extent[i] <- D_extent
}

####################################
# Process original feature rasters 
raster_stack <- readRDS(paste(feature_path,"Var_stack_210m_Block_",Block_ID,".rds",sep=""))
# Convert P50s_P50X to P50S
# P50S = P50s_P50X * P50
raster_stack$P50S <- raster_stack$P50s_P50x * raster_stack$P50
# Calculate theta_sd
raster_stack$Theta_sd <- atan(raster_stack$Theta)/pi*180

############################
# Soil hydraulic traits
# Porosity # Unit m3m-3
Porosity <- raster(SHT_path,varname = "mean_theta_s_0cm")
Porosity <- crop(Porosity,raster_stack)
raster_stack$Porosity <- resample(Porosity,raster_stack)
# n_parameter # unit -
n_parameter <- raster(SHT_path,varname = "n_fit_0cm")
n_parameter <- crop(n_parameter,raster_stack)
raster_stack$n_parameter <- resample(n_parameter,raster_stack)
# Ks Average saturated hydraulic conductivity # unit cmd-1
Ks <- raster(SHT_path,varname = "mean_Ks_0cm")
Ks <- crop(Ks,raster_stack)
raster_stack$Ks <- resample(Ks,raster_stack)

################################
# Sand fraction
T_Sand <- crop(T_Sand,raster_stack)
raster_stack$T_Sand <- resample(T_Sand,raster_stack)

################################
# Aridity index
AI <- readRDS(paste0(AI_path,"AI_Block_",Block_ID,".rds"))
raster_stack$AI <- resample(AI,raster_stack)

#########################################################################
# Get diversity variables

# Initialize a raster to store CH_sd
raster_stack$CH_sd <- raster_stack$Theta
values(raster_stack$CH_sd) <- NA
# Initialize a raster to store LC_shannon
raster_stack$LC_shannon <- raster_stack$CH_sd
# Initialize a raster to store DEM_sd
raster_stack$DEM_sd <- raster_stack$CH_sd

# Loop over cells with valid Theta in this raster_stack and extract required variables
Index <- which(!is.na(values(raster_stack$Theta)))
# Randomly take 0.1% sample
set.seed(1)
Index <- sort(sample(Index,size=length(Index)*0.001))
# Only loop over 0.1% of valid Theta cells 
for(i in Index){
  # Extent of this cell
  c_extent <- extentFromCells(raster_stack,i)
  # Get Study_area, which is the area of this cell
  Study_area <- raster(c_extent,resolution=res(raster_stack),crs="+proj=longlat +datum=WGS84")
  
  ##################
  # Canopy height sd
  CH_tmp <- crop(CH,c_extent)
  raster_stack$CH_sd[i] <- sd(values(CH_tmp),na.rm=T)
  
  ##################
  # LC diversity
  LC_area <- Study_area
  suppressWarnings(LC_area <- projectRaster(LC_area,crs=crs(LC)))
  LC_tmp <- crop(LC,LC_area)
  # Get the number of individual elements in LC
  LC_tmp <- table(values(LC_tmp))
  # Calculate Shannon Index of LC
  raster_stack$LC_shannon[i] <- diversity(LC_tmp,index="shannon")
  
  ##################
  # Get DEM sd
  # Find DEM layer that contains this pixel
  DEM_index <- NA
  for(j in 1:length(DEM_file)){
    if(!is.null(raster::intersect(DEM_extent[[j]],Study_area))){
      DEM_index <- j
    }
  }
  # If there is no such DEM file, give NA
  if(is.na(DEM_index)){
    raster_stack$DEM_sd[i] <- NA
  }else{
    # Read the DEM raster
    DEM <- raster(paste(DEM_path,DEM_file[DEM_index],sep=""))
    # Crop to this pixel
    DEM <- crop(DEM,Study_area)
    # Get DEM sd
    raster_stack$DEM_sd[i] <- sd(values(DEM),na.rm=T)  
  }
  
  print(paste("Complete",round(i/max(Index),4)*100,"%"))
}

# Output this entire raster stack
saveRDS(raster_stack,paste0(Output_path,"RF_raster_all_Block_",Block_ID,".rds"))
print("All done")



