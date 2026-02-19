# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.2.27

# This code is to process variables for RF at 0.25 degree

########
# Global
########

dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)

library(raster)
library(sf)
library(tidyr)
library(ncdf4)
library(stringr)
library(vegan)

# Determines which block to process
#Block_ID <- 6
Block_ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID')) 

# Data path to Liu's MDF hydraulic traits
PHT_path <- "/fs/ess/PAS2204/Data/Liu_MDF_Hydraulic_traits/"
# HWSD soil texture sand fraction
T_Sand <- raster("/fs/ess/PAS2204/Data/HWSD_1247/data/T_SAND.nc4")
# Data path to canopy height (https://glad.umd.edu/dataset/gedi)
CH <- raster("/fs/ess/PAS2204/SharedData/GLAD2019/Forest_height_2019_NAM.tif")
# Data path to rooting depth (https://wci.earth2observe.eu/thredds/catalog/usc/root-depth/catalog.html)
RD <- raster("/fs/ess/PAS2204/Data/Rooting_Depth/maxroot_North_America_CF.nc")
# This is the path to DEM data
DEM_path <- "/fs/ess/PAS2376/terrain_DEM/"
# NLCD Land cover
LC_30m <-raster("/fs/ess/PAS2204/SharedData/nlcd_2019/nlcd_2019_land_cover_l48_20210604.img")
# Processed 0.25 degree LC
LC <- readRDS(paste0("/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Raw_var_rasters/LC/LC_Block_",Block_ID,".rds"))
# Soil hydraulic traits path (use surface parameter)
SHT_path <- "/fs/ess/PAS2204/SharedData/Soil_Hydraulics/Hydraul_Param_SoilGrids_Schaap_sl1.nc"
# Shapefile of CONUS 50 grids
CONUS <- st_read("/fs/ess/PAS2204/SharedData/CONUS_50_shp/CONUS_50_grids.shp")

Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Raw_var_rasters/var_all/"

############
# Functions
############

# This function is to crop and resample the variable
var_process <- function(raster){
  raster <- crop(raster,Block_area)
  raster <- resample(raster,Block_area,method="ngb")
  return(raster)
}

# This function is to calculate sd for high resolution variables in each 0.25D pixel
# Loop over each cell and extract sd of raster
# Input is the raster at the original scale
raster_sd <- function(raster){
  sd <- Block_area
  for(i in 1:ncell(Block_area)){
    # Extent of this cell
    c_extent <- extentFromCells(Block_area,i)
    # raster in this cell at the original resolution
    c_raster <- crop(raster,c_extent)
    # Calculate sd
    sd[i] <- sd(values(c_raster),na.rm=T)
  }
  return(sd)
}

########
# Main
########

# Get Block_area
Block_area <- st_sf(CONUS[[2]][Block_ID])
Block_area <- raster(extent(Block_area),resolution = c(0.25,0.25),crs="+proj=longlat +datum=WGS84")

# Soil hydraulic traits (SHT) ==============================================
# n_paramter # unit -
n_parameter <- raster(SHT_path,varname = "n_fit_0cm")
n_parameter <- var_process(n_parameter)
# Ks Average saturated hydraulic conductivity # unit cmd-1
Ks <- raster(SHT_path,varname = "mean_Ks_0cm")
Ks <- var_process(Ks)
print("Complete Soil hydraulic traits")

# Sand fraction ============================================================
T_Sand <- var_process(T_Sand)
print("Complete Sand fraction")

# Plant hydraulic traits (PHT) =============================================
PHT_ls <- c("g1","C","gpmax","P50","P50s_P50x")
PHT_stack <- stack()
for(i in 1:length(PHT_ls)){
  PHT_tmp <- raster(paste0(PHT_path,"MDF_",PHT_ls[i],".nc"),
                    varname = paste0(PHT_ls[i],"_50"))
  PHT_tmp <- var_process(PHT_tmp)
  names(PHT_tmp) <- PHT_ls[i]
  PHT_stack <- stack(PHT_stack,PHT_tmp)
}
print("Complete Plant hydraulic traits")

# Rooting depth ============================================================
RD <- var_process(RD)
print("Complete Rooting depth")

# Canopy height ============================================================
CH_tmp <- crop(CH,Block_area)
# Pixel values 0-60 Forest canopy height in meters
# https://glad.umd.edu/dataset/gedi
values(CH_tmp)[values(CH_tmp)>=60] <- NA
# Mean CH after resampling to 0.25D
CH_mean <- var_process(CH_tmp)
# CH_sd in each pixel
CH_sd <- raster_sd(CH)
print("Complete Canopy Height nad CH sd")

# LC ==============================================
# LC at 0.25D (modal) is from previous first step in 1_0.25D_Full_df.R
# Initiate a LC_shannon raster
LC_shannon <- Block_area
for(i in 1:ncell(Block_area)){
  # Extent of this cell
  c_extent <- extentFromCells(Block_area,i)
  # Get Study_area, which is this cell
  c_area <- raster(c_extent,resolution=res(Block_area),crs=crs(Block_area))
  suppressWarnings(c_area <- projectRaster(c_area,crs=crs(LC_30m)))
  if(tryCatch(!is.null(raster::intersect(LC_30m,c_area)),error=function(e) return(FALSE))){
    LC_tmp <- crop(LC_30m,c_area)
    # Get the number of individual elements in LC_tmp
    LC_tmp <- table(values(LC_tmp))
    LC_shannon[i] <- diversity(LC_tmp,index="shannon")
  }else{
    LC_shannon[i] <- NA
  }
}
print("Complete LC and LC shannon")

# DEM ==============================================
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

# Initiate a DEM raster and DEM_sd raster
DEM <- Block_area
DEM_sd <- Block_area

for(i in 1:ncell(Block_area)){
  # Extent of this cell
  c_extent <- extentFromCells(Block_area,i)
  # Find DEM layer that contains this pixel
  DEM_index <- NA
  for(j in 1:length(DEM_file)){
    if(!is.null(raster::intersect(DEM_extent[[j]],c_extent))){
      DEM_index <- j
    }
  }
  # If there is no such DEM file, give NA to DEM and DEM_sd
  if(is.na(DEM_index)){
    DEM[i] <- NA
    DEM_sd[i] <- NA
  }else{
    # Read the DEM raster
    DEM_tmp <- raster(paste0(DEM_path,DEM_file[DEM_index]))
    # Crop to this pixel
    DEM_tmp <- crop(DEM_tmp,c_extent)
    DEM[i] <- mean(values(DEM_tmp),na.rm=T)
    DEM_sd[i] <- sd(values(DEM_tmp),na.rm=T)
  }
}
print("Complete DEM and DEM_sd")

# Output all rasters ===============================================================
# Stack these rasters together
raster_all <- stack(n_parameter,Ks,T_Sand,RD,PHT_stack,LC,LC_shannon,CH_mean,CH_sd,DEM,DEM_sd)
names(raster_all) <- c("n_fit","Ks","T_Sand","RD",names(PHT_stack),"LC","LC_shannon","CH","CH_sd","DEM","DEM_sd")

saveRDS(raster_all,paste0(Output_path,"var_all_Block_",Block_ID,".rds"))
print("All done")


