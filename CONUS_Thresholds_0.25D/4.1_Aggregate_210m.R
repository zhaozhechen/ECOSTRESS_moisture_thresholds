# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.3.4

# This code is to aggregated 210m thresholds to 0.25D and get sd

#########
# Global
#########

dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)

library(raster)
library(sf)

# Input path of raw 0.25 degree threshold rasters
raster_0.25D_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Feature_rasters/"
# Input path of raw 210m threshold rasters
raster_210m_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/RF_Samples/210m_all/"
# Input path of 0.25 degree RF variables
var_0.25D_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Raw_var_rasters/var_all/"

# Shapefile of CONUS 50 grids
CONUS <- st_read("/fs/ess/PAS2204/SharedData/CONUS_50_shp/CONUS_50_grids.shp")

# Output of combined thresholds
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Raw_var_rasters/var_all_combined/"

########
# Main
########
# Initialize a raster to store 0.25D theta in all blocks
raster_0.25D <- raster(extent(CONUS),res=c(0.25,0.25),crs="+proj=longlat +datum=WGS84")
# Initialize a raster to store 0.25D RF variables in all blocks
var_0.25D <- raster_0.25D
# Initialize a raster to store aggregated theta from 210m
raster_Theta_agg <- raster_0.25D
# Initialize a raster to store aggregated theta sd from 210m
raster_Theta_agg_sd <- raster_0.25D

for(i in 1:length(dir(raster_0.25D_path))){
  # Get Block_ID 
  Block_ID <- as.numeric(sub(".*Block_*","",sub("*.rds","",dir(raster_0.25D_path)[i])))
  # Read in 0.25degree thresholds
  raster_0.25D_tmp <- readRDS(paste0(raster_0.25D_path,"features_Block_",Block_ID,".rds"))
  # Merge theta
  raster_0.25D <- merge(raster_0.25D,raster_0.25D_tmp)
  
  # Read in 0.25D RF variables
  var_0.25D_tmp <- readRDS(paste0(var_0.25D_path,"var_all_Block_",Block_ID,".rds"))
  # Merge them
  var_0.25D <- merge(var_0.25D,var_0.25D_tmp)
  
  # Read in 210m thresholds
  raster_210m <- readRDS(paste0(raster_210m_path,"Var_stack_210m_Block_",Block_ID,".rds"))
  # Only keep theta
  raster_210m <- raster_210m[["Theta"]]
  
  # Aggregate 210m theta to 0.25D
  Theta_agg <- raster::aggregate(raster_210m,fact = res(raster_0.25D)/res(raster_210m),fun=mean,na.rm=TRUE)
  # Resample to the same resolution as 0.25D raster
  Theta_agg <- resample(Theta_agg,raster_0.25D,method="ngb")
  # Merge theta
  raster_Theta_agg <- merge(raster_Theta_agg,Theta_agg)
  
  # Aggregate 210m theta to 0.25D, get sd
  Theta_agg_sd <- raster::aggregate(raster_210m,fact = res(raster_0.25D)/res(raster_210m),fun=sd,na.rm=TRUE)
  # Resample to the same resolution as 0.25D raster
  Theta_agg_sd <- resample(Theta_agg_sd,raster_0.25D,method="ngb")
  # Merge theta sd
  raster_Theta_agg_sd <- merge(raster_Theta_agg_sd,Theta_agg_sd)
  
  print(paste("Complete",round(i/length(dir(raster_0.25D_path)),2)*100,"%"))
}

# Stack these raster together
names(raster_0.25D) <- names(raster_0.25D_tmp)
names(var_0.25D) <- names(var_0.25D_tmp)
names(raster_Theta_agg) <- "Theta_agg210m"
names(raster_Theta_agg_sd) <- "Theta_agg210m_sd"

raster_all <- stack(raster_Theta_agg,raster_Theta_agg_sd,raster_0.25D,var_0.25D)

# Output this raster stack
saveRDS(raster_all,paste0(Output_path,"raster_all.rds"))





