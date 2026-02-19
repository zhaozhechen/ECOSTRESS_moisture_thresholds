# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.3.10

# This code is to combine 0.25D and 210m rasters from all 50 blocks to 1 raster, respectively
# Note: takes too long to combine 210m. Did not use this code to combine 210m

#########
# Global
#########

dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)

library(raster)
library(sf)

# Manually change 0.25D or 210m
Resolution <- "0.25D"

# Input path of 0.25D or 210m rasters for 50 blocks
Input_path <- paste0("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/",Resolution,"_rasters_all/")

# Shapefile of CONUS 50 grids
CONUS <- st_read("/fs/ess/PAS2204/SharedData/CONUS_50_shp/CONUS_50_grids.shp")

# Output of combined rasters
Output_path <- paste0("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/",Resolution,"_raster_combined/")

#######
# Main
#######
# Initialize a raster stack to store all raster stacks
if(Resolution == "0.25D"){
  res <- c(0.25,0.25)
}else if(Resolution == "210m"){
  res <- c(0.00189,0.00189)
}
raster_all <- raster(extent(CONUS),res=res,crs="+proj=longlat +datum=WGS84")

for(Block_ID in 1:50){
  raster_stack <- readRDS(paste0(Input_path,"RF_raster_all_Block_",Block_ID,".rds"))
  raster_stack <- raster_stack[[var_ls]]
  # For LC, use ngb to resample
  LC <- resample(raster_stack$LC,raster_all,method="ngb")
  raster_stack <- resample(raster_stack,raster_all)
  raster_stack$LC <- LC
  raster_all <- merge(raster_all,raster_stack)
  print(paste("Complete",Block_ID/50*100,"%"))
}

names(raster_all) <- names(raster_stack)

# Output this combined raster
saveRDS(raster_all,paste0(Output_path,"Combined_raster.rds"))








