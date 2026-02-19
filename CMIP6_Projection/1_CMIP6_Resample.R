# Author: Zhaozhe Chen (chen.8926@ous.edu)
# Date: 2024.6.18

# This code is to get daily maxVPD, daily T, and daily SM 
# For the 15 models
# For the two ssp
# For each year

# The first digit in arrayid is for SSP
# 1: SSP245
# 2: SSP585
# The later two digits are for year index in the year list 

#########
# Global
#########
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
library(raster)
library(ncdf4)
library(sf)
library(stringr)

# Input path of CMIP6 daily
CMIP_path <- "/fs/ess/PAS2204/Data/CMIP6/Data/Daily/"
# Input path of combined 0.25D rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")
# This is the output for lists of raster stacks
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/CMIP6_raw/"
# List of 15 models
models_ls <- c("CMCC-ESM2","CanESM5","CanESM5-1","EC-Earth3","INM-CM4-8",
               "INM-CM5-0","IPSL-CM6A-LR","KACE-1-0-G","MIROC6","MPI-ESM1-2-HR",
               "MPI-ESM1-2-LR","MRI-ESM2-0","NorESM2-LM","NorESM2-MM","TaiESM1")
# Year to process
year_ls <- c(2018:2022,2046:2050,2096:2100)
# Decide which ssp and which year to process
arrayid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

############
# Functions
############
# Input includes 
# var_name: variable name
# model_name: model name
# ssp: ssp245 or ssp585
# year: the target year
Get_file <- function(var_name,model_name,ssp,year){
  # Use 05-01 of each year as a reference time
  time <- as.Date(paste0(year,"-05-01"))
  # Path of this var at this ssp
  nc_path <- paste0(CMIP_path,var_name,"_",ssp,"/")
  # File names of all these files for this model
  nc_files <- dir(nc_path)[grepl(model_name,dir(nc_path))]
  # Find the file name that contains the target date (time)
  for(j in 1:length(nc_files)){
    # Get the time of this file
    nc_time <- str_extract_all(nc_files,"\\d{8}")[[j]]
    nc_time <- as.Date(paste(substr(nc_time,1,4),substr(nc_time,5,6),substr(nc_time,7,8),sep="-"))
    # Check if this nc_time contains the target date
    if(time >= nc_time[1] & time <= nc_time[2]){
      nc_path <- paste0(nc_path,nc_files[j])
      break
    }
  }
  return(nc_path)  
}  

# This function is to get raster stacks of growing season for this year
# For this variable, at this ssp, for the target model
# And resample it to 0.25D
Stack_range <- function(var_name,model_name,ssp,year){
  nc_file <- Get_file(var_name,model_name,ssp,year)
  nc_stack <- stack(nc_file)
  # Get the time of these nc files
  nc_year <- substr(names(nc_stack),2,5)
  nc_month <- substr(names(nc_stack),7,8)
  nc_day <- substr(names(nc_stack),10,11)
  nc_time <- as.Date(paste(nc_year,nc_month,nc_day,sep="-"))
  # Only keep growing season of this year
  index <- which(nc_time >= as.Date(paste0(year,'-05-01')) &
                   nc_time <= as.Date(paste0(year,"-09-30")))
  nc_stack <- nc_stack[[index]]
  # Crop to CONUS and resample to the resolution of 0.25D
  extent(nc_stack)[1] <- extent(nc_stack)[1] - 360
  extent(nc_stack)[2] <- extent(nc_stack)[2] - 360
  nc_stack <- crop(nc_stack,raster_all)
  nc_stack <- resample(nc_stack,raster_all)
  return(nc_stack)
}

# This function is to calculate VPD
# Using daily maximum T to approximate daily maximum VPD
# Input is tasmax, unit:K
# hurs, unit: %
Cal_VPD <- function(temp,hurs){
  temp <- temp - 273.15
  # Relative humidity greater than 100 -> 100
  hurs[hurs>100] <- 100
  # Calculate es [kPa]
  es <- 0.6108*exp(17.27*temp/(temp+237.3))
  # Calculate VPD [kPa]
  VPD <- es*(1-hurs/100)
  VPD[VPD<0] <- 0
  names(VPD) <- names(temp)
  return(VPD)
}

#######
# Main
#######
# SSP for the current job
ssp_ls <- c("ssp245","ssp585")
ssp <- ssp_ls[as.numeric(substr(arrayid,1,1))]
# Year for the current job
year <- year_ls[as.numeric(substr(arrayid,2,3))]

# Loop over the 15 models
for(i in 1:length(models_ls)){
  # Initialize a list to store output rasters
  out <- list()
  model_name <- models_ls[i]
  # Daily maximum T [K]
  suppressWarnings(tasmax <- Stack_range("tasmax",model_name,ssp,year))  
  # Relative humidity [%]
  suppressWarnings(hurs <- Stack_range("hurs",model_name,ssp,year))
  # Calculate daily maximum VPD
  print("Calcualting daily maximum VPD ...")
  daily_maxVPD <- Cal_VPD(tasmax,hurs) 
  out[[1]] <- daily_maxVPD
  rm(tasmax,hurs)
  # Daily average T [K]
  suppressWarnings(tas <- Stack_range("tas",model_name,ssp,year))
  # Convert temperature unit to C
  tas <- tas - 273.15  
  out[[2]] <- tas
  # Daily average SM [kg m-2]
  suppressWarnings(SM <- Stack_range("mrsos",model_name,ssp,year))
  # Convert unit of SM to mm3/mm3
  SM <- SM/100  
  out[[3]] <- SM
  names(out) <- paste0(c("daily_maxVPD","daily_T","daily_SM"),"_",ssp)
  # Output this raster stack list
  saveRDS(out,paste0(Output_path,model_name,"_",year,"_",ssp,".rds"))
  print(paste("Complete models",round(i/15,2)*100,"%"))
}
print("All done !!!")






