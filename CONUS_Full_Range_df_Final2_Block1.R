# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2023.10.28

# This code is to get full-range df for CONUS data
# Output includes:
# 1. A list of Full-range df, including normalized daily SM, normalized daily VPD, and daily ESI (RData)
# Only overlapped SM, VPD, ESI points were kept to minimize storage space
# 2. Maximum and minimum of original SM, VPD, ESI, to recover their scales later (RData)
# 3. Coordinates, and LC of the pixels (csv) (Only required LC were kept)
# 4. ERA 5 climate data, including 
# mean(5-year growing season t2m monthly max), std(5-year growing season t2m monthly max)
# mean(5-year growing season t2m monthly min), std(5-year growing season t2m monthly min)
# mean(5-year growing season t2m monthly average),std(5-year growing season t2m monthly average)
# mean(5-year growing season netrad monthly average), std(5-year growing season netrad monthly average)
# mean(5-year growing season monthly accumulated precipitation), std(5-year growing season monthly accumulated precipitation)
# mean(5-year growing season monthly average daily max VPD), std(5-year growing season monthly average daily max VPD)
# mean(5-year growing season monthly average SM), std(5-year monthly average growing season SM)
# median(ESI) (From full-range df)

# Note: The output is in 210 m (3x3 pixels)

# Filters applied
# 1. Only keep 41,42,43,51,52,71,72,81,82 for LC
# 2. ETdaily filtered by ET QC_Flag (only 0,1,4,5 kept)
# 3. ETdaily = 0.01 removed
# 4. T < 0 filtered out
# 5. Non-growing season filtered out

# NO ETuncertainty used for filter

##########
# Global
##########

dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)

library(raster)
library(ncdf4)
library(sf)
library(scales)
library(lubridate)

# Determine which block to process (50 blocks in CONUS)
Block_ID  <- 1
# Determine which grid to process (100 grids in each block)
#Grid_ID   <-1
Grid_ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
# number of grids in one block
n_grid    <- 100

# Match IGBP to gs
#IGBP_ls <- c("DBF","ENF","MF","OSH","CSH","GRA","GRA","CRO","CRO")
NLCD_ls <- c(41,42,43,51,52,71,72,81,82)
GSlist <- 1/c(100,125,125,170,300,40,40,40,40)
# m/s, corresponding stomatal conductance from  MODIFIED_IGBP_MODIS_NOAH

# Dataset paths
# ERA5
ERA5_path <- paste("/fs/ess/PAS2204/SharedData/ERA5/","grid",Block_ID,"_",2018:2022,".nc",sep="")
# Shapfile of CONUS 50 grids
CONUS <- st_read("/fs/ess/PAS2204/SharedData/CONUS_50_shp/CONUS_50_grids.shp")
# NLCD path
LC_path    <-"/fs/ess/PAS2204/SharedData/nlcd_2019/nlcd_2019_land_cover_l48_20210604.img"
# ECOSTRESS ET path
ET_path <- paste("/fs/ess/PAS2204/SharedData/ECOSTRESS_ET/ET",Block_ID,sep="")
# ECOSTRESS ET QC path
ET_QC_path <- paste("/fs/ess/PAS2204/SharedData/ECOSTRESS_ET_uncertainty/ET_uncertainty_",Block_ID,sep="")
ET_QC_file <- dir(ET_QC_path)[grepl("QualityFlag_doy",dir(ET_QC_path))]
# SMAP path
SMAP_path <- paste("/fs/ess/PAS2204/Data/SMAP_L4_Root_Zone_Moisture/SMAP_L4_Grided/","SMAP_L4_Grid_",Block_ID,"/",sep="")

# Output path for full range df and meta data
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/Full_range_df/"
# Output path is named "./Block_ID/"
# Create Block path
Block_path <- paste(Output_path,"Block_",Block_ID,sep="")
if(!file.exists(Block_path)){
  dir.create(Block_path)  
}

#################
# Get Study Area
#################

# Get study area, which is one grid in one block
# Current Block
Block_area <- st_sf(CONUS[[2]][Block_ID])
# Divide the Block into 100 grids
Block_area <- st_make_grid(Block_area,n=c(sqrt(n_grid),sqrt(n_grid)))
# Get the extent of the current grid
Grid_area <- st_bbox(Block_area[Grid_ID])
Grid_extent <- extent(c(Grid_area[["xmin"]],Grid_area[["xmax"]],
                        Grid_area[["ymin"]],Grid_area[["ymax"]]))

############
# Functions
############

# This function stacks all ERA layers from 2018 to 2022,
# and crops the required variable
# Temporal resolution is hourly
# Input is variable name
# u10: 10 m u-component wind speed (m/s)
# v10: 10 m v-component wind speed (m/s)
# d2m: 2m dew point T (K)
# t2m: 2m T (K)
# zust: friction_velocity (m/s), equal to USTAR in AMF
# skt: skin temperature (K)
# ssr: surface_net_solar_radiation (J m**-2)
# sshf: surface_sensible_heat_flux (J m**-2)
# tp: total precipitation (m)

Read_ERA5 <- function(variable){
  print(paste("Reading ERA5",variable,"..."))
  # Read in all ERA5 layers for the variable
  ERA5 <- stack(ERA5_path[1:5],varname = variable)
  print(paste("Cropping ERA5",variable,"..."))
  # Crop it to the grid area
  # Use snap="out" to cover the whole grid
  ERA5 <- crop(ERA5,Grid_extent,snap="out")
  # Output the cropped ERA5 variable
  return(ERA5)
}

# Saturated water pressure, kPa
T2ES <- function(T){
  es <- 0.6108*exp(17.27*T/(T+237.3))
  return(es)
}

# This function calculates VPD
# Input is t2m and d2m
Cal_VPD <- function(t2m,d2m){
  ea <- T2ES(d2m-273.15)
  es <- T2ES(t2m-273.15)
  VPD <- es-ea
  return(VPD)
}

# Delta
PM_delta <- function(temp){
  tc <- temp - 273.15
  es <- 0.6108*exp(17.27*tc/(tc+237.3))*1000 # Unit Pa
  Delta <- 4098*es/(237.3+tc)^2 # Unit Pa/K
  return(Delta)
}

# Get growing season variables median and std
# input is daily variable raster and function
# fun is how to get monthly variable, e.g., monthly average, monthly max, or monthly min
gs_var <- function(daily_var,fun){
  # Stores output median and sd rasters
  out <- list()
  # Get growing season index
  gs_indices <- as.numeric(substr(ERA5_time,6,7))
  gs_indices[gs_indices <5 | gs_indices >9] <- NA
  # Get growing season monthly variable
  gs <- stackApply(daily_var,gs_indices,fun=fun,na.rm=T)
  # Get mean and std of monthly variable
  # This index includes NA and monthly 5:9
  indices <- rep(1,6)
  gs_mean <- stackApply(gs,indices,fun=mean,na.rm = T)
  gs_std <- stackApply(gs,indices,fun=sd,na.rm=T)
  out[[1]] <- gs_mean
  out[[2]] <- gs_std
  return(out)
}

# Reads ECOSTRESS ET related files
# varname is a string, should be "ETdaily_doy","QualityFlag_doy","ETdailyUncertainty_doy".
Read_file<-function(path,varname,i){
  path  <- path
  dir   <- dir(path)
  dir   <- dir[grepl(varname,dir)]
  file_path <- paste(path,dir[i],sep="/")
  raster <- raster(file_path)
  return(raster)
}

#########
# Main
#########

###############
# 1. Get LC
# Read in CONUS LC
LC <- raster(LC_path)

# Get an area a little larger than the Grid_extent to crop LC to cover the whole grid
LC_area <- Grid_extent + c(-1,1,-1,1)
LC_area <- raster(LC_area,resolution = c(0.25,0.25),
                  crs="+proj=longlat +datum=WGS84")
values(LC_area) <- NA
suppressWarnings(LC_area <- projectRaster(LC_area,crs = crs(LC)))

# Crop LC to LC_area first to reduce computational time
suppressWarnings(LC <- crop(LC,LC_area))
# Convert the projection of LC back
print("Reprojecting LC")
suppressWarnings(LC <- projectRaster(LC,crs="+proj=longlat +datum=WGS84",method="ngb"))
# Remove all LC that are not needed
LC[LC!=41 & LC!=42 & LC!=43 & LC!=51 & LC!=52 & 
     LC!=71 & LC!= 72 & LC!=81 & LC!=82] <- NA
# This is original 30m resolution LC in this grid

print("1. Complete LC pre-processing")
print("################################")

###################
# 2. Calculate PET 
# And get daily average PET,
# daily average T,
# daily maximum VPD

tryCatch({
  # Netrad,varname "ssr"
  netrad <- Read_ERA5("ssr") # Unit: J/m2
  # convert unit
  netrad <- netrad/3600 # Unit: W/m2
  # netrad < 0 equals 0
  netrad[netrad<0] <- 0
  
  # total precipitation, varname "tp"
  tp <- Read_ERA5("tp") # Unit: m
  
  # Surface T, varname "t2m"
  t2m <- Read_ERA5("t2m") # Unit: K
  
  # dewpoint T, varname "d2m"
  d2m <- Read_ERA5("d2m") # Unit:K
  
  # VPD unit kPa
  vpd <- Cal_VPD(t2m,d2m)
  # VPD < 0 equals 0
  vpd[vpd<0] <- NA
  
  # Delta unit Pa/K
  DEL <- PM_delta(t2m) 
  
  # Calculate PET unit mm/day
  rhoa  <- 1.225  # kg/m3
  Cp    <- 1005   # J/kg/K
  gamma <- 66     # Pa/K
  Lv    <- 2453e6 # J/m3
  
  # Wind speed
  # horizontal speed of air moving towards the east
  u10 <- Read_ERA5("u10") # unit m/s
  # horizontal speed of air moving towards the north
  v10 <- Read_ERA5("v10") # unit m/s
  # Get total horizontal ws
  u <- sqrt(u10^2 + v10^2)
  # friction_velocity unit m/s
  ustar <- Read_ERA5("zust")
  ga  <- 1/(u/ustar^2+6.2*ustar^(-2/3)) # m/s
  
  # Get PFT for the pixels
  # Aggregate LC to this ERA5 low resolution, using modal (the most frequent value)
  LC_ERA5 <- suppressWarnings(raster::aggregate(LC,fact=floor(res(netrad)/res(LC)),fun=modal,expand=TRUE))
  # Resample LC to this ERA5 resolution, using ngb
  LC_ERA5 <- resample(LC_ERA5,netrad,method="ngb")
  # Match NLCD to gs
  PFT <- match(LC_ERA5,NLCD_ls)
  # Get gs values
  values(PFT) <- GSlist[values(PFT)]
  gs <- PFT
  
  # Calulcate hourly PET
  PET <- (DEL*netrad+rhoa*Cp*ga*vpd*1000)/(DEL+gamma*(1+ga/gs))/Lv*1000*3600*24 # mm/day
  
  # Index to get daily variables
  daily_indices <- rep(1:(nlayers(netrad)/24),each = 24)
  # Get daily average PET
  PET_daily <- stackApply(PET,daily_indices,fun=mean,na.rm=T)
  # Get daily average netrad
  netrad_daily <- stackApply(netrad,daily_indices,fun=mean,na.rm = T)
  # Get daily average T in degree C
  T_daily <- stackApply(t2m,daily_indices,fun=mean,na.rm=T) - 273.15
  # Get daily maximum VPD
  VPD_daily_max <- stackApply(vpd,daily_indices,fun=max,na.rm=T)
  # Get daily total P
  p_daily <- stackApply(tp,daily_indices,fun=sum,na.rm=T)
  
  # Get ERA5 time
  ERA5_time <- seq(from=as.Date("2018-01-01"),to=as.Date("2022-12-31"),by=1)
  
  # Get median and std for climate data
  # Growing season monthly average T
  T_gs <- gs_var(T_daily,mean)
  # Growing season monthly max T
  T_max_gs <- gs_var(T_daily,max)
  # Growing season monthly min T
  T_min_gs <- gs_var(T_daily,min)
  # Growing season monthly accumulated P
  p_gs <- gs_var(p_daily,sum)
  # Growing season monthly average daily max VPD
  vpd_gs <- gs_var(VPD_daily_max,mean)
  # Growing season monthly average netrad
  netrad_gs <- gs_var(netrad_daily,mean)
  
  rm(PET,t2m,d2m,netrad,vpd,tp,v10,u10,u,ustar,ga,DEL,PFT,p_daily,netrad_daily)
  
  print("2. Complete daily PET, T, and VPD")
  print("################################")
  
  #######################
  # 3. Get 210m ETdaily
  
  # Use the resolution of the first ETdaily layer for the rest
  ET_area <- Read_file(ET_path,"ETdaily_doy",1)
  # Aggregate this to 210m footprint
  ET_area <- suppressWarnings(raster::aggregate(ET_area,fact=3,method=mean))
  ET_res <- res(ET_area)
  # Get study area for ET resample
  Study_area <- raster(Grid_extent,resolution = ET_res,
                       crs="+proj=longlat +datum=WGS84")
  
  # Initialize ET_name_all to store all original ET names
  ET_name_all <- c()
  # Initialize ET_daily_all to store all ETdaily
  ETdaily_all <- stack()
  
  for(i in 1:length(dir(ET_path))){
    # Read ETdaily rasters in turn
    ETdaily <- Read_file(ET_path,"ETdaily_doy",i)
    # If the current raster overlaps with Grid area, proceed
    if(suppressWarnings(!is.null(raster::intersect(extent(ETdaily),Grid_extent)))){
      tryCatch({
        # Crop ETdaily
        ETdaily <- crop(ETdaily,Grid_extent)
        # Record the name of the original ETdaily
        ET_name <- names(ETdaily)
        ET_name_all <- c(ET_name_all,ET_name)
        # layer name to find corresponding ET_QC and uncertainty file
        layer_name <- sub(".*doy","",ET_name)
        # Get corresponding ET_QC if exists
        if(sum(grepl(layer_name,ET_QC_file))==1){
          ET_QC <- ET_QC_file[grepl(layer_name,ET_QC_file)]
          ET_QC <- raster(paste(ET_QC_path,ET_QC,sep="/"))
          # Crop ET_QC
          ET_QC <- crop(ET_QC,Grid_extent)
          # QC filter, only keep QC=0,1,4,5
          values(ET_QC)[values(ET_QC)!=0 & 
                          values(ET_QC)!=1 & 
                          values(ET_QC)!=4 & 
                          values(ET_QC)!=5] <- NA
        }else{
          # If no corresponding QC file
          # Make a QC file of NA
          ET_QC <- ETdaily
          values(ET_QC) <- NA
        }
        # Mask ETdaily with good QC
        ETdaily <- mask(ETdaily,mask=ET_QC)
        # Remove ETdaily == 0.01
        ETdaily[abs(ETdaily - 0.01) < 10^-4] <- NA
        # Aggregate ETdaily to 210m footprint
        ETdaily <- suppressWarnings(raster::aggregate(ETdaily,fact=3,fun=mean,expand=TRUE))
        # Resample to the study area
        ETdaily <- resample(ETdaily,Study_area,method="bilinear")
        # Store all ETdaily layers
        ETdaily_all <- stack(ETdaily_all,ETdaily)
        },error=function(e){cat("ERROR: This grid and ET layer",i,conditionMessage(e),"\n")})
      }else{
        # If no overlap, jump to the next loop
        next()
      }
    print(paste("Complete ET",round(i/length(dir(ET_path))*100,2),"%"))
    }
  
  # Get time for ETdaily
  ET_name <- sub(".*doy","",ET_name_all)
  ET_year <- substr(ET_name,1,4)
  ET_day  <- as.numeric(substr(ET_name,5,7))-1
  ET_time <- as.Date(ET_day,
                     origin = as.Date(paste(ET_year,"-01-01",sep="")))
  # Aggregate same days
  ET_daily <- stackApply(ETdaily_all,indices = ET_time,fun=mean,na.rm=T)
  ET_time <- ET_time[!duplicated(ET_time)]
  
  rm(ETdaily_all,ET_name_all,ET_QC,ETdaily)
  #object_size(list=ls())
  print("3. Complete daily ET")
  print("################################")
  
  ####################
  # 4. Get SMAP daily
  
  # Initialize SMAP_name_all to sto all original names
  SMAP_name_all <- c()
  # Initialize SMAPdaily_all to store all SMAP layers
  SMAPdaily_all <- stack()
  
  for(i in 1:13){
    SMAP <- stack(paste(SMAP_path,i,"_raster_stack.tif",sep=""))
    # Crop to the Grid extent
    SMAP <- crop(SMAP,Study_area)
    # Store all SMAP rasters
    SMAPdaily_all <- stack(SMAPdaily_all,SMAP)
    # Read original SMAP name
    SMAP_name <- read.csv(paste(SMAP_path,i,"_names.csv",sep=""))$x
    SMAP_name_all <- c(SMAP_name_all,SMAP_name)
    print(paste("Complete SM",round(i/13*100,2),"%"))
  }
  
  # Get time for SMAP
  SM_name <- sub(".*gph_","",SMAP_name_all)
  SM_year <- substr(SM_name,1,4)
  SM_month <- substr(SM_name,5,6)
  SM_day <- substr(SM_name,7,8)
  SM_time <- as.Date(paste(SM_year,SM_month,SM_day,sep="-"))
  
  # Aggregate same days
  SM_daily <- stackApply(SMAPdaily_all,indices = SM_time,fun=mean,na.rm=T)
  SM_time <- SM_time[!duplicated(SM_time)]
  
  # Get mean and std for growing season monthly average SMAP
  # Indices for monthly average SMAP
  SM_month_indicies <- as.numeric(substr(SM_time,6,7))
  # Only keep growing season
  SM_month_indicies[SM_month_indicies<5 | SM_month_indicies>9] <- NA
  # Get monthly average SMAP
  SM_gs_monthly <- stackApply(SM_daily,SM_month_indicies,fun=mean,na.rm=T)
  # Get mean and std of growing season average monthly SMAP
  indices <- rep(1,nlayers(SM_gs_monthly))
  SM_gs_mean  <- stackApply(SM_gs_monthly,indices,fun=mean,na.rm=T)
  SM_gs_std   <- stackApply(SM_gs_monthly,indices,fun=sd,na.rm=T)
  
  rm(SMAPdaily_all,SMAP_name_all)
  print("4. Complete daily SM")
  print("################################")
  
  ##########################################
  # 5. Process raster values for data frame
  # Resample PET_daily to the same as ET_daily
  PET_daily <- resample(PET_daily,Study_area,method="ngb")
  print("Complete PET resampling")
  
  # Resample T_daily to the same as ET_daily
  T_daily <- resample(T_daily,Study_area,method="ngb")
  print("Complete T resampling")
  
  # Resample VPD_daily_max to the same as ET_daily
  VPD_daily_max <- resample(VPD_daily_max,Study_area,method="ngb")
  print("Complete VPD resampling")
  
  # Resample SM_daily to the same as ET_daily
  SM_daily <- resample(SM_daily,Study_area,method="ngb")
  print("Complete SM resampling")
  
  # Put their values in matrix, each column is one layer,each row is a pixel
  ET_daily <- values(ET_daily)
  PET_daily <- values(PET_daily)
  T_daily <- values(T_daily)
  VPD_daily_max <- values(VPD_daily_max)
  SM_daily <- values(SM_daily)
  
  # Resample climate data
  print("Resampling climate data")
  T_gs_mean <- resample(T_gs[[1]],Study_area,method="ngb")
  T_gs_std  <- resample(T_gs[[2]],Study_area,method="ngb")
  T_max_gs_mean <- resample(T_max_gs[[1]],Study_area,method="ngb")
  T_max_gs_std  <- resample(T_max_gs[[2]],Study_area,method="ngb")
  T_min_gs_mean <- resample(T_min_gs[[1]],Study_area,method="ngb")
  T_min_gs_std  <- resample(T_min_gs[[2]],Study_area,method="ngb")
  p_gs_mean <- resample(p_gs[[1]],Study_area,method="ngb")
  p_gs_std  <- resample(p_gs[[2]],Study_area,method="ngb")
  vpd_gs_mean <- resample(vpd_gs[[1]],Study_area,method="ngb")
  vpd_gs_std  <- resample(vpd_gs[[2]],Study_area,method="ngb")
  netrad_gs_mean  <- resample(netrad_gs[[1]],Study_area,method="ngb")
  netrad_gs_std   <- resample(netrad_gs[[2]],Study_area,method="ngb")
  SM_gs_mean  <- resample(SM_gs_mean,Study_area,method="ngb")
  SM_gs_std   <- resample(SM_gs_std,Study_area,method="ngb") 
  
  print("5. Complete raster values processing")
  print("################################")

  #########################
  # 6. Make full-range df
  
  print("################################")
  print("6. Start full-range df")
  
  # Output coordiates of all the pixels and their LC
  LC_ET <- suppressWarnings(raster::aggregate(LC,fact=floor(ET_res/res(LC)),fun=modal,expand=TRUE))
  # Resample LC to the ET resolution, using ngb
  LC_ET <- resample(LC_ET,Study_area,method="ngb")
  meta_df <- as.data.frame(LC_ET,xy=TRUE)
  names(meta_df) <- c("cor_x","cor_y","LC")
  meta_df$index <- rep(1:nrow(meta_df))
  # Put climate data together with LC and cor
  meta_df$T_gs_mean <- values(T_gs_mean)
  meta_df$T_gs_std  <- values(T_gs_std)
  meta_df$T_max_gs_mean <- values(T_max_gs_mean)
  meta_df$T_max_gs_std  <- values(T_max_gs_std)
  meta_df$T_min_gs_mean <- values(T_min_gs_mean)
  meta_df$T_min_gs_std  <- values(T_min_gs_std)
  meta_df$p_gs_mean <- values(p_gs_mean)
  meta_df$p_gs_std  <- values(p_gs_std)
  meta_df$vpd_gs_mean <- values(vpd_gs_mean)
  meta_df$vpd_gs_std  <- values(vpd_gs_std)
  meta_df$netrad_gs_mean  <- values(netrad_gs_mean)
  meta_df$netrad_gs_std   <- values(netrad_gs_std)
  meta_df$SM_gs_mean  <- values(SM_gs_mean)
  meta_df$SM_gs_std   <- values(SM_gs_std)
  
  # For LC = NA, no need to proceed
  # Record index of LC = NA
  index <- which(!is.na(meta_df$LC))

  # Initialize a list to store all full-range df
  Full_df_list <- list()
  # Initialize vecotrs to store scales
  SM_max <- rep(NA,nrow(meta_df))
  SM_min <- rep(NA,nrow(meta_df))
  VPD_max <- rep(NA,nrow(meta_df))
  VPD_min <- rep(NA,nrow(meta_df))
  ESI_max <- rep(NA,nrow(meta_df))
  ESI_min <- rep(NA,nrow(meta_df))
  ESI_median <- rep(NA,nrow(meta_df))
  # Loop over all pixels with LC!=NA, make df for each pixel
  for(i in index){
    # Get df for the ith pixel
    ET_df <- data.frame(time = ET_time,
                        ET_daily = ET_daily[i,])
    ERA5_df <- data.frame(time = ERA5_time,
                          PET_daily = PET_daily[i,],
                          T_daily = T_daily[i,],
                          VPD_daily_max = VPD_daily_max[i,])
    SM_df <- data.frame(time = SM_time,
                        SM_daily = SM_daily[i,])
    # Combine these df
    df <- merge(ERA5_df,ET_df,by="time",all=T)
    df <- merge(df,SM_df,by="time",all=T)
    
    # Calculate ESI
    df$ESI_daily <- df$ET_daily/df$PET_daily
    # ESI > 1 equals 1
    df$ESI_daily[df$ESI_daily>1] <- 1
    
    # Apply filters
    # Filter out T < 0
    filter_T <- which(df$T_daily < 0)
    # Filter out non-growing season
    month <- as.numeric(substr(df$time,6,7))
    filter_S <- which(month < 5 | month > 9)
    df[filter_T,] <- NA
    df[filter_S,] <- NA
    # Remove filtered rows
    df <- df[!is.na(df$time),]
    
    # Record the maximum and minimum of SM, VPD, ESI to scale back
    SM_max[i] <- max(df$SM_daily,na.rm=T)
    SM_min[i] <- min(df$SM_daily,na.rm=T)
    VPD_max[i] <- max(df$VPD_daily_max,na.rm=T)
    VPD_min[i] <- min(df$VPD_daily_max,na.rm=T)
    ESI_max[i] <- max(df$ESI_daily,na.rm=T)
    ESI_min[i] <- min(df$ESI_daily,na.rm=T)
    ESI_median[i] <- median(df$ESI_daily,na.rm=T)
    
    # Normalize SM and VPD
    VPD <- rescale(df$VPD_daily_max)*100
    SM <- rescale(df$SM_daily)*100
    
    # Only keep SM and VPD overlapped with ESI
    df <- data.frame(ESI = df$ESI_daily,
                     VPD = VPD,
                     SM = SM)
    df <- na.omit(df)
    # Store this df into list
    Full_df_list[[i]] <- df
    print(paste("Complete df",round(i/max(index),4)*100,"%"))
  }
  
  # Combine meta_df
  meta_df$SM_max <- SM_max
  meta_df$SM_min <- SM_min
  meta_df$VPD_max <- VPD_max
  meta_df$VPD_min <- VPD_min
  meta_df$ESI_max <- ESI_max
  meta_df$ESI_min <- ESI_min
  meta_df$ESI_median <- ESI_median
  
  # Only keep LC ! =0
  meta_df <- meta_df[index,]
  
  # Output the list of Full-range df
  saveRDS(Full_df_list,file=paste(Block_path,"/Full_df_ls_Grid_",Grid_ID,".rds",sep=""))
  # Output meta df
  write.csv(meta_df,paste(Block_path,"/Meta_df_Grid_",Grid_ID,".csv",sep=""))
  
  print(paste("Complete full-range df and meta data for Block",Block_ID,"Grid",Grid_ID))
  
  # This error function is for the whole code
},error=function(e){
    cat("ERROR: ERA5 and this grid",
        conditionMessage(e),"\n")})

print(paste("Finish Block",Block_ID,"Grid",Grid_ID))







