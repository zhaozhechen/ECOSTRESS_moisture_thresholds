# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2026.8.17

# This code is to calculate how many projected growing-season days fall
# outside the historical range of daily maximum VPD and daily SM
# For the 15 CMIP6 models
# For mid-century and end-century
# For SSP245 and SSP585

# Historical ranges are calculated for each 0.25D pixel using the same
# historical observation files used in the final CMIP6 projection workflow

#########
# Global
#########
library(raster)

# Path to the historical SM empirical distributions for each pixel
Hist_SM_range_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/00_Data/Processed/Future_range/Reference/Hist_SM_range_ls.rds"
# Path to daily adjusted CMIP6 projections
CMIP6_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/00_Data/Processed/Future_range/CMIP6_adjusted/"
# Path to the final 0.25D threshold raster stack
Raster_stack_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/00_Data/Processed/Future_range/Reference/Combined_0.25D_Intercepts.rds"
# Output path
Output_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/02_Results/Future_outside_historical_range/"
# Output path for summary tables
Table_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/02_Results/Tables/"

# List of 15 models
models_ls <- c("CMCC-ESM2","CanESM5","CanESM5-1","EC-Earth3","INM-CM4-8",
               "INM-CM5-0","IPSL-CM6A-LR","KACE-1-0-G","MIROC6","MPI-ESM1-2-HR",
               "MPI-ESM1-2-LR","MRI-ESM2-0","NorESM2-LM","NorESM2-MM","TaiESM1")

time_name_ls <- c("End")
ssp_ls <- c("ssp245","ssp585")

# Metrics to calculate
metric_names <- c("VPD_below","VPD_above","VPD_outside",
                  "SM_below","SM_above","SM_outside","Either_outside")

############
# Functions
############

# Check if all required input files are available
Check_inputs <- function(){
  if(!file.exists(Hist_SM_range_path)){
    stop("The historical SM range list was not found: ",Hist_SM_range_path)
  }
  if(!file.exists(Raster_stack_path)){
    stop("The reference raster stack was not found: ",Raster_stack_path)
  }
  missing_CMIP <- c()
  for(model_name in models_ls){
    for(time_name in time_name_ls){
      for(ssp in ssp_ls){
        file_name <- paste0(CMIP6_path,model_name,"_",time_name,"_",ssp,".rds")
        if(!file.exists(file_name)){
          missing_CMIP <- c(missing_CMIP,file_name)
        }
      }
    }
  }
  if(length(missing_CMIP)>0){
    stop("Missing ",length(missing_CMIP)," adjusted CMIP6 files. First missing file: ",missing_CMIP[1])
  }
}

# Calculate historical minimum and maximum for each pixel
Get_historical_range <- function(raster_stack,analysis_mask){
  # Historical VPD minimum and maximum were saved with the final threshold rasters
  Hist_VPD_min <- raster_stack$VPD_min
  Hist_VPD_max <- raster_stack$VPD_max
  
  # Get historical SM minimum and maximum from the saved empirical distributions
  Hist_SM_range_ls <- readRDS(Hist_SM_range_path)
  if(length(Hist_SM_range_ls)!=ncell(analysis_mask)){
    stop("Hist_SM_range_ls and the threshold raster do not have the same number of cells")
  }
  SM_min_values <- rep(NA,ncell(analysis_mask))
  SM_max_values <- rep(NA,ncell(analysis_mask))
  for(i in 1:length(Hist_SM_range_ls)){
    if(is.function(Hist_SM_range_ls[[i]])){
      SM_values <- knots(Hist_SM_range_ls[[i]])
      SM_min_values[i] <- min(SM_values,na.rm=TRUE)
      SM_max_values[i] <- max(SM_values,na.rm=TRUE)
    }
    if(i%%1000==0){
      print(paste("Complete historical SM range",i,"out of",length(Hist_SM_range_ls)))
    }
  }
  Hist_SM_min <- raster(analysis_mask)
  Hist_SM_max <- raster(analysis_mask)
  values(Hist_SM_min) <- SM_min_values
  values(Hist_SM_max) <- SM_max_values
  
  Historical_range <- stack(Hist_VPD_min,Hist_VPD_max,
                            Hist_SM_min,Hist_SM_max)
  Historical_range <- mask(Historical_range,analysis_mask)
  names(Historical_range) <- c("VPD_min","VPD_max","SM_min","SM_max")
  return(Historical_range)
}

# Calculate number and fraction of days meeting a condition
Count_days <- function(condition,valid,analysis_mask){
  # Exclude days without valid input data from both numerator and denominator
  condition <- mask(condition,valid,maskvalue=0)
  valid_days <- calc(valid,fun=sum,na.rm=TRUE)
  outside_days <- calc(condition,fun=sum,na.rm=TRUE)
  outside_days[valid_days==0] <- NA
  outside_fraction <- outside_days/valid_days
  outside_days <- mask(outside_days,analysis_mask)
  outside_fraction <- mask(outside_fraction,analysis_mask)
  return(list(days=outside_days,fraction=outside_fraction))
}

# Calculate all out-of-range metrics for one model, time, and scenario
Cal_outside_range <- function(CMIP_file,Historical_range,analysis_mask){
  data <- readRDS(CMIP_file)
  VPD <- data$daily_maxVPD
  SM <- data$daily_SM
  
  VPD_valid <- !is.na(VPD)
  SM_valid <- !is.na(SM)
  Joint_valid <- VPD_valid & SM_valid
  
  # Get the number of days that VPD or SM fall beyond historical range across the 5 years
  VPD_below <- VPD < Historical_range$VPD_min
  VPD_above <- VPD > Historical_range$VPD_max
  VPD_outside <- VPD_below | VPD_above
  SM_below <- SM < Historical_range$SM_min
  SM_above <- SM > Historical_range$SM_max
  SM_outside <- SM_below | SM_above
  Either_outside <- VPD_outside | SM_outside
  
  condition_ls <- list(VPD_below,VPD_above,VPD_outside,
                       SM_below,SM_above,SM_outside,Either_outside)
  valid_ls <- list(VPD_valid,VPD_valid,VPD_valid,
                   SM_valid,SM_valid,SM_valid,Joint_valid)
  
  day_stack <- stack()
  fraction_stack <- stack()
  for(i in 1:length(condition_ls)){
    out <- Count_days(condition_ls[[i]],valid_ls[[i]],analysis_mask)
    day_stack <- stack(day_stack,out$days)
    fraction_stack <- stack(fraction_stack,out$fraction)
  }
  names(day_stack) <- metric_names
  names(fraction_stack) <- metric_names
  return(list(days=day_stack,fraction=fraction_stack))
}

# Calculate area-weighted mean for one raster
Weighted_mean <- function(r,area_raster){
  return(weighted.mean(values(r),values(area_raster),na.rm=TRUE))
}

#######
# Main
#######
dir.create(Output_path,recursive=TRUE,showWarnings=FALSE)
dir.create(paste0(Output_path,"Model_results/"),recursive=TRUE,showWarnings=FALSE)
dir.create(Table_path,recursive=TRUE,showWarnings=FALSE)

# Check required inputs
Check_inputs()

# Read the final threshold raster stack and use slope as the analysis mask
raster_stack <- readRDS(Raster_stack_path)
analysis_mask <- raster_stack$Slope
area_raster <- area(analysis_mask)

# Get historical ranges from the final saved threshold products
Historical_range <- Get_historical_range(raster_stack,analysis_mask)
saveRDS(Historical_range,paste0(Output_path,"Historical_VPD_SM_ranges.rds"))

# Initialize a list to store model statistics
statistics_ls <- list()
statistics_id <- 1

for(time_name in time_name_ls){
  for(ssp in ssp_ls){
    # Store out-of-range fractions for all 15 models
    model_fraction_all <- stack()
    
    for(model_name in models_ls){
      CMIP_file <- paste0(CMIP6_path,model_name,"_",time_name,"_",ssp,".rds")
      out <- Cal_outside_range(CMIP_file,Historical_range,analysis_mask)
      
      # Save the compact result for this model
      saveRDS(out,paste0(Output_path,"Model_results/",model_name,"_",
                         time_name,"_",ssp,"_outside_range.rds"))
      model_fraction_all <- stack(model_fraction_all,out$fraction)
      
      # Calculate CONUS statistics for each metric
      for(metric in metric_names){
        day_raster <- out$days[[metric]]
        fraction_raster <- out$fraction[[metric]]
        statistics_ls[[statistics_id]] <- data.frame(
          Model = model_name,
          Time = time_name,
          SSP = ssp,
          Metric = metric,
          Mean_days = Weighted_mean(day_raster,area_raster),
          Median_days = median(values(day_raster),na.rm=TRUE),
          Mean_fraction = Weighted_mean(fraction_raster,area_raster),
          Median_fraction = median(values(fraction_raster),na.rm=TRUE)
        )
        statistics_id <- statistics_id+1
      }
      print(paste("Complete",model_name,time_name,ssp))
    }
    
    # Calculate multi-model mean and standard deviation for each metric
    model_index <- rep(1:length(metric_names),times=length(models_ls))
    Model_mean <- stackApply(model_fraction_all,model_index,fun=mean,na.rm=TRUE)
    Model_sd <- stackApply(model_fraction_all,model_index,fun=sd,na.rm=TRUE)
    names(Model_mean) <- metric_names
    names(Model_sd) <- metric_names
    saveRDS(list(Mean=Model_mean,SD=Model_sd),
            paste0(Output_path,"Model_mean_",time_name,"_",ssp,"_outside_range.rds"))
  }
}

# Output statistics for all models
statistics_all <- do.call(rbind,statistics_ls)
write.csv(statistics_all,paste0(Table_path,"Future_outside_historical_range_by_model.csv"),
          row.names=FALSE)

# Summarize statistics across 15 models
statistics_mean <- aggregate(cbind(Mean_days,Median_days,
                                   Mean_fraction,Median_fraction)~Time+SSP+Metric,
                             data=statistics_all,FUN=mean,na.rm=TRUE)
statistics_sd <- aggregate(Mean_fraction~Time+SSP+Metric,
                           data=statistics_all,FUN=sd,na.rm=TRUE)
names(statistics_sd)[names(statistics_sd)=="Mean_fraction"] <- "SD_mean_fraction"
statistics_summary <- merge(statistics_mean,statistics_sd,
                            by=c("Time","SSP","Metric"))
write.csv(statistics_summary,paste0(Table_path,"Future_outside_historical_range_summary.csv"),
          row.names=FALSE)

print("All done !!!")
