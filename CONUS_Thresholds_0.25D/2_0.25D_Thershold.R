# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.2.27

# This code is to calculate thresholds at 0.25 degree
# All input raw data are at 0.25 degree

# Filters applied to get slopes and theta:
# 1. At least n(ESI!=1) >= 10 for the pixel
# 2. At least n_WLS number of points in each slice to fit a slice WLS
# 3. At least Slice_R2 of slice WLS to fit 0 point
# 4. Significantly wrong slice slopes are removed
# 5. Overextrapolated Fit 0 are removed, cannot exceed SM and VPD range in their normalized full-range df
# 6. At least n_Fit0 number of Fit0 points to fit a contour line
# 7. At most Slope_CV_th for the uncertainty of slopes (Slope_sd/Slope)
# 8. At least Slope_R2_th for the final contour slope
# 9. Theta < theta_th degree removed

# Filters applied to get thresholds
# 1. SM and VPD threshold quantiles must be in [0,100], which means they cannot exceed the original full-range scale
# 2. VPD threshold absolute values when fixing SM level must be in [0,6]
# 3. SM threshold absolute values when fixing VPD level must be in [0,1]
# These three filters applied together

########
# Global
########

dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)

library(raster)
library(dplyr)
library(scales)
library(sf)

Block_ID <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# This path is the output folder of full-range df for all blocks
Input_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Full_range_df/"
# Shapefile of CONUS 50 grids
CONUS <- st_read("/fs/ess/PAS2204/SharedData/CONUS_50_shp/CONUS_50_grids.shp")
# This path is the output folder for contour features
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Feature_rasters/"

# GS parameters
n_grid  <- 7
n_k     <- 3
sigma   <- 1

# bootstrapping iteration for slope sd
num_iteration <- 100

# WLS and contour line parameters
# At least n_WLS number of points in each slice to fit a slice WLS
n_WLS <- 4
# Slice R2 at least greater or equal to Slice_R2 to Fit0
Slice_R2 <- 0.8
# Least number of Fit0 point to fit a contour slope
n_Fit0 <- 3

# Slopes filter, Slope_sd/Slope at most Slope_CV_th
Slope_CV_th <- 100
# Slopes filter, R2 at least Slope_R2_th
Slope_R2_th <- 0.5
# Slope angle filter, angle theta < theta_th removed
theta_th <- -20

###########
# Functions
###########

# Gaussian filtration
gaussian.kernel <- function(sigma=2, s=5) {
  m <- matrix(ncol=s, nrow=s)
  mcol <- rep(1:s, s)
  mrow <- rep(1:s, each=s)
  x <- mcol - ceiling(s/2)
  y <- mrow - ceiling(s/2)
  m[cbind(mrow, mcol)] <- 1/(2*pi*sigma^2) * exp(-(x^2+y^2)/(2*sigma^2))
  m / sum(m)
}

gaussian.smooth <- function(x, sigma = sigma, n=n,type = mean, ...) {  
  if (!inherits(x, "RasterLayer")) stop("MUST BE RasterLayer OBJECT")
  gm <- gaussian.kernel(sigma=sigma, s=n)
  return(raster::focal(x, w = gm, fun = type, na.rm=TRUE, pad=FALSE, ...) )
} 

GSfiltration <- function(df){
  data_matrix   <- cbind(df$SM,df$VPD,df$Z)
  extent        <- extent(data_matrix[,1:2])
  data_raster   <- raster(extent,ncol=n_grid,nrow=n_grid)
  data_raster   <- raster::rasterize(data_matrix[,1:2],data_raster,data_matrix[,3],fun=mean)
  raster_GS     <- gaussian.smooth(data_raster,sigma=sigma,n=n_k)
  return(raster_GS)
}

# This function is to get normalized df
# Normalize ESI to 0-100
# Z = ESI - median(ESI)
Normalize_df <- function(df){
  ESI <- rescale(df$ESI)*100
  Z <- ESI-median(ESI)
  # Get df for GS filtration
  df <- data.frame(SM = df$SM,
                   VPD = df$VPD,
                   Z = Z)
  return(df)
}

# Weighted least square regression (WLS)
# Get weight matrix, which is the smoothed number of points
# Input is original full-range df
WLS_w_mat <- function(df){
  # Z values are 1 for all points
  data_matrix <- cbind(df$SM,df$VPD,rep(1,nrow(df)))
  extent        <- extent(data_matrix[,1:2])
  data_raster   <- raster(extent,ncol=n_grid,nrow=n_grid)
  # Sum of numbers of data points in each pixel
  data_raster   <- raster::rasterize(data_matrix[,1:2],data_raster,data_matrix[,3],fun=sum)
  data_raster   <- gaussian.smooth(data_raster,sigma=sigma,n=n_k,type = "sum")
  w_mat <- raster::as.matrix(data_raster)
  return(w_mat)
}

# Fit WLS for each slice and output slope, intercept, R2, and p-value
WLS_fit <- function(slice,coor,w){
  # Only proceed if there are at least n_WLS in the slice
  if(sum(!is.na(slice))>=n_WLS){
    WLS_df <- data.frame(coor,slice)
    # WLS fit
    lm <- lm(data=WLS_df,slice~coor,weights = w)
    # p-value of model
    p_value <- summary(lm)$coefficients[2,4]
    # R sqaure of model
    R_square  <- summary(lm)$r.squared
    # Slope of model
    Slope     <- lm$coefficients[2]
    # Intercept of model
    Intercept <- lm$coefficients[1]
  }else{
    #If less than n_WLS valid values in the slice, record NA
    p_value   <- NA
    R_square  <- NA
    Slope     <- NA
    Intercept <- NA
  }
  return(c(p_value,R_square,Slope,Intercept))
}

# WLS fit for each slice
# Input is raster_GS and original full-range df
# Output is a feature table including fixed level,p-value of the fit,
# R2 of the fit, Slope, and Intercept
# For all 5 VPD levels and 5 SM level
Slice_WLS <- function(raster_GS,df){
  # GS matrix
  GS_mat <- raster::as.matrix(raster_GS)
  # Get coordinates of SM and VPD levels
  raster_GS_df <- as.data.frame(raster_GS,xy=TRUE)
  SM_levels   <- unique(raster_GS_df$x)
  VPD_levels  <- unique(raster_GS_df$y)
  
  # Weight matrix for WLS
  w_mat <- WLS_w_mat(df)
  
  # Store all WLS features for this pixel
  WLS_features_all <- c()
  # Don't need to consider the first and last slice
  # Get features for rows first
  for(i in 2:(n_grid-1)){
    # Get a row slice at given VPD levels
    row_slice <- GS_mat[i,]
    # SM coordinates of each pixel
    coor <- SM_levels
    # This level is VPD level for this row
    level <- VPD_levels[i]
    # Get the WLS weight for this slice
    w <- w_mat[i,]
    # Get WLS features for this row slice
    features <- WLS_fit(row_slice,coor,w)
    # All features for this slice include level,p_value,R2,Slope,Intercept
    features <- c(level,features)
    WLS_features_all <- rbind(WLS_features_all,features)
  }
  # Do the same for columns
  for(i in 2:(n_grid-1)){
    # Get a column slice at given SM levels
    col_slice <- GS_mat[,i]
    # VPD coordinates of each pixel
    coor <- VPD_levels
    # This level is SM level for this column
    level <- SM_levels[i]
    # Get the WLS weight for this slice
    w <- w_mat[,i]
    # Get WLS features for this column slice
    features <- WLS_fit(col_slice,coor,w)
    # All features for this slice include level,p_value,R2,Slope,Intercept
    features <- c(level,features)
    WLS_features_all <- rbind(WLS_features_all,features)
  }
  WLS_features_df <- as.data.frame(WLS_features_all)
  names(WLS_features_df) <- c("Level","p_value","R2","Slope","Intercept")
  # Fix means either fix VPD level or SM level
  WLS_features_df$Fix <- c(rep("VPD",n_grid-2),rep("SM",n_grid-2))
  rownames(WLS_features_df) <- NULL
  return(WLS_features_df)
}

# This function is to extract WLS slopes and intercepts, and fitness of WLS
# For each slice
# Input data is the original full-range df
WLS_features <- function(df){
  # Get normalized df
  df <- Normalize_df(df)
  # Apply GS filtration
  raster_GS <- GSfiltration(df)
  # Apply WLS to slices to get slice features
  WLS_features_df <- Slice_WLS(raster_GS,df)
  return(WLS_features_df)
}

# This function is to filter features table
# Input is WLS_features_df
Filter_WLS_feature <- function(features_df){
  # Remove slices with low R2
  # Remove slices with low R2
  features_df <- features_df %>% filter(R2 >= Slice_R2)
  # Remove significantly incorrect slices
  filter_1 <- which((features_df$p_value < 0.05) & 
                      (features_df$Fix == "VPD") &
                      (features_df$Slope < 0))
  filter_2 <- which((features_df$p_value < 0.05) & 
                      (features_df$Fix == "SM") &
                      (features_df$Slope > 0))
  features_df[filter_1,] <- NA
  features_df[filter_2,] <- NA
  
  features_df <- na.omit(features_df)
  
  return(features_df)
}

# Get uncertainty of contour line slope and theta using bootstrap
# Input is Fit0_df
SD_contour <- function(Fit0_df){
  set.seed(1)
  slope_all <- c()
  theta_all <- c()
  for(i in 1:num_iteration){
    # Randomly sample the same number of selected points
    size <- nrow(Fit0_df)
    sample_indices <- sample(1:size,size=size,replace=TRUE)
    sample_data <- Fit0_df[sample_indices,]
    # Fit a line using the sampled points
    lm <- lm(data=sample_data,VPD_cor~SM_cor)
    slope <- lm$coefficients[2]
    slope_all <- c(slope_all,slope)
    theta <- atan(slope)/pi*180
    theta_all <- c(theta_all,theta)
  }  
  # get sd of slope
  slope_sd <- sd(slope_all,na.rm = T)
  # get sd of theta
  theta_sd <- sd(theta_all,na.rm = T)
  out <- c(slope_sd,theta_sd)
  return(out)
}

# This function gets Fitted 0 for each slice
# And get contour line features
# Input is the feature table and the original full-range df
Contour_slope <- function(features_df,df){
  # Fit0
  features_df$Fit0 <- -features_df$Intercept/features_df$Slope
  # Remove Fit0 that are over extrapolated
  # Which means they cannot exceed the original range of SM and VPD
  # For conservative reason
  min_SM <- min(df$SM,na.rm=T)
  max_SM <- max(df$SM,na.rm=T)
  min_VPD <- min(df$VPD,na.rm=T)
  max_VPD <- max(df$VPD,na.rm=T)
  filter1 <- which((features_df$Fix=="VPD" & features_df$Fit0 > max_SM)|
                     (features_df$Fix=="VPD" & features_df$Fit0 < min_SM))
  filter2 <- which((features_df$Fix=="SM" & features_df$Fit0 > max_VPD)|
                     (features_df$Fix=="SM" & features_df$Fit0 < min_VPD))
  features_df$Fit0[filter1] <- NA
  features_df$Fit0[filter2] <- NA
  
  # Get coordinates of all these Fit0
  SM_cor <- c(features_df$Fit0[features_df$Fix == "VPD"],
              features_df$Level[features_df$Fix == "SM"])
  VPD_cor <- c(features_df$Level[features_df$Fix == "VPD"],
               features_df$Fit0[features_df$Fix == "SM"])
  
  Fit0_df <- data.frame(SM_cor,VPD_cor)
  Fit0_df <- na.omit(Fit0_df)
  
  # If number of Fit 0 < n_Fit0, slope = NA
  if(nrow(Fit0_df) < n_Fit0){
    Slope <- NA
    R2 <- NA
    Slope_sd <- NA
    theta_sd <- NA
    Intercept <- NA
    theta <- NA
  }else{
    # Fit a straight line for all Fit 0 points
    lm <- lm(data=Fit0_df,VPD_cor~SM_cor)
    Slope <- lm$coefficients[2]
    # Get theta as arctan(slope)/pi*180 in degree
    theta <- atan(Slope)/pi*180
    R2 <- summary(lm)$r.squared
    Slope_sd <- SD_contour(Fit0_df)[1]
    theta_sd <- SD_contour(Fit0_df)[2]
    Intercept <- lm$coefficients[1]
  }
  # Return Slope,theta,Slope_sd,R2,Intercept
  Out <- c(Slope,Slope_sd,theta,theta_sd,R2,Intercept)
  return(Out)
}
#######
# Main
#######
tryCatch({
  # Get Block_area
  Block_area <- st_sf(CONUS[[2]][Block_ID])
  Block_area <- raster(extent(Block_area),resolution = c(0.25,0.25),crs="+proj=longlat +datum=WGS84")
  
  df_list <- readRDS(paste0(Input_path,"Full_df_ls_Block_",Block_ID,".rds"))
  meta_df <- read.csv(paste0(Input_path,"Meta_df_Block_",Block_ID,".csv"))
  
  # Get the number of ESI!=1 for each df in the list
  # Stores all ESI!=1
  n_ESI_no1 <- c()
  for(i in 1:length(df_list)){
    df <- df_list[[i]]
    n_ESI_no1 <- c(n_ESI_no1,sum(df$ESI!=1))
  }
  # Only keep pixels with n(ESI!=1) >= 10
  Index <- which(n_ESI_no1>=10)
  
  # Initiate vector to store all features
  Contour_features_all <- c()
  for(arrayid in Index){
    # Get the df for this pixel
    full_df <- df_list[[arrayid]]
    # Get slice WLS features
    WLS_features_df <- WLS_features(full_df)
    # Apply filters to features
    WLS_features_df <- Filter_WLS_feature(WLS_features_df)
    # Get contour line features
    Contour_features <- Contour_slope(WLS_features_df,full_df)
    # Combine all contour features
    Contour_features_all <- rbind(Contour_features_all,Contour_features)
    print(paste("Complete",round(arrayid/max(Index)*100,4),"%"))
  }
  
  print("Making and filtering final Contour feature df")
  
  # Convert to data frame
  feature_df <- as.data.frame(Contour_features_all)
  # Add Index for later match with meta df
  feature_df$Index <- Index
  names(feature_df) <- c("Slope","Slope_sd","Theta","Theta_sd","R2","Intercept","index")
  rownames(feature_df) <- NULL
  
  # Filter Contour_features_df
  # Remove Slope NA
  # Remove Slope_sd > Slope_CV_th
  # Remove R2 < Slope_R2_th
  feature_df <- feature_df %>% 
    filter(!is.na(Slope)) %>%
    filter(Slope_sd<=Slope_CV_th) %>%
    filter(R2>=Slope_R2_th) %>%
    filter(Theta>theta_th)
  
  # Match this with meta df
  feature_df <- merge(feature_df,
                      meta_df[,-1],
                      by="index",
                      all.x=TRUE)
  
  # Calculate required threshold features
  # Get SM and VPD quantile thresholds
  # This means the corresponding VPD quantile when fixing SM at 25 quantile
  feature_df$SM_25th <- 25*feature_df$Slope+feature_df$Intercept
  feature_df$SM_50th <- 50*feature_df$Slope+feature_df$Intercept
  feature_df$SM_75th <- 75*feature_df$Slope+feature_df$Intercept
  
  # Similarly, this means the corresponding SM quantile when fixing VPD at 25 quantile
  feature_df$VPD_25th <- (25-feature_df$Intercept)/feature_df$Slope
  feature_df$VPD_50th <- (50-feature_df$Intercept)/feature_df$Slope
  feature_df$VPD_75th <- (75-feature_df$Intercept)/feature_df$Slope
  
  # Get SM and VPD thresholds at original scale
  VPD_range <- feature_df$VPD_max-feature_df$VPD_min
  # This means the corresponding absolute value of VPD threshold at the original scale when fixing SM at 25 quantile
  feature_df$SM_25_abs <- feature_df$SM_25th/100*VPD_range + feature_df$VPD_min
  feature_df$SM_50_abs <- feature_df$SM_50th/100*VPD_range + feature_df$VPD_min
  feature_df$SM_75_abs <- feature_df$SM_75th/100*VPD_range + feature_df$VPD_min
  
  SM_range <- feature_df$SM_max - feature_df$SM_min
  # This means the corresponding absolute value of SM threshold at the original scale when fixing VPD at 25 quantile
  feature_df$VPD_25_abs <- feature_df$VPD_25th/100*SM_range + feature_df$SM_min
  feature_df$VPD_50_abs <- feature_df$VPD_50th/100*SM_range + feature_df$SM_min
  feature_df$VPD_75_abs <- feature_df$VPD_75th/100*SM_range + feature_df$SM_min
  
  # Apply filters to thresholds
  # Remove any quantile thresholds less than 0 or greater than 100,
  # And absolute VPD thresholds when fixing SM levels less than 0 or greater than 6 kPa.
  # And absolute SM thresholds when fixing VPD levels less than 0 or greater than 1.
  
  feature_df$SM_25th[feature_df$SM_25th<0|feature_df$SM_25th>100|feature_df$SM_25_abs<0|feature_df$SM_25_abs>6] <- NA
  feature_df$SM_50th[feature_df$SM_50th<0|feature_df$SM_50th>100|feature_df$SM_50_abs<0|feature_df$SM_50_abs>6] <- NA
  feature_df$SM_75th[feature_df$SM_75th<0|feature_df$SM_75th>100|feature_df$SM_75_abs<0|feature_df$SM_75_abs>6] <- NA
  
  feature_df$SM_25_abs[feature_df$SM_25th<0|feature_df$SM_25th>100|feature_df$SM_25_abs<0|feature_df$SM_25_abs>6] <- NA
  feature_df$SM_50_abs[feature_df$SM_50th<0|feature_df$SM_50th>100|feature_df$SM_50_abs<0|feature_df$SM_50_abs>6] <- NA
  feature_df$SM_75_abs[feature_df$SM_75th<0|feature_df$SM_75th>100|feature_df$SM_75_abs<0|feature_df$SM_75_abs>6] <- NA
  
  feature_df$VPD_25th[feature_df$VPD_25th<0|feature_df$VPD_25th>100|feature_df$VPD_25_abs<0|feature_df$VPD_25_abs>1] <- NA
  feature_df$VPD_50th[feature_df$VPD_50th<0|feature_df$VPD_50th>100|feature_df$VPD_50_abs<0|feature_df$VPD_50_abs>1] <- NA
  feature_df$VPD_75th[feature_df$VPD_75th<0|feature_df$VPD_75th>100|feature_df$VPD_75_abs<0|feature_df$VPD_75_abs>1] <- NA
  
  feature_df$VPD_25_abs[feature_df$VPD_25th<0|feature_df$VPD_25th>100|feature_df$VPD_25_abs<0|feature_df$VPD_25_abs>1] <- NA
  feature_df$VPD_50_abs[feature_df$VPD_50th<0|feature_df$VPD_50th>100|feature_df$VPD_50_abs<0|feature_df$VPD_50_abs>1] <- NA
  feature_df$VPD_75_abs[feature_df$VPD_75th<0|feature_df$VPD_75th>100|feature_df$VPD_75_abs<0|feature_df$VPD_75_abs>1] <- NA
  
  if(nrow(feature_df)==0){
    stop("nrow of feature_df is 0")
  }
  
  print("Put feature values into rasters")
  # Put the values of these features into rasters, same resolution as other variables
  # Initialize a raster stack to store all feature rasters
  feature_raster <- Block_area
  
  for(i in 1:ncol(feature_df)){
    raster_tmp <- Block_area
    for(j in 1:nrow(feature_df)){
      # Index of the pixel
      p_index <- feature_df$index[j]
      raster_tmp[p_index] <- feature_df[j,i]
    }
    names(raster_tmp) <- names(feature_df)[i]
    feature_raster <- stack(feature_raster,raster_tmp)
    print(paste("Complete",round(i/ncol(feature_df),4)*100,"%"))
  }
  
  # Output this feature_raster
  saveRDS(feature_raster,paste0(Output_path,"features_Block_",Block_ID,".rds"))
  
  # This error function is for the whole code, when there is no df input    
},error=function(e){cat("ERROR: No valid df for this grid:",conditionMessage(e),"\n")})

print(paste("Complete Block",Block_ID))










