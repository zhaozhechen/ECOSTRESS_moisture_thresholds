# Functions copied from Validation_thresholds_20250202.R
# Compatibility change: summarise -> transmute in Normalize_df preserves the
# original multi-row output under current dplyr; calculations are unchanged.
# Original author: Zhaozhe Chen. Original server code remains unchanged.

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
# Z = ESI - median(ESI)
# SM and VPD are quantiles of full-range
# Input is original full-range AMF or RS df
Normalize_df <- function(data_df){
  # Get only ESI, VPD, and SM
  raw_df <- subset(data_df,select=c(ESI_daily,Daily_max_VPD,SM_daily))
  # Normalize ESI, VPD and SM to full-range quantiles
  raw_df <- raw_df %>% transmute(ESI = rescale(ESI_daily)*100,
                                 VPD = rescale(Daily_max_VPD)*100,
                                 SM = rescale(SM_daily)*100)
  raw_df$Z <- raw_df$ESI - median(raw_df$ESI,na.rm=T)
  # Remove NA, at this step, 
  # SM and VPD are kept as their full-range quantiles
  df <- na.omit(raw_df)
  return(df)
}

# Weighted least square regression (WLS)
# Get weight matrix
# Input is original full-range df
WLS_w_mat <- function(data_df){
  # Normalize df
  df <- Normalize_df(data_df)
  # Z values are 1 for all points
  data_matrix   <- cbind(df$SM,df$VPD,rep(1,nrow(df)))
  extent        <- extent(data_matrix[,1:2])
  data_raster   <- raster(extent,ncol=n_grid,nrow=n_grid)
  # Sum of numbers of data points in each pixel
  data_raster   <- raster::rasterize(data_matrix[,1:2],data_raster,data_matrix[,3],fun=sum)
  data_raster   <- gaussian.smooth(data_raster,sigma=sigma,n=n_k,type = "sum")
  w_mat <- raster::as.matrix(data_raster)
  return(w_mat)
}

# Fit SLR and output slope, intercept, R2, and p-value
# Actually WLS fit after adding weight matrix
SLR_fit <- function(slice,coor,w){
  # Only proceed if there are at least n_SLR in the slice
  if(sum(!is.na(slice))>=n_SLR){
    SLR_df <- data.frame(coor,slice)
    # WLS fit
    lm <- lm(data = SLR_df,slice~coor,weights = w)
    # p-value of model
    p_value <- summary(lm)$coefficients[2,4]
    # R sqaure of model
    R_square  <- summary(lm)$r.squared
    # Slope of model
    Slope     <- lm$coefficients[2]
    # Intercept of model
    Intercept <- lm$coefficients[1]
  }else{
    #If less than n_SLR valid values in the slice, record NA
    p_value   <- NA
    R_square  <- NA
    Slope     <- NA
    Intercept <- NA
  }
  return(c(p_value,R_square,Slope,Intercept))
}

# SLR fit for each slice. 
# Input is raster_GS and original full-range df
# Output is a features table including fixed level, slopes, intercepts, R2, and p-value
# For all 5 VPD levels and 5 SM levels
Slice_SLR <- function(raster_GS,data_df){
  # GS matrix
  GS_mat <- raster::as.matrix(raster_GS)
  # Get coordinates of SM and VPD levels
  raster_GS_df <- as.data.frame(raster_GS,xy=TRUE)
  SM_levels   <- unique(raster_GS_df$x)
  VPD_levels  <- unique(raster_GS_df$y)
  
  # Weight matrix for WLS
  w_mat <- WLS_w_mat(data_df)
  
  # Store all SLR features for this Site
  features_all <- c()
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
    
    features <- SLR_fit(row_slice,coor,w)
    # Features include level, p_value,R2,Slope,Intercept
    features <- c(level,features)
    features_all <- rbind(features_all,features)
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
    
    features <- SLR_fit(col_slice,coor,w)
    # Features include level, p_value,R2,Slope,Intercept
    features <- c(level,features)
    features_all <- rbind(features_all,features)
  }
  features_df <- as.data.frame(features_all)
  names(features_df) <- c("Level","p_value","R2","Slope","Intercept")
  features_df$Fix <- c(rep("VPD",n_grid-2),rep("SM",n_grid-2))
  rownames(features_df) <- NULL
  return(features_df)
}

# This function is to extract SLR slopes and intercepts for each slice
# Input data is the original full-range df
SLR_features <- function(data_df){
  # Get normalized df
  df <- Normalize_df(data_df)
  # Apply GS filtration
  raster_GS <- GSfiltration(df)
  # Apply SLR to slices
  features_df <- Slice_SLR(raster_GS,data_df)
  return(features_df)
}

# This function is to filter features table
# Input is features_df
Filter_SLR_feature <- function(features_df){
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

# This function is to calculate theta from Slope
# Input is k (slope)
k2alpha <- function(k){
  theta <- atan(k)*180/pi
  return(theta)
}

# This function is to remove outlier first, then calculate sd
# Input is the vector of numbers to calculate
rmoutlier_sd <- function(vector){
  # Compute Q1,Q3,and IQR
  Q1 <- quantile(vector,0.25,na.rm=T)
  Q3 <- quantile(vector,0.75,na.rm=T)
  IQR <- Q3 - Q1
  # Define lower and upper bounds for outliers
  l_bound <- Q1 - 1.5*IQR
  u_bound <- Q3 + 1.5*IQR
  # Remove outliers
  filtered_vector <- vector[vector >= l_bound & vector <= u_bound]
  std <- sd(filtered_vector,na.rm=T)
  return(std)
}

# Get uncertainty of contour line slope using bootstrap
# Input is Fit0_df
SD_contour <- function(Fit0_df,VPD_range,VPD_min,SM_range,SM_min){
  set.seed(1)
  slope_all <- c()
  theta_all <- c()
  VPD50_qt_all <- c()
  SM50_qt_all <- c()
  VPD50_abs_all <- c()
  SM50_abs_all <- c()
  for(i in 1:num_iteration){
    # Randomly sample the same number of selected points
    size <- nrow(Fit0_df)
    sample_indices <- sample(1:size,size=size,replace=TRUE)
    sample_data <- Fit0_df[sample_indices,]
    # Fit a line using the sampled points
    lm <- lm(data=sample_data,VPD_cor~SM_cor)
    slope <- lm$coefficients[2]
    Intercept <- lm$coefficients[1]
    slope_all <- c(slope_all,slope)
    # Calculate theta in degree
    theta <- k2alpha(slope)
    theta_all <- c(theta_all,theta)
    # Calculate SM50 and VPD50 in quantiles
    VPD50_qt <- 50*slope + Intercept
    # If quantile is greater than 100, then equal to 100
    # If quantile is less than 0, then equal to 0
    #if(!is.na(VPD50_qt) & VPD50_qt > 100){VPD50_qt <- 100}
    #if(!is.na(VPD50_qt) & VPD50_qt < 0){VPD50_qt <- 0}

    VPD50_qt_all <- c(VPD50_qt_all,VPD50_qt)
    SM50_qt <- (50-Intercept)/slope
    #if(!is.na(SM50_qt) & SM50_qt > 100){SM50_qt <- 100}
    #if(!is.na(SM50_qt) & SM50_qt < 0){SM50_qt <- 0}
    SM50_qt_all <- c(SM50_qt_all,SM50_qt)

    # Calculate SM50 and VPD50 in absolute values
    VPD50_abs <- VPD50_qt/100*VPD_range + VPD_min
    VPD50_abs_all <- c(VPD50_abs_all,VPD50_abs)
    SM50_abs <- SM50_qt/100*SM_range + SM_min
    SM50_abs_all <- c(SM50_abs_all,SM50_abs)
  }  
  # get sd of slope
  slope_sd <- rmoutlier_sd(slope_all)
  # get sd of theta
  theta_sd <- rmoutlier_sd(theta_all)
  # get sd of SM50_qt and VPD50_qt
  VPD50_qt_sd <- rmoutlier_sd(VPD50_qt_all)
  SM50_qt_sd <- rmoutlier_sd(SM50_qt_all)
  # get sd of SM50_abs and VPD50_abs
  VPD50_abs_sd <- rmoutlier_sd(VPD50_abs_all)
  SM50_abs_sd <- rmoutlier_sd(SM50_abs_all)
  return(c(slope_sd,theta_sd,VPD50_qt_sd,SM50_qt_sd,VPD50_abs_sd,SM50_abs_sd))
}

# This function gets Fitted 0 for each slice
# And get contour line features
# Input is the feature table
Contour_slope <- function(features_df,full_df){
  # Fit0
  features_df$Fit0 <- -features_df$Intercept/features_df$Slope
  # Remove Fit0 that are over extrapolated
  # Get maximum and minimum SM and maximum and minimum VPD as boundaries
  df <- Normalize_df(full_df)
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
  
  # Get VPD and SM range and min, for the calculation of SM50 and VPD50 in their original scale
  VPD_range <- max(full_df$Daily_max_VPD,na.rm=T) - min(full_df$Daily_max_VPD,na.rm=T)
  VPD_min <- min(full_df$Daily_max_VPD,na.rm=T)
  SM_range <- max(full_df$SM_daily,na.rm=T) - min(full_df$SM_daily,na.rm=T)
  SM_min <- min(full_df$SM_daily,na.rm=T)
  
  # If number of Fit 0 < n_Fit0, slope = NA
  if(nrow(Fit0_df) < n_Fit0){
    Slope <- NA
    R2 <- NA
    Slope_sd <- NA
    Intercept <- NA
    theta <- NA
    theta_sd <- NA
    VPD50_qt <- NA
    VPD50_qt_sd <- NA
    SM50_qt <- NA
    SM50_qt_sd <- NA
    VPD50_abs <- NA
    VPD50_abs_sd <- NA
    SM50_abs <- NA
    SM50_abs_sd <- NA
  }else{
    # Fit a straight line for all Fit 0 points
    lm <- lm(data=Fit0_df,VPD_cor~SM_cor)
    Slope <- lm$coefficients[2]
    R2 <- summary(lm)$r.squared
    Slope_sd <- SD_contour(Fit0_df,VPD_range,VPD_min,SM_range,SM_min)[1]
    Intercept <- lm$coefficients[1]
    theta <- k2alpha(Slope)
    theta_sd <- SD_contour(Fit0_df,VPD_range,VPD_min,SM_range,SM_min)[2]
    # SM50_qt means the SM threshold when fixing VPD at the 50 quantile
    VPD50_qt <- 50*Slope + Intercept
    # If quantile is greater than 100, then equal to 100
    # If quantile is less than 0, then equal to 0
    if(!is.na(VPD50_qt) & VPD50_qt > 100){VPD50_qt <- 100}
    if(!is.na(VPD50_qt) & VPD50_qt < 0){VPD50_qt <- 0}
    VPD50_qt_sd <- SD_contour(Fit0_df,VPD_range,VPD_min,SM_range,SM_min)[3]
    # Thresholds in quantiles
    # VPD50_qt means the VPD threshold when fixing SM at the 50 quantile
    SM50_qt <- (50-Intercept)/Slope
    if(!is.na(SM50_qt) & SM50_qt > 100){SM50_qt <- 100}
    if(!is.na(SM50_qt) & SM50_qt < 0){SM50_qt <- 0}
    SM50_qt_sd <- SD_contour(Fit0_df,VPD_range,VPD_min,SM_range,SM_min)[4]
    # Thresholds in absolute values 
    # SM50_abs means the SM thresholds when fixing VPD at the 50 quantile
    VPD50_abs <- VPD50_qt/100*VPD_range + VPD_min
    VPD50_abs_sd <- SD_contour(Fit0_df,VPD_range,VPD_min,SM_range,SM_min)[5]
    SM50_abs <- SM50_qt/100*SM_range + SM_min
    SM50_abs_sd <- SD_contour(Fit0_df,VPD_range,VPD_min,SM_range,SM_min)[6]
  }

  # Return Slope,theta,Slope_sd,R2,Intercept
  out <- c(Slope,Slope_sd,theta,theta_sd,R2,Intercept,
           VPD50_qt,VPD50_qt_sd,SM50_qt,SM50_qt_sd,
           VPD50_abs,VPD50_abs_sd,SM50_abs,SM50_abs_sd)
  names(out) <- c("Slope","Slope_sd","theta","theta_sd","R2","Intercept",
                  "VPD50_qt","VPD50_qt_sd","SM50_qt","SM50_qt_sd",
                  "VPD50_abs","VPD50_abs_sd","SM50_abs","SM50_abs_sd")
  return(out)
}

# This function is to get all contour line features and SM and VPD range and min
# Input is the original full df
SM_VPD <- function(full_df){
  # Get slice WLS features
  features_df <- SLR_features(full_df)
  # Apply filters to features
  features_df <- Filter_SLR_feature(features_df)
  # Get contour line features
  features_df <- Contour_slope(features_df,full_df)
  features_df <- as.data.frame(t(features_df))
  return(features_df)
}



