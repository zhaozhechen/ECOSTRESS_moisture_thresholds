# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2026.8.19

# This code is to compare different meshgrid sizes for quantifying
# remote sensing-derived thresholds at US-A32 and US-CF3

# Reference code: /fs/ess/PAS2204/Code/CONUS_Threshold/Validation_Final_2/Validation_GS_Plots_20240327.R

# ---- Global --------
library(raster)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)
library(RColorBrewer)

# Paths to full-range df
Input_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/00_Data/Processed/Full_range_df/"
# Path to output figures
Output_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/02_Results/Figures"
# Path to output tables
Table_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/02_Results/Tables"

# Sites and meshgrid sizes to test
Site_ls <- c("US-A32","US-CF3")
n_grid_ls <- c(3,5,7,9,11)

# GS parameters
n_grid  <- 7
n_k     <- 3
sigma   <- 1

# SLR and contour line parameters
# At least n_SLR number of points in each slice to fit a slice SLR
n_SLR <- 4
# Slice R2 at least greater or equal to Slice_R2 to Fit0
Slice_R2 <- 0.8
# Least number of Fit0 point to fit a contour slope
n_Fit0 <- 3

# ---- Functions --------
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
    # R square of model
    R_square  <- summary(lm)$r.squared
    # Slope of model
    Slope     <- lm$coefficients[2]
    # Intercept of model
    Intercept <- lm$coefficients[1]
  }else{
    # If less than n_SLR valid values in the slice, record NA
    p_value   <- NA
    R_square  <- NA
    Slope     <- NA
    Intercept <- NA
  }
  return(c(p_value,R_square,Slope,Intercept))
}

# SLR fit for each slice
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

# This function gets fitted 0 for each slice
# And gets contour line features
Contour_slope <- function(features_df,data_df){
  # Fit0
  features_df$Fit0 <- -features_df$Intercept/features_df$Slope
  # Remove Fit0 that are over extrapolated
  df <- Normalize_df(data_df)
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
    Intercept <- NA
  }else{
    # Fit a straight line for all Fit 0 points
    lm <- lm(data=Fit0_df,VPD_cor~SM_cor)
    Slope <- lm$coefficients[2]
    R2 <- summary(lm)$r.squared
    Intercept <- lm$coefficients[1]
  }
  out <- list(Fit0_df=Fit0_df,
              Contour_feature=c(Slope=unname(Slope),
                                  R2=unname(R2),
                                  Intercept=unname(Intercept)))
  return(out)
}

my_theme <- theme(
  axis.line=element_line(color="black"),
  panel.background = element_blank(),
  text = element_text(size=14),
  panel.border = element_rect(colour="black",fill=NA),
  legend.key = element_blank(),
  aspect.ratio = 1/1,
  legend.title=element_text(size=12),
  axis.text = element_text(size=14),
  legend.text = element_text(size=14)
)

# This function is to get range of GS plot
GS_range <- function(data_df){
  # Get normalized df
  df <- Normalize_df(data_df)
  # Apply GS filtration
  raster_GS <- GSfiltration(df)
  # Normalize values on the GS plot as delta_ESI/sd(delta_ESI)
  values(raster_GS) <- values(raster_GS)/sd(values(raster_GS),na.rm=T)
  GS_df <- as.data.frame(raster_GS,xy=TRUE)
  # Get range of heat map
  SM_min <- min(GS_df$x)
  SM_max <- max(GS_df$x)
  VPD_min <- min(GS_df$y)
  VPD_max <- max(GS_df$y)

  # Output these ranges
  range <- c(SM_min,SM_max,VPD_min,VPD_max)
  return(range)
}

# Make GS plots with the manuscript styling
GS_plot <- function(Contour_feature,data_df,range,status,arrayid,n_grid){
  # Slope of the contour line
  k <- Contour_feature["Slope"]
  # Intercept of the contour line
  b <- Contour_feature["Intercept"]
  # Get normalized df
  df <- Normalize_df(data_df)
  # Apply GS filtration
  raster_GS <- GSfiltration(df)
  # Retain the original mesh boundaries before smoothing the raster
  grid_extent <- extent(raster_GS)
  grid_x <- seq(xmin(grid_extent),xmax(grid_extent),length.out=n_grid+1)
  grid_y <- seq(ymin(grid_extent),ymax(grid_extent),length.out=n_grid+1)
  # Normalize values on the GS plot as delta_ESI/sd(delta_ESI)
  values(raster_GS) <- values(raster_GS)/sd(values(raster_GS),na.rm=T)
  # Make this raster smoother
  # Use as a reference
  extent      <- extent(raster_GS)
  raster_rf   <- raster(extent,ncol=201,nrow=201)

  raster_GS <- resample(raster_GS,raster_rf,method="bilinear")
  GS_df <- as.data.frame(raster_GS,xy=TRUE)
  GS_df$layer[GS_df$layer>=2] <- 2

  # Normalize original data points
  df$Z <- df$Z/sd(df$Z)
  RdBu_color <- RColorBrewer::brewer.pal(11, "RdBu")

  # Manually set ticks as in the manuscript figure
  if(arrayid == 1){
    x_min <- 10
    x_max <- 80
    y_min <- 30
    y_max <- 80
  }else{
    x_min <- 30
    x_max <- 80
    y_min <- 30
    y_max <- 70
  }

  g <- ggplot(data=df)+
    geom_tile(data=GS_df,aes(x,y,fill=layer))+
    geom_vline(xintercept=grid_x[-c(1,length(grid_x))],
               color="grey35",linetype=2,linewidth=0.25,alpha=0.65)+
    geom_hline(yintercept=grid_y[-c(1,length(grid_y))],
               color="grey35",linetype=2,linewidth=0.25,alpha=0.65)+
    geom_point(aes(x=SM,y=VPD,fill=Z),size=2.3,color="black",pch=21,stroke = 0.3)+
    scale_fill_gradient2(na.value="white",
                         low = "#DF0101",
                         high = RdBu_color[10],
                         midpoint = 0,
                         guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"),
                         limits = c(-2,2))+
    my_theme+
    scale_x_continuous(breaks=seq(x_min,x_max,10),limits = c(range[1],range[2]))+
    scale_y_continuous(breaks=seq(y_min,y_max,10),limits = c(range[3],range[4]))+
    labs(x="SM quantile",y="VPD quantile",fill="ESI")

  # Only add a threshold line when all original criteria are met
  if(!is.na(k)){
    g <- g+geom_abline(slope=k,intercept=b,color="black",linetype=2,linewidth=0.8)
  }else{
    g <- g+annotate("label",
                    x=mean(range[1:2]),
                    y=range[4]-0.08*(range[4]-range[3]),
                    label=status,
                    size=3,
                    fill="white")
  }

  return(g)
}

# ---- Main --------
# Get the same plotting ranges used for the 7 x 7 manuscript panels
Site_range_ls <- list()
n_grid <- 7

for(arrayid in 1:length(Site_ls)){
  Site_ID <- Site_ls[arrayid]
  # Read in AMF full-range df only to reproduce the original plotting range
  AMF_df <- read.csv(paste0(Input_path,"AMF/AMF_Full_range_df_",Site_ID,".csv"))
  names(AMF_df)[names(AMF_df)=="VPD_daily_max"] <- "Daily_max_VPD"
  # Read in RS full-range df
  RS_df <- read.csv(paste0(Input_path,"RS/RS_Full_range_df_",Site_ID,".csv"))

  AMF_range <- GS_range(AMF_df)
  RS_range  <- GS_range(RS_df)
  Site_range_ls[[arrayid]] <- c(max(AMF_range[1],RS_range[1])-5,
                                min(AMF_range[2],RS_range[2]),
                                max(AMF_range[3],RS_range[3]),
                                min(AMF_range[4],RS_range[4])-5)
}

# Initialize outputs
g_ls <- list()
summary_df <- data.frame()
legend <- NULL

for(arrayid in 1:length(Site_ls)){
  Site_ID <- Site_ls[arrayid]
  RS_df <- read.csv(paste0(Input_path,"RS/RS_Full_range_df_",Site_ID,".csv"))
  range <- Site_range_ls[[arrayid]]

  for(grid_id in 1:length(n_grid_ls)){
    n_grid <- n_grid_ls[grid_id]

    # Get slice SLR features
    RS_features_all <- SLR_features(RS_df)
    n_available_slices <- sum(!is.na(RS_features_all$Slope))

    # Apply filters to features
    RS_features <- Filter_SLR_feature(RS_features_all)

    # Get contour line slopes
    Contour_out <- Contour_slope(RS_features,RS_df)
    RS_contour_features <- Contour_out$Contour_feature
    n_Fit0_points <- nrow(Contour_out$Fit0_df)

    # Apply the same final boundary criteria used in the manuscript analysis
    theta <- atan(RS_contour_features["Slope"])/pi*180

    if(n_grid <= 5){
      status <- "Insufficient cells per slice"
    }else if(n_Fit0_points < n_Fit0){
      status <- "Insufficient threshold points"
    }else if(RS_contour_features["R2"] < 0.5){
      status <- "Threshold points are inconsistent"
    }else if(theta <= -20){
      status <- "Boundary failed angle criterion"
    }else{
      status <- "Threshold estimated"
    }

    # Do not draw a line when the threshold does not pass all criteria
    Plot_contour_features <- RS_contour_features
    if(status != "Threshold estimated"){
      Plot_contour_features[c("Slope","Intercept")] <- NA
    }

    suppressWarnings(g <- GS_plot(Plot_contour_features,
                                  RS_df,
                                  range,
                                  status,
                                  arrayid,
                                  n_grid))

    if(is.null(legend) & n_grid==7){
      legend <- get_legend(g)
    }
    g <- g+theme(legend.position="none")
    g_ls[[paste(Site_ID,n_grid,sep="_")]] <- g

    summary_df <- rbind(summary_df,
                        data.frame(
                          Site_ID=Site_ID,
                          Meshgrid=paste0(n_grid,"x",n_grid),
                          Complete_observations=nrow(Normalize_df(RS_df)),
                          Available_slices=n_available_slices,
                          Retained_slices=nrow(RS_features),
                          Threshold_points=n_Fit0_points,
                          Slope=RS_contour_features["Slope"],
                          Theta=theta,
                          R2=RS_contour_features["R2"],
                          Intercept=RS_contour_features["Intercept"],
                          Status=status,
                          row.names=NULL
                        ))
  }
}

# Column headings
title_ls <- lapply(n_grid_ls,function(x){
  ggdraw()+draw_label(paste0(x," x ",x),fontface="plain",size=14)
})
title_row <- plot_grid(plotlist=title_ls,nrow=1)

# Plot rows
row_1 <- plot_grid(plotlist=g_ls[paste("US-A32",n_grid_ls,sep="_")],nrow=1,align="hv",axis="lrbt")
row_2 <- plot_grid(plotlist=g_ls[paste("US-CF3",n_grid_ls,sep="_")],nrow=1,align="hv",axis="lrbt")
main <- plot_grid(title_row,row_1,row_2,ncol=1,rel_heights=c(0.08,1,1))

# Row headings
row_heading <- plot_grid(
  NULL,
  ggdraw()+draw_label("US-A32",angle=90,size=14),
  ggdraw()+draw_label("US-CF3",angle=90,size=14),
  ncol=1,
  rel_heights=c(0.08,1,1)
)
main <- plot_grid(row_heading,main,nrow=1,rel_widths=c(0.035,1))

# Add the shared legend
g_all <- plot_grid(main,legend,nrow=1,rel_widths=c(1,0.08))

pdf(paste0(Output_path,"/Meshgrid_size_comparison.pdf"),
    height=7,
    width=18)
suppressWarnings(print(g_all))
dev.off()

png(paste0(Output_path,"/Meshgrid_size_comparison.png"),
    height=7,
    width=18,
    units="in",
    res=600)
suppressWarnings(print(g_all))
dev.off()

write.csv(summary_df,
          paste0(Table_path,"/Meshgrid_size_comparison.csv"),
          row.names=FALSE)

print(summary_df)
