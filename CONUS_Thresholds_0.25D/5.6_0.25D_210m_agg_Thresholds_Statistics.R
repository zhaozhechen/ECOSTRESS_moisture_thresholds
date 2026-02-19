# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2025.2.12

# This code is to get alpha (theta) and thresholds statistics for manuscript

# AI classification: Accelerated dryland expansion under climate change
# Hyper-arid: AI < 0.05
# Arid: 0.05 <= AI < 0.2
# Semiarid: 0.2 <= AI < 0.5
# Dry subhumid: 0.5 <= AI < 0.65
# Humid: AI >= 0.65

# Group hyper-arid and arid together
# Arid: AI < 0.2
# Semiarid: 0.2 <= AI < 0.5
# Semihumid: 0.5 <= AI < 0.65
# Humid: AI >= 0.65

#########
# Global
#########
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)

library(raster)
library(sf)
library(dplyr)

# Input path of combined 0.25D rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")

############
# Functions
############

# This function is to get the average within-PFT variance and cross-PFT variance
# Input is the var name to be calculated
var_PFT <- function(var){
  # Within-PFT variance
  within_var <- df %>%
    group_by(PFT) %>%
    summarise(within_var = var(get(var)))
  # Average within-PFT variance
  within_var_mean <- mean(within_var$within_var)
  # Cross-PFT variance
  cross_var <- var(df[[var]])
  print(c(within_var_mean,cross_var))
}

#######
# Main
#######

# Pre-processing of dataset ===============================
# Mask all layers with Theta
raster_all <- mask(raster_all,mask = raster_all$Theta)

# Get a df of three metrics and AI, and PFT
df <- data.frame(Alpha = values(raster_all$Theta),
                 VPD50qt = values(raster_all$SM_50th),
                 SM50qt = values(raster_all$VPD_50th),
                 VPD50abs = values(raster_all$SM_50_abs),
                 SM50abs = values(raster_all$VPD_50_abs),
                 AI = values(raster_all$AI),
                 PFT = values(raster_all$LC))
df <- na.omit(df)
df$Type <- rep(NA,nrow(df))
# Get dryland type classification
df$Type[df$AI < 0.2] <- "1"
df$Type[df$AI >= 0.2 & df$AI < 0.5] <- "2"
df$Type[df$AI >= 0.5 & df$AI < 0.65] <- "3"
df$Type[df$AI >= 0.65] <- "4"

df$Type <- as.factor(df$Type)

levels(df$Type) <- c("Arid","Semiarid","Semihumid","Humid")

# Define Dryland or Humid region
df$Dry <- rep(NA,nrow(df))
df <- df %>%
  mutate(Dry = replace(Dry,AI < 0.2,"Arid")) %>%
  mutate(Dry = replace(Dry,AI >= 0.2,"Humid"))
df$Dry <- as.factor(df$Dry)

###################################################################
# Get thresholds and alpha statistics ==========================
# Proportion of alpha > 30 in arid regions ========
# Number of arid pixels
n_arid <- nrow(df[df$Type == "Arid",])
# Number of alpha > 30 in arid pixels
n_arid_alpha30up <- nrow(df[df$Alpha > 30 & df$Type == "Arid",])
# Proportion
n_arid_alpha30up/n_arid

# Proportion of alpha < 30 in humid and semihumid regions ========
# Number of humid pixels
n_humid <- nrow(df[df$Type == "Humid",])
# Number of alpha < 30 in humid pixels
n_humid_alpha30down <- nrow(df[df$Alpha < 30 & df$Type == "Humid",])
# Proportion
n_humid_alpha30down/n_humid

# Number of semihumid pixels
n_semihumid <- nrow(df[df$Type == "Semihumid",])
# Number of alpha < 30 in humid pixels
n_semihumid_alpha30down <- nrow(df[df$Alpha < 30 & df$Type == "Semihumid",])
# Proportion
n_semihumid_alpha30down/n_semihumid

# Median of alpha in different regions ==============
median(df$Alpha[df$Type == "Arid"])
median(df$Alpha[df$Type == "Semiarid"])
median(df$Alpha[df$Type == "Semihumid"])
median(df$Alpha[df$Type == "Humid"])

# Median of SM50 in different regions ==============
median(df$SM50abs[df$Type == "Arid"])
median(df$SM50abs[df$Type == "Semiarid"])
median(df$SM50abs[df$Type == "Semihumid"])
median(df$SM50abs[df$Type == "Humid"])

# Median of VPD50 in different regions =============
median(df$VPD50abs[df$Type == "Arid"])
median(df$VPD50abs[df$Type == "Semiarid"])
median(df$VPD50abs[df$Type == "Semihumid"])
median(df$VPD50abs[df$Type == "Humid"])


# Average within-PFT and cross-PFT variance ===============
var_PFT("Alpha")
var_PFT("VPD50qt")
var_PFT("SM50qt")
var_PFT("VPD50abs")
var_PFT("SM50abs")





