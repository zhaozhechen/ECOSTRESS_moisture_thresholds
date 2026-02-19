# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2025.6.23

# This code is to calculate statistics related to Geographical divisions

#########
# Global
#########

#dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
#dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
library(raster)
library(ncdf4)
library(sf)
library(ggplot2)
library(cowplot)
library(wesanderson)
library(dplyr)
library(tidyr)
library(ggpattern)

# Input path of projected frequency changing only SM or VPD
Input_path1 <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_SM_VPD/"
# Input path of historical frequency and changing both SM and VPD
Input_path2 <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_V2/"
# raster_all with Geographical information
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Figures_V3/raster_all_with_GEO.rds")

#######
# Main
#######

# Read in data and pre-processing ============================================
time_name <- "End"
ssp <- "ssp585"

# Historical frequency
Hist_freq <- readRDS(paste0(Input_path2,"Hist_obs_freq.rds"))$Freq
Hist_freq[Hist_freq == 1] <- NA
# Mask by theta
Hist_freq <- mask(Hist_freq,raster_all$Theta)

# Model mean frequency, changing both VPD and SM
Proj_freq <- readRDS(paste0(Input_path2,"Model_mean_freq_",time_name,"_",ssp,".rds"))$Freq_mean
# Mask by theta
Proj_freq_SM_VPD <- mask(Proj_freq,raster_all$Theta)
# Calculate delta_frequency,using both projected SM and VPD
delta_freq_SM_VPD <- Proj_freq_SM_VPD - Hist_freq
# The corresponding p-value raster
p_value_SM_VPD <- readRDS(paste0(Input_path2,"p_value_rasters.rds"))
p_value_SM_VPD <- p_value_SM_VPD[[which(names(p_value_SM_VPD) == paste0("p_value_",time_name,"_",ssp))]]

# Model mean frequency, changing only SM
Proj_freq <- readRDS(paste0(Input_path1,"Model_mean_",time_name,"_",ssp,"_SM_chg_freq.rds"))$Freq_mean
# Mask by theta
Proj_freq_SM <- mask(Proj_freq,raster_all$Theta)
# Calculate delta_frequency,using only projected SM
delta_freq_SM <- Proj_freq_SM - Hist_freq
# The corresponding p-value raster
p_value_SM <- readRDS(paste0(Input_path1,"p_value_rasters.rds"))$p_value_End_ssp585_SM_chg

# Model mean frequency, changing only VPD
Proj_freq <- readRDS(paste0(Input_path1,"Model_mean_",time_name,"_",ssp,"_VPD_chg_freq.rds"))$Freq_mean
# Mask by theta
Proj_freq_VPD <- mask(Proj_freq,raster_all$Theta)
# Calculate delta_frequency,using only projected VPD
delta_freq_VPD <- Proj_freq_VPD - Hist_freq
# The corresponding p-value raster
p_value_VPD <- readRDS(paste0(Input_path1,"p_value_rasters.rds"))$p_value_End_ssp585_VPD_chg

# Compound effect of VPD + SM, which is delta_freq_SM_VPD - (delta_freq_SM + delta_freq_VPD)
delta_freq_compound <- delta_freq_SM_VPD - (delta_freq_SM + delta_freq_VPD)
p_value_compound <- p_value_SM_VPD

# Do statistics -----------------
# Get mean frequency induced by VPD in West North Central (region_id = 8)
stat1 <- mean(delta_freq_VPD[raster_all$division == 8])
# Mean Historical frequency in West North Central
stat2 <- mean(Hist_freq[raster_all$division == 8])

# Get mean frequency induced by VPD in West South Central (region_id = 9)
stat3 <- mean(delta_freq_VPD[raster_all$division == 9],na.rm=TRUE)
# Mean Historical frequency in West South Central
stat4 <- mean(Hist_freq[raster_all$division == 9],na.rm=TRUE)

# Get mean frequency induced by SM in the West US (region_id = 6 and 4)
stat5 <- mean(delta_freq_SM[raster_all$division == 6|raster_all$division == 4],na.rm=TRUE)
# Mean historical frequency in the West US (region_id = 6 and 4)
stat6 <- mean(Hist_freq[raster_all$division == 6|raster_all$division == 4],na.rm=TRUE)

# Get mean frequency induced by SM in Midwest (region_id = 8 and 1)
stat7 <- mean(delta_freq_SM[raster_all$division == 8|raster_all$division == 1],na.rm=TRUE)
# Mean historical frequency in Midwest
stat8 <- mean(Hist_freq[raster_all$division == 8|raster_all$division == 1],na.rm=TRUE)

# Combined mean frequency induced by both SM and VPD in CONUS
stat9 <- mean(values(delta_freq_SM_VPD),na.rm=TRUE)
# Mean historical frequency in CONUS
stat10 <- mean(values(Hist_freq),na.rm=TRUE)
stat9/stat10
