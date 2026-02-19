# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.6.27

# This code is to calculate spatial correlations between delta_frequency vs
# alpha, delta_VPD, and delta_SM

#########
# Global
#########

dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
library(raster)
library(ncdf4)
library(sf)
library(ggcorrplot)
library(sensemakr)

# Input path for historical observation
hist_obs_gs_mean <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Hist_obs/Hist_obs_gs_mean.rds")
# Input path for adjusted CMIP projection
CMIP_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/CMIP6_adjusted/"
# Input path of frequency and delta_VPD rasters
Freq_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_V2/"
# Input path of combined 0.25D rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")

#######
# Main
#######

time_name <- "End"
ssp <- "ssp585"

# Projected CMIP VPD
CMIP_VPD <- readRDS(paste0(CMIP_path,"Mean_models_gs_mean_",time_name,"_",ssp,"_daily_maxVPD.rds"))
# Projected CMIP SM
CMIP_SM <- readRDS(paste0(CMIP_path,"Mean_models_gs_mean_",time_name,"_",ssp,"_daily_SM.rds"))
# Change in VPD
delta_VPD <- CMIP_VPD - hist_obs_gs_mean$daily_maxVPD
# Change in SM
delta_SM <- CMIP_SM - hist_obs_gs_mean$daily_SM

# Historical frequency
Hist_freq <- readRDS(paste0(Freq_path,"Hist_obs_freq.rds"))$Freq
# Model mean frequency
Proj_freq <- readRDS(paste0(Freq_path,"Model_mean_freq_",time_name,"_",ssp,".rds"))$Freq_mean
# Change in frequency
Delta_freq <- Proj_freq - Hist_freq
# p-value for delta_frequency
p_value <- readRDS(paste0(Freq_path,"p_value_rasters.rds"))
p_value <- p_value[[which(names(p_value) == paste0("p_value_",time_name,"_",ssp))]]
# Only keep p-value <=0.05
p_value[p_value>0.05] <- NA
Delta_freq <- mask(Delta_freq,p_value)

# Alpha
Alpha <- raster_all$Theta

# Get a full df for correlation matrix
df <- data.frame(Delta_freq = values(Delta_freq),
                 Delta_SM = values(delta_SM),
                 Delta_VPD = values(delta_VPD),
                 Alpha = values(Alpha),
                 Delta_VPD_SM = -values(delta_SM * delta_VPD))
                 #alpha_Delta_VPD = values(delta_VPD * Alpha),
                 #alpha_Delta_SM = values(delta_SM * Alpha),
                 #alpha_Delta_VPD_SM = values(delta_SM * Alpha * delta_VPD))
df <- na.omit(df)
write.csv(df,paste0("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/delta_freq_End_ssp585.csv"))
# Output this df and use 
# D:\OneDrive - The Ohio State University\Research\ECOSTRESS\Results\CMIP_Projection\CMIP Projection Final 3 Figures\delta_freq_End_ssp585_CM.R 
# to make the plot because cannot install GGally on server

# Do a MLR =========================
# Full model
lm <- lm(data=df,Delta_freq~.)
summary(lm)

# Get partial R2 ===================
partial_r2 <- partial_r2(lm)[2:8]
# Normalize partial r2
partial_r2_norm <- partial_r2/sum(partial_r2)
partial_r2_df <- data.frame(partial_r2 = partial_r2_norm,
                            var = names(partial_r2_norm))
ggplot(partial_r2_df,aes(x=partial_r2,y=reorder(var,partial_r2)))+
  geom_bar(stat='identity')+
  labs(y="",x="Contribution")
  






