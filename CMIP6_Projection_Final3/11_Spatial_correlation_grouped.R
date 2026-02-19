# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.8.1

# This code is to calculate spatial correlations between delta_frequency vs
# alpha, delta_VPD, and delta_SM, grouped by binned alpha

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
library(dplyr)
library(tidyr)

# Input path for historical observation
hist_obs_gs_mean <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Hist_obs/Hist_obs_gs_mean.rds")
# Input path for adjusted CMIP projection
CMIP_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/CMIP6_adjusted/"
# Input path of frequency and delta_VPD rasters
Freq_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Frequency_stacks_V2/"
# Input path of combined 0.25D rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")

###########
# Functions
###########
# Define a function to calculate correlations within each bin
calculate_correlations <- function(data) {
  result <- data %>%
    summarise(
      count = n(),
      cor_Delta_freq_Delta_SM = cor(Delta_freq, Delta_SM, use = "complete.obs"),
      cor_Delta_freq_Delta_VPD = cor(Delta_freq, Delta_VPD, use = "complete.obs"),
      cor_Delta_freq_Delta_VPD_SM = cor(Delta_freq, Delta_VPD_SM, use = "complete.obs")
    )
  return(result)
}

my_theme <- theme(
  axis.line=element_line(color="black"),
  panel.background = element_blank(),
  #text = element_text(size=100),
  panel.border = element_rect(colour="black",fill=NA),
  legend.key = element_blank(),
  #legend.key.size = unit(6,"cm"),
  #aspect.ratio = 1/1,
  #legend.key.size = unit(0.3,'cm'),
  legend.text = element_text(size=16),
  plot.title = element_text(size=16),
  legend.title = element_text(size=16),
  #plot.margin = margin(80,80,80,80),
  axis.text.x = element_text(size=16),
  axis.text.y = element_text(size=16),
  axis.title = element_text(size=16)
  #axis.ticks = element_text(size=16)
  #legend.margin = margin(0,80,0,80)
  #axis.ticks.length = unit(-0.8,"cm")
)

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
df <- na.omit(df)

# Bin alpha (every 10 degree)
bin_width <- 10

# Bin Alpha with the bin_width
df <- df %>%
  mutate(Alpha_bin = cut(Alpha,breaks=seq(min(Alpha,na.rm=T),
                                          max(Alpha,na.rm=T),
                                          by = bin_width),
                         include.lowest = TRUE))
# Apply the function to calculate correlation to each bin
# Get the small bins (with small sample sizes)
small_bins <- df %>%
  group_by(Alpha_bin) %>%
  do(calculate_correlations(.)) %>%
  ungroup() %>%
  filter(count < 350) %>%
  pull(Alpha_bin)

# Combine the last four small bins into one category if there are enough small bins
combined_bin_label <- "Combined_Last_Four_Bins"

# Combine small bins
df <- df %>%
  mutate(Alpha_bin = ifelse(Alpha_bin %in% small_bins, combined_bin_label, as.character(Alpha_bin)))

# Recalculate statistics with the combined bins
final_statistics_by_bin <- df %>%
  group_by(Alpha_bin) %>%
  do(calculate_correlations(.))

# Rename the levels of the bins
final_statistics_by_bin$levels <- c(2,4,3,5,6,7,8,1)
final_statistics_by_bin <- final_statistics_by_bin %>%
  pivot_longer(cols = starts_with("cor_"),
               names_to = "correlation_type",
               values_to = "correlation_value")
final_statistics_by_bin$correlation_value <- abs(final_statistics_by_bin$correlation_value)

# Make plots =================================================
x_lables <- c("1" = "[-20,-10]",
              "2" = "(-10,0]",
              "3" = "(0,10]",
              "4" = "(10,20]",
              "5" = "(20,30]",
              "6" = "(30,40]",
              "7" = "(40,50]",
              "8" = "(50,90]")

# Count df
counts_df <- final_statistics_by_bin %>%
  filter(correlation_type == "cor_Delta_freq_Delta_VPD")

g <- ggplot(final_statistics_by_bin,
            aes(x=levels,y=correlation_value,fill=correlation_type))+
  geom_bar(stat="identity",position="dodge",color="black")+
  scale_x_continuous(breaks = c(1:8),
    labels = x_lables)+
  labs(x=expression(alpha~"("~degree~")"),y="|r|",fill="")+
  scale_fill_manual(values = c("#C8AF81","#F4931A","#A5D0E4"),
    labels = c("cor_Delta_freq_Delta_SM" = expression(r~"("~Delta~Frequency~"~"~Delta~SM~")"),
                               "cor_Delta_freq_Delta_VPD" = expression(r~"("~Delta~Frequency~"~"~Delta~VPD~")"),
                               "cor_Delta_freq_Delta_VPD_SM" = expression(r~"("~Delta~Frequency~"~"~Delta~SM~"*"~Delta~VPD~")")))+
  geom_text(data=counts_df,aes(x=levels,y=correlation_value,label = paste0("n =",count)),vjust=-1.5)+
  my_theme+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Adjust size and angle
    legend.position = c(0.8,0.9)
  )











