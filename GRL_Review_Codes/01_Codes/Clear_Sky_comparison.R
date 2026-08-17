# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2026.8.17

# This code it to compare ESI values within lower-quantile and higher-quantile
# of net radiation groups, at the two example AMF sites (US-A32 and US-CF3)

# Reference code: /fs/ess/PAS2204/Code/Validation_ET_ESI_All_AMF/AMF_Full_range_df_20230930.R

# ---- Global --------
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(scales)

# Path to raw AMF dataset
AMF_path <- "/fs/ess/PAS2204/SharedData/AmeriFlux_All_Sites/"
# Path to full-range df of AMF and RS
Full_df_path <- "/fs/ess/PAS2204/Results/Validation_ET_ESI_All_AMF/Full_range_df/"

# Source plotting functions
source("/fs/ess/PAS2204/Code/GRL_Review_Codes/Plotting_functions.R")
# Sites to test
Site_ID_ls <- c("US-A32","US-CF3")

my_color <- brewer.pal(6,"Set2")
############
# Functions
############

Daily_mean <- function(data){
  DM <- matrix(data=data,nrow=window_size)
  DM <- colMeans(DM,na.rm=T)
  return(DM)
}

# Get required variable after QC
# Input is the variable name
Var_QC <- function(variable){
  varlist <- colnames(AMF)
  # Variable name in the dataset
  var <- varlist[grepl(variable,varlist)&!grepl("QC",varlist)]
  var <- AMF[var][,1]
  # QC for this variable
  var_QC <- varlist[grepl(variable,varlist)&grepl("QC",varlist)]
  if(length(var_QC!=0)){
    # Apply QC, only keep QC = 0
    var_QC <- AMF[var_QC][,1]
    var[var_QC!=0] <- NA
  }
  return(var)
}

# ------ Main -------
# Initialize a list to store figures
g_ls <- list()

for(arrayid in 1:length(Site_ID_ls)){
  Site_ID <- Site_ID_ls[arrayid]
  # Read in raw AMF data to get net radiation data ==========
  # Find the folder for this site
  folder_name <- dir(AMF_path)[grepl(Site_ID,dir(AMF_path))]
  folder_name <- paste(AMF_path,folder_name,sep="") 
  # Check if the data is half-hourly or hourly
  if(sum(grepl("SUBSET_HH_",dir(folder_name)))==1){
    # If there is half-hourly dataset
    # Window_size is 48
    window_size <- 48
    # Read the data
    file_name <- dir(folder_name)[grepl("SUBSET_HH_",dir(folder_name))]
  }else if(sum(grepl("SUBSET_HR_",dir(folder_name)))==1){
    # If there is hourly dataset
    # Window_size is 24
    window_size <- 24
    # Read the data
    file_name <- dir(folder_name)[grepl("SUBSET_HR_",dir(folder_name))]
  }
  
  AMF <- read.csv(paste(folder_name,file_name,sep="/"))
  
  # Check if the dataset has more than 4 years
  n_years <- nrow(AMF)/window_size/365
  if(n_years > 4){
    # If there are records of more than 4 years, take the most recent 4 years
    AMF <- AMF[(nrow(AMF)-365*4*window_size+1):nrow(AMF),]
  }
  
  # Get strat and end time
  t_start <- AMF$TIMESTAMP_START[1]
  t_start <- as.Date(paste(substr(t_start,1,4),substr(t_start,5,6),substr(t_start,7,8),sep="-"))
  t_end <- tail(AMF$TIMESTAMP_START,1)
  t_end <- as.Date(paste(substr(t_end,1,4),substr(t_end,5,6),substr(t_end,7,8),sep="-"))
  time <- seq(from=t_start,to=t_end,by="day")
  
  # Netrad
  netrad <- Var_QC("NETRAD")
  # netrad < 0 equals 0
  netrad[netrad<0] <- 0
  
  # Daily mean net radiation
  NETRAD_daily <- Daily_mean(netrad)
  
  # Output dataframe
  NETRAD_df <- data.frame(
    time = time,
    NETRAD_daily = NETRAD_daily
  )
  
  # Merge with processed data frames ===========
  # Read in processed AMF_df
  AMF_df <- read.csv(paste0(Full_df_path,"AMF/AMF_Full_range_df_",Site_ID,".csv"))
  
  # Merge df with AMF_Netrad
  AMF_df <- AMF_df %>%
    merge(NETRAD_df,by="time") %>%
    # Remove rows with ESI_daily as NA
    filter(!is.na(ESI_daily))
  
  # Calculate ESI following the same normalization used
  # in the threshold extraction
  ESI_rescaled <- rescale(AMF_df$ESI_daily) * 100
  AMF_df$ESI_Z <- ESI_rescaled - median(ESI_rescaled, na.rm = TRUE)
  
  # Get 20th and 80th percentile thresholds
  NETRAD_q20 <- quantile(AMF_df$NETRAD_daily, 0.2, na.rm = TRUE)
  NETRAD_q80 <- quantile(AMF_df$NETRAD_daily, 0.8, na.rm = TRUE)
  
  # Assign low/high net radiation groups
  AMF_df <- AMF_df %>%
    mutate(
      NETRAD_group = case_when(
        NETRAD_daily <= NETRAD_q20 ~ "Low",
        NETRAD_daily >= NETRAD_q80 ~ "High",
        TRUE ~ NA_character_
      )
    )
  
  # Keep only low and high groups
  AMF_compare <- AMF_df %>%
    filter(!is.na(NETRAD_group),
           !is.na(ESI_Z))
  
  # Summary statistics
  ESI_summary <- AMF_compare %>%
    group_by(NETRAD_group) %>%
    summarise(
      n = n(),
      mean_ESI = mean(ESI_Z, na.rm = TRUE),
      median_ESI = median(ESI_Z, na.rm = TRUE),
      sd_ESI = sd(ESI_Z, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Compare ESI between low and high net radiation groups
  t_test_result <- t.test(
    ESI_Z ~ NETRAD_group,
    data = AMF_compare
  )
  
  p_value <- t_test_result$p.value
  
  # Make boxplots of daily ESI between the two net radiation groups
  g <- ggplot(data=AMF_compare,aes(x=NETRAD_group,y=ESI_Z,fill=NETRAD_group))+
    geom_boxplot(outlier.size=0.2)+
    scale_fill_manual(values = c("High" = my_color[2],
                                 "Low" = my_color[3]))+
    my_theme+
    annotate("text",
             x=Inf, y=Inf,
             label=paste0("p = ", format.pval(p_value, digits=2, eps=0.001)),
             hjust=1.1, vjust=1.5)+
    labs(x = "Net radiation group",y="Daily ESI")+
    ggtitle(Site_ID)
  
  # Store the figure
  g_ls[[arrayid]] <- g
  message(arrayid)
}

# Combine the two plots
g_all <- plot_grid(plotlist = g_ls,nrow=1,labels = "auto")




