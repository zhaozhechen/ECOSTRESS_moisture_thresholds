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

# Path to raw AMF dataset
AMF_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/00_Data/Raw/AmeriFlux_All_Sites/"
# Path to full-range df of AMF and RS
Full_df_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/00_Data/Processed/Full_range_df/"
# Path to output figures
Output_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/02_Results/Figures"

# Source plotting functions
source("D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/01_Codes/Plotting_functions.R")
source("D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/01_Codes/General_functions.R")

# Sites to test
Site_ID_ls <- c("US-A32","US-CF3")

my_color <- brewer.pal(6,"Set2")

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
  # Change missing values to NA
  netrad[netrad==-9999] <- NA
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
  
  # Get 25th and 75th percentile thresholds for each month
  # Assign low/high net radiation groups based on monthly thresholds
  AMF_df <- AMF_df %>%
    mutate(month = format(as.Date(time), "%m"),
           year_month = format(as.Date(time), "%Y-%m")) %>%
    group_by(year_month) %>%
    mutate(
      NETRAD_q25 = quantile(NETRAD_daily, 0.25, na.rm = TRUE),
      NETRAD_q75 = quantile(NETRAD_daily, 0.75, na.rm = TRUE),
      NETRAD_group = case_when(
        NETRAD_daily <= NETRAD_q25 ~ "Low",
        NETRAD_daily >= NETRAD_q75 ~ "High",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup()
  
  # Keep only low and high groups
  AMF_compare <- AMF_df %>%
    filter(!is.na(NETRAD_group),
           !is.na(ESI_daily))
  AMF_compare$NETRAD_group <- factor(AMF_compare$NETRAD_group,
                                     levels=c("Low","High"))
  
  # Calculate monthly mean ESI for each net radiation group
  ESI_monthly_low <- AMF_compare %>%
    filter(NETRAD_group=="Low") %>%
    group_by(year_month) %>%
    summarise(
      ESI_low = mean(ESI_daily, na.rm = TRUE),
      .groups = "drop"
    )
  ESI_monthly_high <- AMF_compare %>%
    filter(NETRAD_group=="High") %>%
    group_by(year_month) %>%
    summarise(
      ESI_high = mean(ESI_daily, na.rm = TRUE),
      .groups = "drop"
    )
  ESI_monthly <- merge(ESI_monthly_low,ESI_monthly_high,by="year_month")
  
  # Compare monthly mean ESI between low and high net radiation groups
  t_test_result <- t.test(ESI_monthly$ESI_high,
                          ESI_monthly$ESI_low,
                          paired=TRUE)
  
  p_value <- t_test_result$p.value
  
  # Make boxplots of daily ESI between the two net radiation groups
  g <- ggplot(data=AMF_compare,aes(x=NETRAD_group,y=ESI_daily,fill=NETRAD_group))+
    geom_boxplot(outlier.size=0.2)+
    scale_fill_manual(values = c("High" = my_color[2],
                                 "Low" = my_color[3]))+
    my_theme+
    annotate("text",
             x=Inf, y=Inf,
             label=paste0("p = ",format.pval(p_value, digits=2, eps=0.001)),
             hjust=1.1, vjust=1.5)+
    labs(x = "Net radiation group",y="Daily ESI")+
    ggtitle(Site_ID)
  
  # Store the figure
  g_ls[[arrayid]] <- g
  message(arrayid)
}

# Combine the two plots
g_all <- plot_grid(plotlist = g_ls,nrow=1,labels = "auto")
print_g(g_all,"Clear_Sky_comparison",8,4)


