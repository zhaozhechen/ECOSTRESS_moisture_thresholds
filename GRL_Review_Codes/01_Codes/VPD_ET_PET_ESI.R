# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2026.8.18

# This code is to show the responses of ET, PET, and ESI to VPD
# at the two example AMF sites (US-A32 and US-CF3)

# Reference code: /fs/ess/PAS2204/Code/Validation_ET_ESI_All_AMF/AMF_Full_range_df_20230930.R

# ---- Global --------
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

# Path to raw AMF dataset
AMF_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/00_Data/Raw/AmeriFlux_All_Sites/"
# Path to output figures
Output_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/02_Results/Figures"

# Source plotting functions
source("D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/01_Codes/Plotting_functions.R")

# Sites to test
Site_ID_ls <- c("US-A32","US-CF3")

my_color <- brewer.pal(6,"Set2")
variable_color <- c("ET" = my_color[3],
                    "PET" = my_color[2],
                    "ESI" = my_color[1])

# ---- Functions --------
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
  if(length(var_QC)!=0){
    # Apply QC, only keep QC = 0
    var_QC <- AMF[var_QC][,1]
    var[var_QC!=0] <- NA
  }
  return(var)
}

# Delta
PM_delta <- function(temp){
  tc <- temp - 273.15
  es <- 0.6108*exp(17.27*tc/(tc+237.3))*1000 # Unit Pa
  Delta <- 4098*es/(237.3+tc)^2 # Unit Pa/K
  return(Delta)
}

# ---- Main --------
AMF_daily_all <- c()

for(arrayid in 1:length(Site_ID_ls)){
  Site_ID <- Site_ID_ls[arrayid]
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

  # VPD unit kPa
  vpd <- Var_QC("VPD")/10
  vpd[vpd < 0] <- NA

  # Temperature unit K
  temp <- Var_QC("TA_") + 273.15

  # Delta unit Pa/K
  DEL <- PM_delta(temp) # Unit Pa/K

  # Calculate PET unit mm/day
  rhoa  <- 1.225  # kg/m3
  Cp    <- 1005   # J/kg/K
  gamma <- 66     # Pa/K
  Lv    <- 2453e6 # J/m3
  # Wind speed
  u     <- Var_QC("WS")
  ustar <- Var_QC("USTAR")
  ga    <- 1/(u/ustar^2+6.2*ustar^(-2/3)) # m/s
  # Both example sites have gs = 1/40 m/s based on their PFT
  gs    <- 1/40
  PET   <- (DEL*netrad+rhoa*Cp*ga*vpd*1000)/(DEL+gamma*(1+ga/gs))/Lv*1000*3600*24 # mm/day

  # Calculate ET
  LE <- Var_QC("LE_F")
  ET <- LE*3600*24/(2.45*10^6) # Convert unit from W/m2 to mm/day
  ET[ET<0] <- 0

  # Get daily average variables
  ET_daily <- Daily_mean(ET)
  PET_daily <- Daily_mean(PET)
  ESI_daily <- ET_daily/PET_daily
  # ESI > 1 equals 1
  ESI_daily[ESI_daily>1] <- 1

  # Get daily maximum VPD
  VPD_daily_max <- matrix(data=vpd,nrow=window_size)
  VPD_daily_max <- apply(VPD_daily_max,MARGIN = 2,
                         function(x) if(all(is.na(x))) NA else max(x,na.rm=T))
  VPD_daily_max[VPD_daily_max<0] <- NA

  # Temperature unit C
  # Get daily average T
  temp <- temp - 273.15
  T_daily <- Daily_mean(temp)

  # Daily accumulated precipitation
  P <- Var_QC("P_F")
  # Get daily accumulated P
  P <- matrix(data=P,nrow=window_size)
  P <- colSums(P,na.rm=T)

  AMF_daily <- data.frame(
    Site_ID = Site_ID,
    time = time,
    VPD = VPD_daily_max,
    ET = ET_daily,
    PET = PET_daily,
    ESI = ESI_daily,
    P = P,
    T = T_daily
  )

  # Apply the same filters as the threshold analysis
  AMF_daily <- AMF_daily %>%
    filter(P <= 1,
           T >= 0,
           as.numeric(format(time,"%m")) >= 5,
           as.numeric(format(time,"%m")) <= 9) %>%
    filter(complete.cases(VPD,ET,PET,ESI),
           PET > 0,
           ESI >= 0,
           ESI <= 1)

  AMF_daily_all <- rbind(AMF_daily_all,AMF_daily)
}

# Make one ET, PET, and ESI plot for each site
g_ls <- list()

for(arrayid in 1:length(Site_ID_ls)){
  Site_ID <- Site_ID_ls[arrayid]
  site_df <- AMF_daily_all %>%
    filter(Site_ID==Site_ID_ls[arrayid])

  flux_df <- site_df %>%
    select(time,VPD,ET,PET) %>%
    pivot_longer(cols=c("ET","PET"),
                 names_to="Variable",
                 values_to="Value")

  # Scale ESI to the left axis and recover its original values on the right axis
  flux_axis_max <- ceiling(max(c(site_df$ET,site_df$PET),na.rm=T)/5)*5
  ESI_scale <- flux_axis_max

  g <- ggplot()+
    geom_point(data=flux_df,
               aes(x=VPD,y=Value,color=Variable,fill=Variable),
               size=0.6,alpha=0.2)+
    geom_smooth(data=flux_df,
                aes(x=VPD,y=Value,color=Variable,fill=Variable),
                method="loess",span=0.75,se=TRUE,linewidth=1)+
    geom_point(data=site_df,
               aes(x=VPD,y=ESI*ESI_scale,color="ESI",fill="ESI"),
               size=0.6,alpha=0.2)+
    geom_smooth(data=site_df,
                aes(x=VPD,y=ESI*ESI_scale,color="ESI",fill="ESI"),
                method="loess",span=0.75,se=TRUE,linewidth=1)+
    scale_color_manual(values=variable_color,
                       breaks=c("ET","PET","ESI"))+
    scale_fill_manual(values=variable_color,
                      breaks=c("ET","PET","ESI"))+
    scale_y_continuous(
      limits=c(0,flux_axis_max*1.03),
      oob=scales::oob_keep,
      name=expression("Daily water flux (mm "*d^{-1}*")"),
      sec.axis=sec_axis(as.formula(paste0("~./",ESI_scale)),
                        name="Daily ESI",
                        breaks=seq(0,1,0.25))
    )+
    my_theme+
    theme(legend.position="top",
          plot.title=element_text(hjust=0))+
    labs(x="Daily maximum VPD (kPa)",
         color=NULL,
         fill=NULL)+
    ggtitle(Site_ID)

  g_ls[[arrayid]] <- g
}

# Get one shared legend
g_legend <- get_legend(g_ls[[1]])
g_ls <- lapply(g_ls,function(g) g+theme(legend.position="none"))

# Combine the two site plots
g_panels <- plot_grid(plotlist=g_ls,nrow=1,labels="auto",align="hv")
g_all <- plot_grid(g_legend,g_panels,ncol=1,rel_heights=c(0.13,1))
print_g(g_all,"VPD_ET_PET_ESI",9,4.5)
