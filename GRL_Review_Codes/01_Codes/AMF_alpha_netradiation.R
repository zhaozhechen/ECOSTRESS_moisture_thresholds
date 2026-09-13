# Date: 2026.9.13
# Compare AMF moisture sensitivity using all observations and excluding low net radiation.
# Reference: Validation_Final_2/Validation_thresholds_20250202.R
# Radiation processing follows Clear_Sky_comparison.R.

# ---- Global --------
library(raster)
library(dplyr)
library(scales)
library(ggplot2)
library(cowplot)

args <- commandArgs(trailingOnly=FALSE)
script_file <- sub("^--file=", "", args[grepl("^--file=", args)])
if(length(script_file)==1){
  Review_path <- dirname(dirname(normalizePath(script_file, winslash="/")))
}else{
  Review_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes"
}
source(file.path(Review_path,"01_Codes/AMF_alpha_functions.R"))
source(file.path(Review_path,"01_Codes/General_functions.R"))
AMF_path <- file.path(Review_path,"00_Data/Raw/AmeriFlux_All_Sites")
Input_path <- file.path(Review_path,"00_Data/Processed/Full_range_df/AMF")
Output_path <- file.path(Review_path,"02_Results/Figures")
Table_path <- file.path(Review_path,"02_Results/Tables")
Audit_path <- file.path(Review_path,"02_Results/Intermediate/AMF_alpha_netradiation")
for(path in c(Output_path,Table_path,Audit_path)) dir.create(path,recursive=TRUE,showWarnings=FALSE)

n_grid <- 7
n_k <- 3
sigma <- 1
n_SLR <- 4
Slice_R2 <- 0.8
n_Fit0 <- 3
num_iteration <- 100
options(lifecycle_verbosity="quiet")
Site_ls <- c("US-A32","US-CF3")
summary_ls <- list()
cutoff_ls <- list()
boundary_ls <- list()
provenance_ls <- list()

# ---- Main --------
for(Site_ID in Site_ls){
  folder_name <- list.dirs(AMF_path,recursive=FALSE,full.names=TRUE)
  folder_name <- folder_name[grepl(Site_ID,basename(folder_name),fixed=TRUE)]
  stopifnot(length(folder_name)==1)
  file_name <- list.files(folder_name,pattern="SUBSET_HH_.*csv$",full.names=TRUE)
  window_size <- 48
  if(length(file_name)==0){
    file_name <- list.files(folder_name,pattern="SUBSET_HR_.*csv$",full.names=TRUE)
    window_size <- 24
  }
  stopifnot(length(file_name)==1)
  AMF <- read.csv(file_name)
  if(nrow(AMF)/window_size/365 > 4){
    AMF <- tail(AMF,365*4*window_size)
  }
  # Match the previous radiation analysis: QC=0, missing to NA, negatives to zero.
  netrad <- Var_QC("NETRAD")
  netrad[netrad==-9999] <- NA
  netrad[netrad<0] <- 0
  dates <- as.Date(substr(as.character(AMF$TIMESTAMP_START),1,8),format="%Y%m%d")
  stopifnot(nrow(AMF) %% window_size==0,
            all(vapply(split(dates,ceiling(seq_along(dates)/window_size)),
                       function(x) length(unique(x))==1,logical(1))))
  NETRAD_df <- data.frame(time=dates[seq(1,length(dates),by=window_size)],
                          NETRAD_daily=Daily_mean(netrad))
  stopifnot(!anyDuplicated(NETRAD_df$time))

  input_file <- file.path(Input_path,paste0("AMF_Full_range_df_",Site_ID,".csv"))
  full_df <- read.csv(input_file)
  names(full_df)[names(full_df)=="VPD_daily_max"] <- "Daily_max_VPD"
  full_df$time <- as.Date(full_df$time)
  stopifnot(!anyDuplicated(full_df$time))
  full_df <- left_join(full_df,NETRAD_df,by="time")
  full_df$year_month <- format(full_df$time,"%Y-%m")

  # Same year-month cutoff population as Clear_Sky_comparison.R: nonmissing ESI.
  cutoffs <- full_df %>% filter(!is.na(ESI_daily)) %>%
    group_by(year_month) %>%
    summarise(NETRAD_q25=quantile(NETRAD_daily,0.25,na.rm=TRUE),
              N_ESI=n(),N_radiation=sum(is.finite(NETRAD_daily)),.groups="drop")
  full_df <- left_join(full_df,cutoffs[,c("year_month","NETRAD_q25")],by="year_month")
  full_df$Keep_radiation <- is.finite(full_df$NETRAD_daily) &
    is.finite(full_df$NETRAD_q25) & full_df$NETRAD_daily >= full_df$NETRAD_q25
  full_df$Complete_alpha <- complete.cases(full_df[,c("ESI_daily","Daily_max_VPD","SM_daily")])
  write.csv(full_df,file.path(Audit_path,paste0(Site_ID,"_daily_selection.csv")),row.names=FALSE)
  cutoff_ls[[Site_ID]] <- mutate(cutoffs,Site_ID=Site_ID)
  datasets <- list(All=full_df,NETRAD_ge_Q25=full_df[full_df$Keep_radiation,])

  for(group in names(datasets)){
    data_df <- datasets[[group]]
    # Apply the original calculations, including normalization within each dataset.
    features_all <- SLR_features(data_df)
    features <- Filter_SLR_feature(features_all)
    result <- SM_VPD(data_df)
    normalized <- Normalize_df(data_df)
    # Save all slice diagnostics and retained zero-crossing coordinates.
    features$Fit0 <- -features$Intercept/features$Slope
    fit0 <- data.frame(SM=c(features$Fit0[features$Fix=="VPD"],features$Level[features$Fix=="SM"]),
                       VPD=c(features$Level[features$Fix=="VPD"],features$Fit0[features$Fix=="SM"]))
    fit0 <- fit0[complete.cases(fit0) & fit0$SM>=min(normalized$SM) &
                   fit0$SM<=max(normalized$SM) & fit0$VPD>=min(normalized$VPD) &
                   fit0$VPD<=max(normalized$VPD),]
    write.csv(features_all,file.path(Audit_path,paste0(Site_ID,"_",group,"_slices.csv")),row.names=FALSE)
    write.csv(fit0,file.path(Audit_path,paste0(Site_ID,"_",group,"_boundary_points.csv")),row.names=FALSE)
    summary_ls[[paste(Site_ID,group)]] <- cbind(
      data.frame(Site_ID=Site_ID,Dataset=group,N_rows=nrow(data_df),N_complete=nrow(normalized),
                 N_complete_missing_radiation=sum(data_df$Complete_alpha & !is.finite(data_df$NETRAD_daily)),
                 N_retained_slices=nrow(features),N_boundary_points=nrow(fit0),
                 First_date=min(data_df$time[data_df$Complete_alpha]),
                 Last_date=max(data_df$time[data_df$Complete_alpha]),
                 SM_min=min(data_df$SM_daily,na.rm=TRUE),SM_max=max(data_df$SM_daily,na.rm=TRUE),
                 VPD_min=min(data_df$Daily_max_VPD,na.rm=TRUE),VPD_max=max(data_df$Daily_max_VPD,na.rm=TRUE),
                 ESI_median=median(data_df$ESI_daily,na.rm=TRUE)),result)
    boundary_ls[[paste(Site_ID,group)]] <- mutate(fit0,Site_ID=Site_ID,Dataset=group)
  }
  provenance_ls[[Site_ID]] <- data.frame(Site_ID=Site_ID,File=c(file_name,input_file),
                                         MD5=unname(tools::md5sum(c(file_name,input_file))))
  message("Completed ",Site_ID)
}

summary_df <- bind_rows(summary_ls)
# Reproduce the saved manuscript-example slopes (rounded in the original script).
reference_slopes <- c("US-A32"=1.2897607,"US-CF3"=8.8470461)
baseline <- summary_df[summary_df$Dataset=="All",]
stopifnot(all(abs(baseline$Slope-reference_slopes[baseline$Site_ID]) < 1e-6),
          all(summary_df$N_complete_missing_radiation==0))
names(summary_df)[names(summary_df)=="theta"] <- "Alpha_degrees"
names(summary_df)[names(summary_df)=="theta_sd"] <- "Alpha_sd_degrees"
summary_df <- summary_df %>% group_by(Site_ID) %>%
  mutate(Alpha_change_degrees=Alpha_degrees-Alpha_degrees[Dataset=="All"],
         Retained_percent=100*N_complete/N_complete[Dataset=="All"]) %>% ungroup()
write.csv(summary_df,file.path(Table_path,"AMF_alpha_netradiation_summary.csv"),row.names=FALSE)
write.csv(bind_rows(cutoff_ls),file.path(Table_path,"AMF_alpha_netradiation_monthly_cutoffs.csv"),row.names=FALSE)
write.csv(bind_rows(provenance_ls),file.path(Audit_path,"input_provenance.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(Audit_path,"sessionInfo.txt"))

# Direct comparison of alpha, using the previous radiation-plot colors.
summary_df$Dataset <- factor(summary_df$Dataset,levels=c("All","NETRAD_ge_Q25"))
g_ls <- lapply(Site_ls,function(site){
  df <- filter(summary_df,Site_ID==site)
  ggplot(df,aes(Dataset,Alpha_degrees,fill=Dataset))+
    geom_col(width=0.6,color="black")+
    geom_text(aes(label=sprintf("%.1f",Alpha_degrees)),vjust=-0.5,size=4.5)+
    scale_fill_manual(values=c("All"="#8DA0CB","NETRAD_ge_Q25"="#FC8D62"))+
    scale_x_discrete(labels=c("All","Excluding low\nradiation"))+
    scale_y_continuous(limits=c(0,100),breaks=seq(0,90,30),expand=expansion(mult=c(0,0)))+
    labs(x=NULL,y=expression(alpha~(degree)),title=site)+
    theme_classic(base_size=14)+theme(panel.border=element_rect(color="black",fill=NA),
                                    legend.position="none",plot.title=element_text(hjust=0))
})
g <- plot_grid(plotlist=g_ls,nrow=1,labels=c("a","b"))
ggsave(file.path(Output_path,"AMF_alpha_netradiation.png"),g,width=8,height=4,dpi=600,bg="white")
ggsave(file.path(Output_path,"AMF_alpha_netradiation.pdf"),g,width=8,height=4,bg="white")
print(as.data.frame(summary_df[,c("Site_ID","Dataset","N_complete","Alpha_degrees","Alpha_sd_degrees","Alpha_change_degrees","R2")]))
