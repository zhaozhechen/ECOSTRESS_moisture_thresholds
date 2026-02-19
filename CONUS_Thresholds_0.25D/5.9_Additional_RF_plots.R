# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2025.6.23

# This code is to plot additional response curves from random forest results

#########
# Global
#########

library(ggplot2)
library(cowplot)
library(dplyr)

# Input path of intermediate files
Input_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/RF_210m_aggregated/IM_files_V2/"
# Output path of RF figures
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Figures/"

############
# Functions
############
my_theme <- theme(
  axis.line=element_line(color="black"),
  panel.background = element_blank(),
  text = element_text(size=14),
  panel.border = element_rect(colour="black",fill=NA),
  legend.key = element_blank(),
  #aspect.ratio = 2/1,
  #legend.key.size = unit(0.3,'cm'),
  legend.title=element_text(size=12),
  axis.title = element_text(size=12)
  #plot.margin = margin(0,0,0,0,'cm'),
)

#######
# Main
#######

# Initialize IP_df
IP_df <- c()
# Initialize R2
R2 <- c()

for(i in 1:100){
  # Read in the intermediate file
  out <- readRDS(paste0(Input_path,"RF_IM_",i,".rds"))
  
  # Combine all IP df
  IP_df_tmp <- out[[1]][order(out[[1]]$Feature),]
  IP_df <- cbind(IP_df,IP_df_tmp$Importance)
  rownames(IP_df) <- rownames(IP_df_tmp)
  
  # Combine all R2
  R2 <- rbind(R2,out[[2]])
  print(i)
}

# Get average and sd IP
IP_sd <- apply(IP_df,1,sd)
IP_mean <- apply(IP_df,1,mean)
IP_df <- data.frame(var_ls = rownames(IP_df),
                    IP_mean,
                    IP_sd,
                    Group = rep(NA,nrow(IP_df)))

# RF PDP =============================================================================
# Initialize a list to store all pdp_df
pdp_ls <- list()

# All variable names in the pdp_df input
var_pr <- c("T_gs_mean","Ks","T_Sand","n_parameter",
            "vpd_gs_mean","CH","p_gs_mean","RD","DEM_sd",
            "LC_shannon","LC","SM_gs_mean","g1","DEM",
            "gpmax","C","P50")
var_ls_all <- var_pr[var_pr != "LC"]

# Get mean and sd of 100 runs for each variable pdp
# Loop over each variable
for(index in 1:length(var_ls_all)){
  # Get variable name
  var_name <- var_ls_all[index]
  pdp_df_all <- c()
  # Loop over the 100 times
  for(i in 1:100){
    # Read in the pdp_df
    out <- readRDS(paste0(Input_path,"RF_IM_",i,".rds"))
    pdp_df <- out[[3]][[index]]
    pdp_df_all <- cbind(pdp_df_all,pdp_df$Alpha)
  }
  # Get mean and sd of Alpha
  Alpha_mean <- apply(pdp_df_all,1,mean)
  Alpha_sd <- apply(pdp_df_all,1,sd)
  # Make df
  pdp_df <- data.frame(var = pdp_df$x,
                       Alpha = Alpha_mean,
                       sd = Alpha_sd)
  names(pdp_df)[1] <- var_name
  # Store this pdp_df to pdp_ls
  pdp_ls[[index]] <- pdp_df
  print(index)
}

# Read in full df to calculate data density
df_all <- read.csv(paste0(Input_path,"df_all.csv"))

# Get the 8 least important variable names
var_ls <- rownames(IP_df)[order(IP_df$IP_mean,decreasing = TRUE)][9:17]
var_ls <- var_ls[var_ls != "LC"]

# Loop over the least 8 most important variables
# Get range of Alpha
# Initialize a vector to store all max Alpha
Alpha_max <- c()
Alpha_min <- c()
for(i in 1:8){
  # variable name
  var_name <- var_ls[i]
  # Get index of the pdp_df in the pdp_df list
  pdp_index <- which(var_ls_all == var_name)
  # Get pdp_df
  pdp_df <- pdp_ls[[pdp_index]]
  # Max alpha
  Alpha_max <- c(Alpha_max,max(pdp_df$Alpha+pdp_df$sd))
  Alpha_min <- c(Alpha_min,min(pdp_df$Alpha-pdp_df$sd))
}
# Get max alpha
Alpha_max <- max(Alpha_max)
Alpha_min <- min(Alpha_min)

# Make PDP plots
# Initialize a list to store all pdp plots
pdp_g_ls <- list()
for(i in 1:8){
  # variable name
  var_name <- var_ls[i]
  # Get index of the pdp_df in the pdp_df list
  pdp_index <- which(var_ls_all == var_name)
  # Get pdp_df
  pdp_df <- pdp_ls[[pdp_index]]
  # Get this variable in the full df
  var <- df_all[var_name]
  # Calculate data density
  pdp_df$density <- density(var[,1],n=nrow(pdp_df))$y
  # x axis name
  x_name <- IP_df$var_ls[rownames(IP_df)==var_name]
  
  g <- ggplot(data=pdp_df,aes(x=.data[[var_name]],y=Alpha))+
    geom_line(aes(y=Alpha,color=density),size=0.8)+
    geom_ribbon(aes(ymin=Alpha-sd,ymax=Alpha+sd),alpha=0.2)+
    my_theme+
    theme(aspect.ratio = 1/1.2,
          legend.position = "none")+
    labs(x=x_name,y=bquote(alpha~ (degree)))+
    ylim(c(Alpha_min,Alpha_max))+
    scale_color_gradient2(high="#0228FF")
  pdp_g_ls[[i]] <- g
}

# Combine all plots
g_all <- plot_grid(plotlist=pdp_g_ls,nrow=2)

pdf(paste0(Output_path,"0.25D_210m_agg_additional_RF_Plots.pdf"),
    height=4,
    width=9.3)
print(g_all)
dev.off()

