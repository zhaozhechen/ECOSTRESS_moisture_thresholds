# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.3.4

# This code is to make RF plots

#########
# Global
#########

library(ggplot2)
library(cowplot)
library(dplyr)

# Input path of intermediate files
Input_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/RF_210m_aggregated/IM_files/"
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
  aspect.ratio = 2/1,
  #legend.key.size = unit(0.3,'cm'),
  legend.title=element_text(size=14)
  #plot.margin = margin(0,0,0,0,'cm')
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

# RF IP ================================================
# Get average and sd IP
IP_sd <- apply(IP_df,1,sd)
IP_mean <- apply(IP_df,1,mean)
IP_df <- data.frame(var_ls = rownames(IP_df),
                    IP_mean,
                    IP_sd,
                    Group = rep(NA,nrow(IP_df)))
# Get groups
IP_df <- IP_df %>% 
  mutate(Group = replace(Group,var_ls == "T_gs_mean"|var_ls == "p_gs_mean"|
                           var_ls == "vpd_gs_mean"|var_ls == "SM_gs_mean","Climate"),
         Group = replace(Group,var_ls == "CH"|var_ls == "RD"|var_ls == "LC","Ecology"),
         Group = replace(Group,var_ls == "Ks"|var_ls == "T_Sand"|var_ls == "n_parameter","Soil"),
         Group = replace(Group,var_ls == "g1"|var_ls == "gpmax"|var_ls == "C"|var_ls == "P50","PHT"),
         Group = replace(Group,var_ls == "LC_shannon"|var_ls == "CH_sd"|var_ls == "DEM_sd","Diversity"),
         Group = replace(Group,var_ls == "DEM","Topography"))
# Get average R2
R2_train <- round(mean(R2$R2_train),2)
R2_test <- round(mean(R2$R2_test),2)

g_IP <- ggplot(IP_df,aes(x=IP_mean,y=reorder(var_ls,IP_mean)))+
  geom_bar(stat="identity",color="black",aes(fill=Group))+
  geom_errorbar(aes(xmin = IP_mean - IP_sd,xmax = IP_mean + IP_sd),width=0.2)+
  my_theme+
  annotate("text",x=0.8*max(IP_df$IP_mean),y=3,label = bquote(r^2~""[train]~"="~.(R2_train)))+
  annotate("text",x=0.8*max(IP_df$IP_mean),y=2,label = bquote(r^2~""[test]~"="~.(R2_test)))+
  labs(y = "",fill="",x="Relative Importance")

print("Complete IP")

# RF PDP =============================================================================
# Initialize a list to store all pdp_df
pdp_ls <- list()

# Loop over each variable 
for(index in 1:length(out[[3]])){
  # Get variable name
  var_name <- names(out[[3]][[index]])[2]
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
  pdp_df <- data.frame(var = pdp_df[,2],
                       Alpha = Alpha_mean,
                       sd = Alpha_sd)
  names(pdp_df) <- c(var_name,"Alpha","sd")
  # Store this pdp_df to pdp_ls
  pdp_ls[[index]] <- pdp_df
  print(index)
}

# Make PDP plots
# Read in full df to calculate data density
df_all <- read.csv(paste0(Input_path,"df_all.csv"))

# Initialize a list to store pdp_g
pdp_g_ls <- list()

for(i in 1:length(pdp_ls)){
  pdp_df <- pdp_ls[[i]]
  var_name <- names(pdp_df)[1]
  # Get this variable in the full df
  var <- df_all[var_name]
  # Calculate data density
  density <- density(var[,1],n=50)$y
  pdp_df$density <- density
  
  g <- ggplot(data=pdp_df,aes(x=.data[[var_name]],y=Alpha))+
    geom_line(aes(y=Alpha,color=density),size=2)+
    geom_ribbon(aes(ymin=Alpha-sd,ymax=Alpha+sd),alpha=0.1)+
    my_theme+
    theme(aspect.ratio = 1/1.5,
          legend.position = "none")+
    scale_color_gradient2(low="grey99",high="#0228FF")
    
  pdp_g_ls[[i]] <- g
}

# Combine all plots
g_all <- plot_grid(plotlist = pdp_g_ls,ncol=4)
g_all <- plot_grid(g_IP,g_all,
                   #rel_heights = c(1,2),
                   rel_widths = c(1,2.5),
                   axis="t")

png(paste0(Output_path,"0.25D_210m_agg_RF_Plots.png"),
    height=12,
    width=16,
    unit='in',
    res=600)
print(g_all)
dev.off()




