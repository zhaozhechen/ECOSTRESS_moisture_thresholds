# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.3.26

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

# RF IP =========================================
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
                           var_ls == "vpd_gs_mean"|var_ls == "SM_gs_mean","Climatic"),
         Group = replace(Group,var_ls == "CH"|var_ls == "RD"|var_ls == "LC","Ecological"),
         Group = replace(Group,var_ls == "Ks"|var_ls == "T_Sand"|var_ls == "n_parameter","Soil hydraulic traits"),
         Group = replace(Group,var_ls == "g1"|var_ls == "gpmax"|var_ls == "C"|var_ls == "P50","Plant hydraulic traits"),
         Group = replace(Group,var_ls == "LC_shannon"|var_ls == "CH_sd"|var_ls == "DEM_sd","Diversity"),
         Group = replace(Group,var_ls == "DEM","Topographic"))
# Rename variables
IP_df <- IP_df %>%
  mutate(var_ls = replace(var_ls,var_ls == "T_gs_mean","Long-term Tair"),
         var_ls = replace(var_ls,var_ls == "T_Sand","Sand fraction"),
         var_ls = replace(var_ls,var_ls == "n_parameter","n"),
         var_ls = replace(var_ls,var_ls == "RD","Root depth"),
         var_ls = replace(var_ls,var_ls == "vpd_gs_mean","Long-term VPD"),
         var_ls = replace(var_ls,var_ls == "p_gs_mean","Long-term P"),
         var_ls = replace(var_ls,var_ls == "CH","Canopy height"),
         var_ls = replace(var_ls,var_ls == "DEM_sd","sd(Elevation)"),
         var_ls = replace(var_ls,var_ls == "LC_shannon","Land cover shannon index"),
         var_ls = replace(var_ls,var_ls == "LC","Land cover"),
         var_ls = replace(var_ls,var_ls == "SM_gs_mean","Long-term SM"),
         var_ls = replace(var_ls,var_ls == "DEM","Elevation")
         )
  
# Get average R2
R2_train <- round(mean(R2$R2_train),2)
R2_test <- round(mean(R2$R2_test),2)
# Make IP plots
g_IP <- ggplot(IP_df,aes(x=IP_mean,y=reorder(var_ls,IP_mean)))+
  geom_bar(stat="identity",color="black",aes(fill=Group))+
  geom_errorbar(aes(xmin = IP_mean - IP_sd,xmax = IP_mean + IP_sd),width=0.2)+
  my_theme+
  scale_fill_brewer(palette = "Set3")+
  annotate("text",x=0.8*max(IP_df$IP_mean),y=3,label = bquote(R^2~""[train]~"="~.(R2_train)))+
  annotate("text",x=0.8*max(IP_df$IP_mean),y=2,label = bquote(R^2~""[test]~"="~.(R2_test)))+
  labs(y = "",fill="Category",x="Relative Importance")

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

# Get the first 8 most important variable names
var_ls <- rownames(IP_df)[order(IP_df$IP_mean,decreasing = TRUE)][1:8]

# Loop over the first 8 most important variables
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
    geom_line(aes(y=Alpha,color=density),size=2)+
    geom_ribbon(aes(ymin=Alpha-sd,ymax=Alpha+sd),alpha=0.2)+
    my_theme+
    theme(aspect.ratio = 1/1.2,
          legend.position = "none")+
    labs(x=x_name,y=bquote(alpha))+
    ylim(c(Alpha_min,Alpha_max))+
    scale_color_gradient2(high="#0228FF")
  pdp_g_ls[[i]] <- g
}

# Combine all plots
g_all <- plot_grid(plotlist=pdp_g_ls,ncol=3)
g_all <- plot_grid(g_IP,g_all,
                   align="h",
                   axis="bt",
                   rel_heights=c(1,1),
                   labels=c("a","b"))

pdf(paste0(Output_path,"0.25D_210m_agg_RF_Plots_V2.pdf"),
    height=6,
    width=14)
print(g_all)
dev.off()




