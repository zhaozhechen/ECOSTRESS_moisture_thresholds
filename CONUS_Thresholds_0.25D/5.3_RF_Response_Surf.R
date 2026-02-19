# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.3.5

# This code is to get response surface of Theta, in response to the combinations of 
# Diversity variables and other variables

# Reference code: /fs/ess/PAS2204/Code/CONUS_Threshold/CONUS_RF_Final_3/RF_0.25D_Analysis_20240219.R

#########
# Global
#########

library(randomForest)
library(ggplot2)
library(cowplot)
library(dplyr)

# Input path for 0.25D aggregated from 210m data
Input_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/RF_Samples_20240206/0.25D/"
# Output path for response surface
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Figures/"

# Response variable
var_re <- "Theta"
# List of predicting variables
var_pr <- c("T_gs_mean","p_gs_mean","vpd_gs_mean","SM_gs_mean",
            "CH","CH_sd","RD","LC","LC_shannon","DEM_sd","DEM",
            "Ks","n_parameter","T_Sand",
            "g1","P50","gpmax","C")
var_ls <- c(var_re,var_pr)

# Manually fix LC type
LC_type <- c("Shrub_Scrub","Cultivated_Crops","Grassland_Herbaceous",
             "Evergreen_Forest","Deciduous_Forest","Pasture_Hay","Mixed_Forest")[1]

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

# This function is to plot correlation matrix
plot_CM <- function(Corr_matrix){
  CM <- ggcorrplot::ggcorrplot(Corr_matrix,type="lower",lab=T,
                               colors=c("#619CFF","white","#F8766D"),
                               lab_size = 4,
                               tl.cex = 30)+
    my_theme+
    theme(aspect.ratio = 1/1)
  return(CM)
}

# This function is to output plots
# If type = "IP",this is importance plot
# If type = "PDP",this is PDP plot
# n is number of plots, defalut is 16
plotp <- function(g,type,title,n=16){
  if(type == "IP"){
    x <- 5
    y <- 6
  }else if(type == "PDP"){
    x <- ifelse(n>=4,24,n*6)
    y <- 3*ceiling(n/4)
  }
  png(paste0(Output_path,title,".png"),
      width=x,
      height=y,
      unit='in',
      res=600)
  print(g)
  dev.off()
}

#######
# Main
#######

# Combine all 0.25 D complete data
df_all <- data.frame()
for(Block_ID in 1:50){
  df <- read.csv(paste0(Input_path,"Var_all_df_Block_",Block_ID,".csv"))
  df_all <- rbind(df_all,df)
  print(Block_ID)
}
# Only keep required variables
df_all <- df_all[var_ls]
# Only keep complete observations
df_all <- na.omit(df_all)

# Convert LC to characters
df_all$LC[df_all$LC == 41] <- "Deciduous_Forest"
df_all$LC[df_all$LC == 42] <- "Evergreen_Forest"
df_all$LC[df_all$LC == 43] <- "Mixed_Forest"
df_all$LC[df_all$LC == 51] <- "Dwarf_scrub"
df_all$LC[df_all$LC == 52] <- "Shrub_Scrub"
df_all$LC[df_all$LC == 71] <- "Grassland_Herbaceous"
df_all$LC[df_all$LC == 72] <- "Sedge_Herbaceous"
df_all$LC[df_all$LC == 81] <- "Pasture_Hay"
df_all$LC[df_all$LC == 82] <- "Cultivated_Crops"
df_all$LC <- as.factor(df_all$LC)

# Get RF model ========================================
# Split dataset into Training set 70% and Testing set 30%
set.seed(1)
train_index <- sample(1:nrow(df_all),size=nrow(df_all)*0.7)
df_train  <- df_all[train_index,]
df_test   <- df_all[-train_index,]

print('Fitting RF model ...')
# # of trees: 500
# Get formula
set.seed(1)
f <- as.formula(paste(var_re,"~.",sep=""))
system.time(rf <- randomForest(f,
                               data = df_train,
                               importance=TRUE,
                               type = "regression",
                               ntree=500))
print("Complete RF model")

# Response surface ======================================
# Change two variables at a time, fix others at median, fix LC at the most frequent value
# The first variable is one of diversity variables
var1_ls <- c("LC_shannon","DEM_sd","CH_sd")
# The first variable is all other variables
var2_ls <- setdiff(var_pr,c(var1_ls,"LC"))
# Initialize a list to store all plots
g_surf_ls <- list()
# Loop over the first variable
for(i in 1:length(var1_ls)){
  var_chg_name_1 <- var1_ls[i]
  for(j in 1:length(var2_ls)){
    var_chg_name_2 <- var2_ls[j]
    var_chg_1 <- seq(min(df_train[var_chg_name_1]),
                     max(df_train[var_chg_name_1]),
                     length.out = 50)
    var_chg_2 <- seq(min(df_train[var_chg_name_2]),
                     max(df_train[var_chg_name_2]),
                     length.out = 50)
    # Names of the fixing variables
    var_fix_name <- setdiff(var_pr,c(var_chg_name_1,var_chg_name_2,"LC"))
    # Create a dummy training set, fix all fixing variables at median
    train_dummy <- sapply(df_train[var_fix_name],median)
    # Repeat 2500 times
    train_dummy <- as.data.frame(lapply(train_dummy,rep,2500))
    # fix LC at the most frequent value
    LC <- rep(as.factor(LC_type),2500)
    # A trick to equalize classes of training set for LC levels
    LC <- c(df_train$LC[1],LC)
    LC <- LC[-1]
    train_dummy$LC <- LC
    # Get all combinations of values of the two changing variables
    var_chg_1_2 <- expand.grid(var_chg_1,var_chg_2)
    names(var_chg_1_2) <- c(var_chg_name_1,var_chg_name_2)
    # Add the changing variables
    train_dummy <- cbind(train_dummy,var_chg_1_2)
    # Use this dummy predictor sets to predict theta
    Theta_dummy <- predict(rf,train_dummy)
    
    # Get a data frame of the two changing variables, and predicted Theta
    surf_df <- data.frame(Alpha = Theta_dummy,
                          var_chg_1_2)
    # Make plot of response surface, color coded by Alpha
    g_surf <- ggplot(data=surf_df,aes(x = .data[[var_chg_name_1]],
                                      y = .data[[var_chg_name_2]],
                                      fill = Alpha))+
      geom_tile()+
      scale_fill_distiller(type="div",palette = 'RdBu')+
      my_theme+
      theme(aspect.ratio = 1/1)

    g_surf_ls[[(i-1)*length(var2_ls)+j]] <- g_surf
  }
  print(i)
}

# Combine each of the three sets of plots
for(i in 1:3){
  var_name <- var1_ls[i]
  g <- plot_grid(plotlist=g_surf_ls[c(((i-1)*length(var2_ls)+1):(i*length(var2_ls)))],ncol=4)
  plotp(g,'PDP',paste0("0.25D_RF_Surf_",var_name))
  print(i)
}





