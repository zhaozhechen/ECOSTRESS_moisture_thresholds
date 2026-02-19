# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2025.2.13

# This code is to re-plot alpha PDP conditional on precipitation

##########
# Global
##########

library(ggplot2)
library(randomForest)
library(dplyr)
library(cowplot)
library(grid)
library(gridExtra)

# Input path for 0.25 degree data
Input_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/RF_Samples_20240206/0.25D/"
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Figures/"

# Response variable
var_re <- c("Theta")
# List of 17 predicting variables
var_pr <- c("T_gs_mean","T_Sand","Ks","n_parameter","RD","vpd_gs_mean","p_gs_mean",
            "CH","DEM_sd","LC_shannon","g1","LC","SM_gs_mean","DEM","gpmax","C","P50")
var_ls <- c(var_re,var_pr)

# Manually fix LC type
LC_type <- c("Shrub_Scrub","Cultivated_Crops","Grassland_Herbaceous",
             "Evergreen_Forest","Deciduous_Forest","Pasture_Hay","Mixed_Forest")[4]

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
  legend.title=element_text(size=14),
  axis.title = element_text(size=14)
  #plot.margin = margin(0,0,0,0,'cm'),
)

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

# Get conditional PDP ========================================
# Get PDP, changing one target variable, and fixing a second variable at five levels, fixing others at median
# Fixing a second variable at min,25th,median,75th,max
# Also get the average response curves of these 5 levels

# The first variable that need to be changing
var_chg_name_1 <- "T_gs_mean"
var_chg_1 <- seq(min(df_train[var_chg_name_1]),
                 max(df_train[var_chg_name_1]),
                 length.out = 50)
# The second variable that need to be fixed at 5 levels
var_chg_name_2 <- "p_gs_mean"
# The five levels
var_chg_2_lv <- sapply(df_train[var_chg_name_2], quantile)
# Names of the fixing variables
var_fix_name <- var_pr[var_pr != "LC" & var_pr != var_chg_name_1 & var_pr != var_chg_name_2]


# Initialize a vector to store theta predicted at the five levels
Theta_all <- c()
# Initialize a vector to store mean theta at the five levels
Theta_mean <- c()

# Loop over the five levels
for(i in 1:5){
  # Fix the second variable at the level i
  var_chg_2 <- rep(var_chg_2_lv[i],50)
  
  # Create a dummy training set, fix all fixing variables at median
  train_dummy <- sapply(df_train[var_fix_name],median)
  # Repeat 50 times
  train_dummy <- as.data.frame(lapply(train_dummy,rep,50))
  # fix LC at the most frequent value
  LC <- rep(as.factor(LC_type),50)
  # A trick to equalize classes of training set for LC levels
  LC <- c(df_train$LC[1],LC)
  LC <- LC[-1]
  train_dummy$LC <- LC
  # Combine the first changing variable
  train_dummy$var <- var_chg_1
  colnames(train_dummy)[colnames(train_dummy)=="var"] <- var_chg_name_1
  # Add this second variable to the train_dummy set
  train_dummy$var <- var_chg_2
  colnames(train_dummy)[colnames(train_dummy)=="var"] <- var_chg_name_2
  # Use this dummy predictor sets to predict theta
  Theta_dummy <- predict(rf,train_dummy)
  Theta_all <- c(Theta_all,Theta_dummy)
  # Get mean Theta
  Theta_mean <- cbind(Theta_mean,Theta_dummy)
}

# Make df for plots
Theta_df <- data.frame(Theta = Theta_all,
                       var = rep(var_chg_1,5),
                       Type = rep(c("min","25th","median","75th","max"),each=50))
Theta_df$Type <- factor(Theta_df$Type,levels= c("min","25th","median","75th","max"))

# Make conditional PDP
g_PDP <- ggplot(Theta_df,aes(x=var,y=Theta,color=Type))+
  geom_line()+
  my_theme+
  labs(x=bquote("T"[air]~""~(degree~C)),y=bquote(alpha~(degree)),
       color="Precipitation")+
  scale_color_manual(values=c(min = "#2E64FE",
                              `25th` = "#A9BCF5",
                              median = "#A9F5BC",
                              `75th` = "#FFBF00",
                              max = "#FF0000"))

pdf(paste0(Output_path,"0.25D_210m_agg_RF_PDP_T&P.pdf"),
    height = 3,width=5)
print(g_PDP)
dev.off()


