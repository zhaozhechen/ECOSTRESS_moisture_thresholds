# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.3.4

# This code is to run RF for 100 times
# Response variable is 0.25D alpha aggregated from 210m alpha
# Output include IP df, R2 of RF model, and PDP df for all variables

#########
# Global
#########

library(ggplot2)
library(randomForest)
library(dplyr)
library(cowplot)

# Set the random seed
seed <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID')) 
#seed <- 1

# Input path of 0.25 degree data aggregated from 210m
Input_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/RF_Samples_20240206/0.25D/"
# Output of intermediate files, including IP df, PDP df, R2
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/RF_210m_aggregated/IM_files/"

# Response variable
var_re <- c("Theta")
# List of predicting variables
var_pr <- c("T_gs_mean","Ks","T_Sand","n_parameter",
            "vpd_gs_mean","CH","p_gs_mean","RD","DEM_sd",
            "LC_shannon","LC","SM_gs_mean","g1","DEM","CH_sd",
            "gpmax","C","P50")
var_ls_all <- c(var_re,var_pr)

# Manually fix LC type
LC_type <- c("Shrub_Scrub","Cultivated_Crops","Grassland_Herbaceous",
             "Evergreen_Forest","Deciduous_Forest","Pasture_Hay","Mixed_Forest")[1]

########
# Main
########

# Initialize a list to store all output
out <- list()

# Combine all 0.25 D complete data
df_all <- data.frame()
for(Block_ID in 1:50){
  df <- read.csv(paste0(Input_path,"Var_all_df_Block_",Block_ID,".csv"))
  df_all <- rbind(df_all,df)
  print(Block_ID)
}

# Only keep required variables
df_all <- df_all[var_ls_all]
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

# Run RF  =====================================================
# Split dataset into Training set 70% and Testing set 30%
set.seed(seed)
train_index <- sample(1:nrow(df_all),size=nrow(df_all)*0.7)
df_train  <- df_all[train_index,]
df_test   <- df_all[-train_index,]

# Save df_all to calculate data density
#write.csv(df_all,paste0(Output_path,"df_all.csv"))

print('Fitting RF model ...')
# # of trees: 500
# Get formula
set.seed(seed)
f <- as.formula(paste(var_re,"~.",sep=""))
system.time(rf <- randomForest(f,
                               data = df_train,
                               importance=TRUE,
                               type = "regression",
                               ntree=500))
print("Complete RF model")

# Assess RF performance =======================================
# Record r for testing and training set
r_train <- cor(rf$predicted,df_train[var_re],use="pairwise.complete.obs")
rf.test.pred <- predict(rf,df_test)
r_test  <- cor(rf.test.pred,df_test[var_re],use="pairwise.complete.obs")
# Get relative importance
rf_importance <- as.data.frame(importance(rf))
rf_importance$Feature <- rownames(rf_importance)
rf_importance <- rf_importance %>% rename("Importance" = "%IncMSE")
# Rescale importance to sum = 1
rf_importance$Importance <- rf_importance$Importance/sum(rf_importance$Importance)

# Save this relative importance df
out[[1]] <- rf_importance

# Output R2
R2 <- data.frame(R2_train = r_train^2,
                 R2_test = r_test^2)
names(R2) <- c("R2_train","R2_test")
out[[2]] <- R2

# Approximate Response Curves =========================================
# Reference: https://cran.r-project.org/web/packages/pdp/vignettes/pdp-approximate.pdf
# Change one predictor, fix all others at typical value - median, fix LC at the most frequent value

# List of variables to test
var_ls <- var_pr[var_pr != "LC"]
# List to store all pdp_df
pdp_ls <- list()

# Determine which var to change
for(index in 1:length(var_ls)){
  # Name of the changing variable
  var_chg_name <- var_ls[index]
  # Create a sequence of values for this changing variable
  var_chg <- seq(min(df_all[var_chg_name]),
                 max(df_all[var_chg_name]),
                 length.out = 50)
  # Name of the fixing variables
  var_fix_name <- var_pr[var_pr!=var_chg_name & var_pr!="LC"]
  # Create a dummy training set, fix all fixing variables at median
  train_dummy <- sapply(df_train[var_fix_name],median)
  # Repeat 50 times
  train_dummy <- as.data.frame(lapply(train_dummy,rep,50))
  # Fix LC
  LC <- rep(as.factor(LC_type),50)
  # A trick to equalize classes of training set for LC levels
  LC <- c(df_train$LC[1],LC)
  LC <- LC[-1]
  train_dummy$LC <- LC
  # Add the changing variable
  train_dummy$var <- var_chg
  colnames(train_dummy)[colnames(train_dummy)=="var"] <- var_chg_name
  # Use this dummy predictor set to predict theta
  Theta_dummy <- predict(rf,train_dummy)
  pdp_df <- data.frame(Alpha = Theta_dummy,
                       var_chg)
  names(pdp_df) <- c("Alpha",var_chg_name)
  pdp_ls[[index]] <- pdp_df
  print(index)
}

# Store all pdp_df
out[[3]] <- pdp_ls

# Output
saveRDS(out,paste0(Output_path,"RF_IM_",seed,".rds"))
