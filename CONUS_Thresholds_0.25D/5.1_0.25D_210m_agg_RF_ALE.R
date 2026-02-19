# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2025.10.13

# This code is to run RF for 100 times
# Response variable is 0.25D alpha aggregated from 210m alpha


# Output include IP df, R2 of RF model, and ALE df for all variables
# Note: The output ALE are plotted using iml package
# Predictors have no CH_sd or aridity index (AI), keep CH
# Predictors remove plant hydraulic traits

#########
# Global
#########
#dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
#dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)
library(ggplot2)
library(randomForest)
library(dplyr)
library(cowplot)
library(iml)

# Set the random seed
seed <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID')) 
#seed <- 1

# Input path of 0.25 degree data aggregated from 210m full data frame
df_all <- read.csv("/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/RF_210m_aggregated/IM_files_V3/df_all.csv")
# Output of intermediate files, including IP df, PDP df, R2
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/RF_210m_aggregated/IM_files_ALE/"

# Response variable
var_re <- c("Theta")
# List of predicting variables
var_pr <- c("T_gs_mean","Ks","T_Sand","n_parameter",
            "vpd_gs_mean","CH","p_gs_mean","RD","DEM_sd",
            "LC_shannon","LC","SM_gs_mean","DEM")
var_ls_all <- c(var_re,var_pr)

########
# Main
########
df_all <- df_all %>%
  select(-X)

df_all$LC <- as.factor(df_all$LC)

# Run RF  =====================================================
# Initiliaze a list to store all output
out <- list()

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

# ALE ==============================================
# List of variables to test
var_ls <- var_pr[var_pr != "LC"]
# List to store all ale_df
ale_ls <- list()

# Create a predictor object, wrap the random forest model
predictor <- Predictor$new(
  model = rf,
  data = df_all[var_pr],
  y = df_all[[var_re]]
)

# Calculate mean prediction across all data
mean_pred <- mean(predict(rf, df_all))

for(index in 1:length(var_ls)){
  # Name of the target var (for which the ALE should be plotted)
  var_name <- var_ls[index]
  
  # Calculate ALE
  ale_obj <- FeatureEffect$new(predictor,feature = var_name,method="ale",grid.size = 50)
  
  # Add the mean prediction back to shift ALE to theta scale
  ale_df <- ale_obj$results
  ale_df$Alpha <- ale_df$.value + mean_pred
  
  # Make output competible for plotting
  ale_df <- data.frame(x = ale_df[[var_name]],
                       Alpha = ale_df$Alpha)
  
  ale_ls[[index]] <- ale_df
  print(paste("Complete",index,"out of",length(var_ls)))
}

# Store all ale_df
out[[3]] <- ale_ls

# Output
saveRDS(out,paste0(Output_path,"RF_IM_",seed,".rds"))




