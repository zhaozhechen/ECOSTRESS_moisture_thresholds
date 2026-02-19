# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.3.29

# This code is to run RF for 100 times
# Response variable is 0.25D alpha aggregated from 210m alpha
# Output include IP df, R2 of RF model, and PDP df for all variables
# Note: The output PDP are calculated using partialPlot function in randomForest pacakge
# Predictors have no CH_sd or aridity index (AI), keep CH

#########
# Global
#########

library(ggplot2)
library(randomForest)
library(dplyr)
library(cowplot)
library(raster)

# Set the random seed
seed <- 1

# Input path of 0.25 degree data aggregated from 210m
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")
# Output of intermediate files, including IP df, PDP df, R2
#Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/RF_210m_aggregated/IM_files_V2/"

# Response variable
var_re <- c("Theta")
# List of predicting variables
var_pr <- c("T_gs_mean","Ks","T_Sand","n_parameter",
            "vpd_gs_mean","CH","p_gs_mean","RD","DEM_sd",
            "LC_shannon","LC","SM_gs_mean","g1","DEM",
            "gpmax","C","P50")
var_ls_all <- c(var_re,var_pr)

########
# Main
########

# Only keep required variables
raster_all <- raster::subset(raster_all,var_ls_all)
df_all <- as.data.frame(raster_all)
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

# Only keep low alpha (VPD driven)
df_all <- df_all[df_all$Theta <= median(df_all$Theta),]

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

# RF IP ===========================================
g_IP <- ggplot(rf_importance,aes(x=Importance,y=reorder(Feature,Importance)))+
  geom_bar(stat="identity",color="black")


