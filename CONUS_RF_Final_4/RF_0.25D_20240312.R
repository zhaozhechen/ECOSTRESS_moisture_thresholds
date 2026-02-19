# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.3.12

# This code is to run RF
# 210m aggregated alpha is response variable
# Add aridity index (AI) and remove canopy height and CH_sd for predictors
# Note: use partialPlot function in randomForest package to calculate PDP

#########
# Global
#########

library(ggplot2)
library(randomForest)
library(dplyr)
library(cowplot)

# Set the random seed
seed <- 1
#seed <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID')) 

# Input path of 0.25 degree data aggregated from 210m
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")
# Output path of figures
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/Figures/0.25D/"

# Response variable
var_re <- "Theta"
var_pr <- c("T_gs_mean","Ks","T_Sand","n_parameter",
            "vpd_gs_mean","RD","AI",
            "LC_shannon","LC","SM_gs_mean","g1","DEM","DEM_sd",
            "gpmax","C","P50")
var_ls_all <- c(var_re,var_pr)

###########
# Functions
###########

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

# This function is to output plots
# If type = 1, this is importance plot
# If type = 2, this is PDP plot
plotp <- function(g,type,title){
  if(type == 1){
    x <-4
    y <-6
  }else{
    x <-24
    y <-12
  }
  png(paste(Output_path,"/",title,".png",sep=""),
      width=x,
      height=y,
      unit="in",
      res=600)
  print(g)
  dev.off()
}

########
# Main
########

# Only keep required variables
raster_all <- raster_all[[var_ls_all]]
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

# Run RF  =====================================================
# Split dataset into Training set 70% and Testing set 30%
set.seed(seed)
train_index <- sample(1:nrow(df_all),size=nrow(df_all)*0.7)
df_train  <- df_all[train_index,]
df_test   <- df_all[-train_index,]

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
rf_importance <- as.data.frame(importance(rf))
rf_importance$Feature <- rownames(rf_importance)
rf_importance <- rf_importance %>% rename("Importance" = "%IncMSE")
# Rescale importance to sum = 1
rf_importance$Importance <- rf_importance$Importance/sum(rf_importance$Importance)
# Plot of relative importance
g_imp <- ggplot(data=rf_importance,aes(x=Importance,y=reorder(Feature,Importance)))+
  geom_bar(stat="identity",fill="skyblue",color="black")+
  annotate("text",x=0.8*max(rf_importance$Importance),y=3,label=paste("R2_train =",round(r_train^2,2)))+
  annotate("text",x=0.8*max(rf_importance$Importance),y=2,label=paste("R2_test =",round(r_test^2,2)))+
  my_theme+
  labs(y="")
#plotp(g_imp,1,"RF_IP_AI+P")

# PDP
# Initialize a list to store all all 16 variables (no LC)
pdp_ls <- list()
IP_var_ls <- rownames(rf_importance[order(rf_importance$Importance,decreasing = TRUE),])
IP_var_ls <- IP_var_ls[IP_var_ls!="LC"]
for(i in 1:length(IP_var_ls)){
  var_name <- IP_var_ls[i]
  test <- partialPlot(rf,df_train,IP_var_ls[i])
  test <- data.frame(x=test$x,y=test$y)
  # Calculate data density
  test$density <- density(df_train[var_name][,1],n=nrow(test))$y
  # Make PDP plot
  g_pdp <- ggplot(test,aes(x,y))+
    geom_line(aes(color=density),size=2)+
    my_theme+
    theme(aspect.ratio = 1/2)+
    labs(x=var_name,y=var_re)+
    ylim(c(25.5,33))+
    theme(aspect.ratio = 1/1.5,
          legend.position = "none")+
    scale_color_gradient2(low="grey99",high="#0228FF")
  pdp_ls[[i]] <- g_pdp
  print(paste("Complete PDP",i))
}

# Combine all plots
g_all <- plot_grid(plotlist = pdp_ls,ncol=4)
g_all <- plot_grid(g_imp,g_all,
                   rel_widths = c(1,2.5),
                   axis="t")

png(paste0(Output_path,"RF_AI_noP.png"),
    height=12,
    width=15,
    unit='in',
    res=600)
print(g_all)
dev.off()


