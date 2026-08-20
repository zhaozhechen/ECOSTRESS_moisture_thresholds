# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2026.8.19

# This code is to repeat the random forest attribution analysis using SHAP
# Response variable is 0.25D alpha aggregated from 210m alpha

# ---- Global --------
library(ggplot2)
library(randomForest)

# Path to input df
Input_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/00_Data/Processed/RF_SHAP/df_all.csv"
# Path to output figures
Output_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/02_Results/Figures"
# Path to output tables
Table_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/02_Results/Tables"

# Response variable
var_re <- "Theta"
# List of predicting variables
var_pr <- c("T_gs_mean","Ks","T_Sand","n_parameter",
            "vpd_gs_mean","CH","p_gs_mean","RD","DEM_sd",
            "LC_shannon","LC","SM_gs_mean","DEM")

# Use the same random seed as the first RF run in the submitted analysis
seed <- 1
# Use a subset of held-out testing observations for SHAP
n_eval <- 2000
# Number of Monte Carlo permutations used for approximate SHAP
n_sim <- 10

# ---- Functions --------
my_theme <- theme(
  axis.line=element_line(color="black"),
  panel.background=element_blank(),
  text=element_text(size=14),
  panel.border=element_rect(colour="black",fill=NA),
  legend.key=element_blank(),
  legend.title=element_text(size=12),
  axis.title=element_text(size=12)
)

# Approximate SHAP values using random feature permutations
Approximate_SHAP <- function(rf,background,newdata,n_sim){
  n <- nrow(newdata)
  p <- ncol(newdata)
  shap <- matrix(0,nrow=n,ncol=p)
  colnames(shap) <- names(newdata)

  for(sim in 1:n_sim){
    background_id <- sample(1:nrow(background),n,replace=TRUE)
    background_tmp <- background[background_id,,drop=FALSE]
    random_order <- matrix(runif(n*p),nrow=n,ncol=p)

    for(j in 1:p){
      before <- background_tmp
      for(k in 1:p){
        replace_id <- random_order[,k] < random_order[,j]
        before[replace_id,k] <- newdata[replace_id,k]
      }
      after <- before
      after[,j] <- newdata[,j]
      shap[,j] <- shap[,j] + predict(rf,after) - predict(rf,before)
    }
  }

  shap/n_sim
}

Get_variable_name <- function(x){
  name_ls <- c(T_gs_mean="Tair",
               Ks="Ks",
               T_Sand="Sand fraction",
               n_parameter="n",
               vpd_gs_mean="VPD",
               CH="Canopy height",
               p_gs_mean="P",
               RD="Root depth",
               DEM_sd="Topographic diversity",
               LC_shannon="LC diversity",
               LC="LC",
               SM_gs_mean="SM",
               DEM="Elevation")
  unname(name_ls[x])
}

# ---- Main --------
dir.create(Output_path,recursive=TRUE,showWarnings=FALSE)
dir.create(Table_path,recursive=TRUE,showWarnings=FALSE)

df_all <- read.csv(Input_path)
df_all <- df_all[c(var_re,var_pr)]
df_all <- na.omit(df_all)
df_all$LC <- as.factor(df_all$LC)

# Split dataset into Training set 70% and Testing set 30%
set.seed(seed)
train_index <- sample(1:nrow(df_all),size=floor(nrow(df_all)*0.7))
df_train <- df_all[train_index,]
df_test <- df_all[-train_index,]

# Fit the same 500-tree RF model used in the submitted analysis
set.seed(seed)
f <- as.formula(paste(var_re,"~.",sep=""))
rf <- randomForest(f,
                   data=df_train,
                   importance=TRUE,
                   type="regression",
                   ntree=500)

# Evaluate model performance on the separated testing set
rf_test_pred <- predict(rf,df_test)
R2_train <- cor(rf$predicted,df_train[[var_re]],use="pairwise.complete.obs")^2
R2_test <- cor(rf_test_pred,df_test[[var_re]],use="pairwise.complete.obs")^2

# Select a reproducible subset only from the held-out testing set
set.seed(10000+seed)
eval_index <- sample(1:nrow(df_test),min(n_eval,nrow(df_test)))
df_eval <- df_test[eval_index,]

# Calculate SHAP values for the held-out subset
set.seed(20000+seed)
shap <- Approximate_SHAP(rf,
                         background=df_train[var_pr],
                         newdata=df_eval[var_pr],
                         n_sim=n_sim)

# Mean absolute SHAP value for global feature importance
SHAP_importance <- colMeans(abs(shap))
SHAP_importance <- SHAP_importance/sum(SHAP_importance)

# Original permutation importance from the same fitted model for comparison
Permutation_importance <- importance(rf)[,"%IncMSE"]
Permutation_importance <- Permutation_importance/sum(Permutation_importance)

importance_df <- data.frame(
  Feature=names(SHAP_importance),
  Variable=Get_variable_name(names(SHAP_importance)),
  SHAP_relative_importance=as.numeric(SHAP_importance),
  Permutation_relative_importance=as.numeric(Permutation_importance[names(SHAP_importance)])
)
importance_df$SHAP_rank <- rank(-importance_df$SHAP_relative_importance,ties.method="min")
importance_df$Permutation_rank <- rank(-importance_df$Permutation_relative_importance,ties.method="min")

importance_df$Group <- NA
importance_df$Group[importance_df$Feature %in% c("T_gs_mean","p_gs_mean","vpd_gs_mean","SM_gs_mean")] <- "Long-term climate"
importance_df$Group[importance_df$Feature %in% c("CH","RD","LC","LC_shannon")] <- "Ecological condition"
importance_df$Group[importance_df$Feature %in% c("Ks","T_Sand","n_parameter")] <- "Soil property"
importance_df$Group[importance_df$Feature %in% c("DEM","DEM_sd")] <- "Topography"

# Plot the SHAP relative importance in the same style as Figure 3a
g_SHAP <- ggplot(importance_df,
                 aes(x=SHAP_relative_importance,
                     y=reorder(Variable,SHAP_relative_importance)))+
  geom_bar(stat="identity",color="black",aes(fill=Group))+
  my_theme+
  scale_fill_brewer(palette="Set3")+
  labs(y="",fill="Category",x="Relative mean |SHAP value|")+
  theme(legend.position="inside",
        legend.position.inside=c(0.76,0.17))

pdf(paste0(Output_path,"/RF_SHAP_importance.pdf"),height=6,width=5.5)
print(g_SHAP)
dev.off()

png(paste0(Output_path,"/RF_SHAP_importance.png"),
    height=6,width=5.5,units="in",res=600)
print(g_SHAP)
dev.off()

write.csv(importance_df[order(importance_df$SHAP_rank),],
          paste0(Table_path,"/RF_SHAP_importance.csv"),row.names=FALSE)

group_df <- aggregate(SHAP_relative_importance~Group,
                      data=importance_df,
                      FUN=sum)
group_df <- group_df[order(group_df$SHAP_relative_importance,decreasing=TRUE),]
write.csv(group_df,
          paste0(Table_path,"/RF_SHAP_group_importance.csv"),row.names=FALSE)

summary_df <- data.frame(
  Seed=seed,
  Complete_observations=nrow(df_all),
  Training_observations=nrow(df_train),
  Testing_observations=nrow(df_test),
  SHAP_evaluation_observations=nrow(df_eval),
  SHAP_permutations=n_sim,
  R2_train=R2_train,
  R2_test=R2_test,
  Rank_correlation=cor(importance_df$SHAP_rank,
                       importance_df$Permutation_rank,
                       method="spearman")
)
write.csv(summary_df,
          paste0(Table_path,"/RF_SHAP_summary.csv"),row.names=FALSE)

print(importance_df[order(importance_df$SHAP_rank),])
print(summary_df)
