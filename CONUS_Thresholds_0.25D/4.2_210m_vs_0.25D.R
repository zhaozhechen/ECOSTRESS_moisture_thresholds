# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.3.4

# This code is to make maps for all variables at 0.25D, masked by 0.25D theta
# This code is to compare 0.25D thresholds directly from 0.25D data and aggregated from 210m
# This code is to explain 210m sd with diversity variables

##########
# Global
##########

dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)

library(raster)
library(sf)
library(ggplot2)
library(cowplot)
library(dplyr)
library(randomForest)
library(ggpointdensity)


# Input path of all rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Raw_var_rasters/var_all_combined/raster_all.rds")
# Make CONUS boundary
# Whole US map
CONUS <- st_read("/fs/ess/PAS2204/Data/cb_2018_us_state_20m/cb_2018_us_state_20m.shp")
# CONUS outer boundary map
#CONUS <- st_union(CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",])
CONUS <- CONUS[1][CONUS$STUSPS!="AK"&CONUS$STUSPS!="HI"&CONUS$STUSPS!="PR",]
# Output path of maps
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Figures/"

# Variable list to make maps
var_ls <- c("Theta_0.25D","Theta_agg210m","Theta_agg210m_sd",
            "T_gs","p_gs","vpd_gs","SM_gs",
            "CH","CH_sd","RD","LC","LC_shannon","DEM","DEM_sd",
            "Ks","n_fit","T_Sand",
            "g1","P50","gpmax","C")

###########
# Functions
###########
my_theme <- theme(
  axis.line=element_line(color="black"),
  panel.background = element_blank(),
  #text = element_text(size=100),
  panel.border = element_rect(colour="black",fill=NA),
  legend.key = element_blank(),
  #legend.key.size = unit(6,"cm"),
  #aspect.ratio = 1/1,
  #legend.key.size = unit(0.3,'cm'),
  #legend.title=element_text(size=100),
  #plot.title = element_text(margin=margin(0,0,80,0)),
  #plot.margin = margin(80,80,80,80),
  #axis.text.x = element_text(vjust=-3),
  #axis.text.y = element_text(margin=margin(0,50,0,0)),
  #legend.margin = margin(0,80,0,80)
  #axis.ticks.length = unit(-0.8,"cm")
)

plot_CM <- function(Corr_matrix){
  CM <- ggcorrplot::ggcorrplot(Corr_matrix,type="lower",lab=T,
                               colors=c("#619CFF","white","#F8766D"),
                               lab_size = 4,
                               tl.cex = 12)+
    my_theme+
    theme(aspect.ratio = 1/1)
  return(CM)
}

#######
# Main
#######

# Mask all layers with 0.25D theta
raster_all <- mask(raster_all,mask = raster_all$Theta)
# Rename rasters
names(raster_all)[names(raster_all)=="LC.1"] <- "LC"
names(raster_all)[names(raster_all)=="Theta"] <- "Theta_0.25D"

# Make maps ===========================================================================
g_all <- list()
for(i in 1:length(var_ls)){
  var_name <- var_ls[i]
  raster_df <- as.data.frame(raster_all[[var_name]],xy=TRUE)
  if(var_name != "LC"){
    g <- ggplot()+
      geom_tile(data=raster_df,aes(x=x,y=y,fill=.data[[var_name]]))+
      geom_sf(data=CONUS,fill=NA,color="grey",alpha=0.7)+
      my_theme+
      labs(x="",y="",title = var_name,fill="")
  }else if(var_name == "LC"){
    raster_df <- raster_df %>% 
      mutate(LC = replace(LC,LC == 41,"Deciduous"),
             LC = replace(LC,LC == 42,"Evergreen"),
             LC = replace(LC,LC == 43,"Mixed Forest"),
             LC = replace(LC,LC == 51,"Dwarf scrub"),
             LC = replace(LC,LC == 52,"Shurb/Scrub"),
             LC = replace(LC,LC == 71,"Grassland"),
             LC = replace(LC,LC == 72,"Sedge/Herbaceous"),
             LC = replace(LC,LC == 81,"Pasture/Hay"),
             LC = replace(LC,LC == 82,"Crops"))
    
    raster_df$LC <- as.factor(raster_df$LC)
    
    g <- ggplot()+
      geom_tile(data=raster_df,aes(x=x,y=y,fill=.data[[var_name]]))+
      geom_sf(data=CONUS,fill=NA,color="grey",alpha=0.7)+
      my_theme+
      labs(x="",y="",title=var_name,fill="")+
      scale_fill_discrete(na.value=NA)+
      theme(legend.position = "bottom")
  }
  
  # Use different color scales for different variables
  if(grepl("Theta",var_name)){
    g <- g + 
      scale_fill_gradient(na.value="white",
                          low="#0228FF",
                          high="#FFD902",
                          breaks = c(-10,20,50,80),
                          guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))+
      labs(title=var_name)
  }else if(var_name == "SM_gs_mean"){
    BrBG_color <- RColorBrewer::brewer.pal(11, "BrBG")
    BrBG_color <- c(BrBG_color[2:5],BrBG_color[7:10],rep(BrBG_color[10],3))
    g <- g + 
      scale_fill_gradientn(na.value="white",
                           colours = BrBG_color,
                           guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))
  }else if(var_name == "T_Sand"|var_name == "DEM"){
    BrBG_color <- RColorBrewer::brewer.pal(11, "BrBG")
    g <- g + 
      scale_fill_gradientn(na.value="white",
                           colours = BrBG_color,
                           guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))
  }else if(var_name != "LC"){
    g <- g + 
      scale_fill_gradient(na.value="white",
                          low="#01FFF0",
                          high="#FF0110",
                          guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))
  }
  
  # Add this figure to the figure list
  g_all[[i]] <- g
  print(paste("Complete Map",i))
}

# Put all maps together
g <- plot_grid(plotlist = g_all,ncol=3,align="hv",axis="tblr")
png(paste0(Output_path,"0.25D_Maps_all.png"),
    width=12,
    height=21,
    unit='in',
    res=600)
print(g)
dev.off()

# Compare 0.25D thresholds ======================================================
# Compare 0.25D thresholds directly from 0.25D data vs aggregated from 210m thresholds
# Color coded by 210m sd, LC, LC_shannon, CH_sd, and DEM_sd
var_ls <- c("Theta_agg210m_sd","LC_shannon","CH_sd","DEM_sd")
g_ls <- list()
for(i in 1:length(var_ls)){
  var_name <- var_ls[i]
  # Make a df
  df <- data.frame(Theta_0.25D = values(raster_all$Theta_0.25D),
                   Theta_agg210m = values(raster_all$Theta_agg210m),
                   var = values(raster_all[[var_name]]))
  df <- na.omit(df)
  # Fit a SLR
  lm <- lm(data=df,Theta_agg210m~Theta_0.25D)
  r2 <- round(cor(df$Theta_agg210m,df$Theta_0.25D)^2,2)
  # Make scatter plot, color coded by var
  g <- ggplot(data=df,aes(x=Theta_0.25D,y=Theta_agg210m,col=var))+
    geom_point(alpha=0.8)+
    my_theme+
    theme(aspect.ratio = 1/1)+
    scale_color_gradient(low="#01FFF0",
                         high="#FF0110",
                         guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))+
    geom_abline(intercept = 0,slope=1,linetype="dashed",color="black")+
    geom_abline(intercept = lm$coefficients[1],slope = lm$coefficients[2],color="red",linewidth=1)+
    labs(x="Alpha|0.25D",y="Alpha|210m aggregated",color=var_name)+
    annotate("text",label=paste0("r2 = ",r2),x=70,y=0)
  g_ls[[i]] <- g
}

# Make a scatter plot color coded by data density
# Make a df
df <- data.frame(Theta_0.25D = values(raster_all$Theta_0.25D),
                 Theta_agg210m = values(raster_all$Theta_agg210m))
df <- na.omit(df)
g <- ggplot(data=df,aes(x=Theta_0.25D,y=Theta_agg210m))+
  geom_point(alpha=0.8)+
  my_theme+
  geom_pointdensity()+
  theme(aspect.ratio = 1/1)+
  scale_color_gradient(low="#01FFF0",
                       high="#FF0110",
                       guide = guide_colorbar(frame.colour = "black",ticks.colour = "black"))+
  geom_abline(intercept = 0,slope=1,linetype="dashed",color="black")+
  geom_abline(intercept = lm$coefficients[1],slope = lm$coefficients[2],color="red",linewidth=1)+
  labs(x="Alpha|0.25D",y="Alpha|210m aggregated",color="Density")+
  annotate("text",label=paste0("r2 = ",r2),x=70,y=0)

g_ls[[5]] <- g

# Combine these plots
g <- plot_grid(plotlist = g_ls,nrow=3,align="hv")
png(paste0(Output_path,"0.25D_Theta_comparison.png"),
    height = 8,
    width = 8,
    unit='in',
    res=600)
print(g)
dev.off()

# Test correlation among several variables ===================================================
var_ls <- c("Theta_0.25D","Theta_agg210m","Theta_agg210m_sd",
            "T_gs","p_gs","vpd_gs","SM_gs",
            "CH","CH_sd","RD","LC","LC_shannon","DEM","DEM_sd",
            "Ks","n_fit","T_Sand",
            "g1","P50","gpmax","C")
df <- as.data.frame(raster_all)
df <- df[var_ls]
df <- na.omit(df)
# Calculate difference between 210m aggregated theta and 0.25D theta
df$Dif <- df$Theta_0.25D-df$Theta_agg210m
corr_matrix <- cor(df,use="pairwise.complete.obs",method="pearson")
g <- plot_CM(corr_matrix)
png(paste0(Output_path,"Diversity_CM.png"),
    height = 12,
    width=12,
    unit='in',
    res=600)
print(g)
dev.off()

# RF with Dif as the response variable ===================================================
# Explanatory variables for 210m_sd
df_rf <- df[,-c(1,2,3,11)]
set.seed(1)
# Split dataset into Training set 70% and Testing set 30%
train_index <- sample(1:nrow(df_rf),size=nrow(df_rf)*0.7)
df_train  <- df_rf[train_index,]
df_test   <- df_rf[-train_index,]

print('Fitting RF model ...')
# # of trees: 500
set.seed(1)
system.time(rf <- randomForest(Dif~.,
                               data = df_train,
                               importance=TRUE,
                               type = "regression",
                               ntree=500))
print("Complete RF model")

r_train <- cor(rf$predicted,df_train$Dif,use="pairwise.complete.obs")
rf.test.pred <- predict(rf,df_test)
r_test  <- cor(rf.test.pred,df_test$Dif,use="pairwise.complete.obs")
rf_importance <- as.data.frame(importance(rf))
rf_importance$Feature <- rownames(rf_importance)
rf_importance <- rf_importance %>% rename("Importance" = "%IncMSE")
# Rescale importance to sum = 1
rf_importance$Importance <- rf_importance$Importance/sum(rf_importance$Importance)
# Plot of relative importance
g_imp <- ggplot(data=rf_importance,aes(x=Importance,y=reorder(Feature,Importance)))+
  geom_bar(stat="identity",fill="skyblue",color="black")+
  my_theme+
  labs(y="")+
  annotate("text",label = paste0("r2_train = ",round(r_train^2,3)),x=0.07,3)+
  annotate("text",label = paste0("r2_test = ",round(r_test^2,3)),x=0.07,2)+
  ggtitle("Diff (Alpha|0.25D - Alpha|210m aggregated)")

# Approximate response curves
# Initialize a list to store all pdp
pdp_ls <- list()
# variables to plot
var_ls <- names(df_train)[-18]
for(index in 1:length(var_ls)){
  # Name of the changing variable
  var_chg_name <- var_ls[index]
  # Create a sequence of values for this changing variable
  var_chg <- seq(min(df_train[var_chg_name]),
                 max(df_train[var_chg_name]),
                 length.out = 50)
  # Name of the fixing variables
  var_fix_name <- var_ls[var_ls!=var_chg_name]
  # Create a dummy training set, fix all fixing variables at median
  train_dummy <- sapply(df_train[var_fix_name],median)
  # Repeat 50 times
  train_dummy <- as.data.frame(lapply(train_dummy,rep,50))
  # Add the changing variable
  train_dummy$var <- var_chg
  colnames(train_dummy)[colnames(train_dummy)=="var"] <- var_chg_name
  # Use this dummy predictor set to predict Diff
  Theta_dummy <- predict(rf,train_dummy)
  pdp_df <- data.frame(Diff = Theta_dummy,
                       var_chg)
  # Get data density
  density <- density(df_train[var_chg_name][,1],n=50)$y
  pdp_df$density <- density
  # Make PDP plot
  g_PDP <- ggplot()+
    geom_line(data = pdp_df,aes(x = var_chg,y=Diff,color=density))+
    my_theme+
    theme(aspect.ratio=1/1.5,
          legend.position="none")+
    scale_color_gradient2(low="grey99",high="#0228FF")+
    labs(x=var_chg_name,y="Diff")
  # Put this PDP into the list
  pdp_ls[[index]] <- g_PDP
  print(index)
}
# Combine all plots
g_all <- plot_grid(plotlist = pdp_ls,ncol=4)
g_all <- plot_grid(g_imp,g_all,
                   rel_widths=c(1,2.5),
                   axis="t")
png(paste0(Output_path,"0.25D_Diff_RF_Plots.png"),
    height=12,
    width=16,
    unit='in',
    res=600)
print(g_all)
dev.off()



