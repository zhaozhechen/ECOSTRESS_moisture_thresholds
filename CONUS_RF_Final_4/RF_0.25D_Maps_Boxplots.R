# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.3.28

# This code is to output boxplots grouped by Aridity index (AI) 
# for 0.25D thresholds maps (aggregated from 210m)

# AI classification: Accelerated dryland expansion under climate change
# Hyper-arid: AI < 0.05
# Arid: 0.05 <= AI < 0.2
# Semiarid: 0.2 <= AI < 0.5
# Dry subhumid: 0.5 <= AI < 0.65
# Humid: AI >= 0.65

#########
# Global
#########
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)

library(raster)
library(sf)
library(ggplot2)
library(dplyr)
library(cowplot)
library(viridis)
library(ggsignif)

# Input path of combined 0.25D rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Figures/"

###########
# Functions
###########
my_theme <- theme(
  axis.line=element_line(color="black"),
  panel.background = element_blank(),
  text = element_text(size=14),
  panel.border = element_rect(colour="black",fill=NA),
  legend.key = element_blank(),
  aspect.ratio = 1/1,
  #legend.key.size = unit(0.3,'cm'),
  legend.title=element_text(size=12),
  axis.title = element_text(size=12),
  legend.position = "none"
  #plot.margin = margin(0,0,0,0,'cm'),
)

#######
# Main
#######

# Mask all layers with Theta
raster_all <- mask(raster_all,mask = raster_all$Theta)

# Get a df of three metrics and AI
df <- data.frame(Alpha = values(raster_all$Theta),
                 VPD50 = values(raster_all$SM_50_abs),
                 SM50 = values(raster_all$VPD_50_abs),
                 AI = values(raster_all$AI))
df <- na.omit(df)
df$Type <- rep(NA,nrow(df))
# Get dryland type classification
df <- df %>%
  mutate(Type = replace(Type,AI < 0.05,"Hyper-arid")) %>%
  mutate(Type = replace(Type,AI >= 0.05 & AI < 0.2,"Arid")) %>%
  mutate(Type = replace(Type,AI >= 0.2 & AI < 0.5,"Semiarid")) %>%
  mutate(Type = replace(Type,AI >= 0.5 & AI <= 0.65,"Subhumid")) %>%
  mutate(Type = replace(Type,AI >= 0.65,"Humid"))
df$Type <- as.factor(df$Type)
levels(df$Type) <- c("Hyper-arid","Arid","Semiarid","Subhumid","Humid")

# Define Dryland or Humid region
df$Dry <- rep(NA,nrow(df))
df <- df %>%
  mutate(Dry = replace(Dry,AI<0.5,"Dryland")) %>%
  mutate(Dry = replace(Dry,AI>=0.5,"Humid"))
df$Dry <- as.factor(df$Dry)

# Make boxplot grouped by type ==================
# Initialize a list to store boxplots
BP_ls <- list()
for(i in 1:3){
  var_name <- names(df)[i]
  g <- ggplot(data=df,aes(x=Type,y=.data[[var_name]],fill = Type))+
    geom_boxplot()+
    scale_fill_brewer(palette = "RdYlGn")+
    my_theme+
    labs(x="")
  if(i == 1){
    g <- g+
      labs(y = bquote(alpha~ (degree)))
  }else if(i == 2){
    g <- g+
      labs(y = bquote(VPD[50]~" (kPa)"))
  }else if(i == 3){
    g <- g+
      labs(y = bquote(SM[50]~ (m^3/m^3)))
  }
  BP_ls[[i]] <- g
  print(i)
}

# Make boxplot grouped by Dryland ==================
# Initialize a list to store boxplots
BP_ls2 <- list()
my_color <- RColorBrewer::brewer.pal(11,"RdYlGn")
my_color <- c(my_color[3],my_color[10])
for(i in 1:3){
  var_name <- names(df)[i]
  g <- ggplot(data=df,aes(x=Dry,y=.data[[var_name]],fill = Dry))+
    geom_boxplot(fill=my_color)+
    geom_signif(
      comparisons = list(c("Dryland","Humid")),
      map_signif_level = TRUE,
      test="wilcox.test",
      margin_top = 0.2,
    )+
    my_theme+
    labs(x="")
  if(i == 1){
    g <- g+
      labs(y = bquote(alpha~ (degree)))
  }else if(i == 2){
    g <- g+
      labs(y = bquote(VPD[50]~" (kPa)"))
  }else if(i == 3){
    g <- g+
      labs(y = bquote(SM[50]~ (m^3/m^3)))
  }
  BP_ls2[[i]] <- g
  print(i)
}

# Combine boxplots
g <- plot_grid(plotlist = c(BP_ls,BP_ls2),nrow=2)
pdf(paste0(Output_path,"0.25D_210m_agg_boxplots.pdf"),
    height=6,
    width=14)
print(g)
dev.off()




