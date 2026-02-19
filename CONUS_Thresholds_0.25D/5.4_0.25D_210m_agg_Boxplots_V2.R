# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.4.15
# Update: 2025.2.15

# This code is to output boxplots grouped by Aridity index (AI)
# For 0.25D thresholds maps (aggregated from 210m)
# For Figure 2
# Including theta (alpha), SM50 and VPD50 at both quantiles and their original scales

# Also output boxplots, and half violin plots for metrics grouped by PFT

# AI classification: Accelerated dryland expansion under climate change
# Hyper-arid: AI < 0.05
# Arid: 0.05 <= AI < 0.2
# Semiarid: 0.2 <= AI < 0.5
# Dry subhumid: 0.5 <= AI < 0.65
# Humid: AI >= 0.65

# Group hyper-arid and arid together
# Arid: AI < 0.2
# Semiarid: 0.2 <= AI < 0.5
# Semihumid: 0.5 <= AI < 0.65
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
library(gghalves)

# Input path of combined 0.25D rasters
raster_all <- readRDS("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/0.25D_raster_combined/Combined_raster.rds")
Output_path <- "/fs/ess/PAS2204/Results/CONUS_Threshold_0.25D/Figures/"

###########
# Functions
###########
my_theme <- theme(
  axis.line=element_line(color="black"),
  panel.background = element_blank(),
  text = element_text(size=16),
  panel.border = element_rect(colour="black",fill=NA),
  legend.key = element_blank(),
  #aspect.ratio = 1/2,
  #legend.key.size = unit(0.3,'cm'),
  legend.title=element_text(size=16),
  axis.title = element_text(size=16),
  axis.text = element_text(size=16),
  legend.position = "none"
  #plot.margin = margin(0,0,0,0,'cm'),
)

#######
# Main
#######

# Mask all layers with Theta
raster_all <- mask(raster_all,mask = raster_all$Theta)
# Adjust VPD_50_abs for color scale
raster_all$VPD_50_abs[raster_all$VPD_50_abs>=0.5] <- 0.5

# Get a df of three metrics and AI
df <- data.frame(Alpha = values(raster_all$Theta),
                 VPD50qt = values(raster_all$SM_50th),
                 SM50qt = values(raster_all$VPD_50th),
                 VPD50abs = values(raster_all$SM_50_abs),
                 SM50abs = values(raster_all$VPD_50_abs),
                 AI = values(raster_all$AI))
df <- na.omit(df)
df$Type <- rep(NA,nrow(df))
# Get dryland type classification
df$Type[df$AI < 0.2] <- "1"
df$Type[df$AI >= 0.2 & df$AI < 0.5] <- "2"
df$Type[df$AI >= 0.5 & df$AI < 0.65] <- "3"
df$Type[df$AI >= 0.65] <- "4"

df$Type <- as.factor(df$Type)

levels(df$Type) <- c("Arid","Semiarid","Semihumid","Humid")

# Define Dryland or Humid region
df$Dry <- rep(NA,nrow(df))
df <- df %>%
  mutate(Dry = replace(Dry,AI < 0.2,"Arid")) %>%
  mutate(Dry = replace(Dry,AI >= 0.2,"Humid"))
df$Dry <- as.factor(df$Dry)

# Make boxplot grouped by type (aridify index AI) ==================
# Initialize a list to store boxplots
BP_ls <- list()
for(i in 1:5){
  var_name <- names(df)[i]
  g <- ggplot(data=df,aes(y=Type,x=.data[[var_name]],fill = Type))+
    geom_boxplot(outlier.size = 1,
                 outlier.alpha = 0.3,
                 outlier.shape = 23)+
    scale_fill_brewer(palette = "Blues",direction = -1)+
    my_theme+
    theme(aspect.ratio = 1/1.2)+
    labs(y="")
  #scale_y_discrete(guide = guide_axis(n.dodge = 2))
  if(i == 1){
    g <- g+
      labs(x = bquote(alpha~ (degree)))
  }else if(i == 2){
    g <- g+
      labs(x = bquote(VPD[50]))
  }else if(i == 3){
    g <- g+
      labs(x = bquote(SM[50]))
  }else if(i == 4){
    g <- g+
      labs(x = bquote(VPD[50] (kPa)))
  }else if(i == 5){
    g <- g+
      labs(x = bquote(SM[50] (m^3/m^3)))
  }
  BP_ls[[i]] <- g
  print(i)
}

# Combine boxplots
g <- plot_grid(plotlist = BP_ls,nrow=2)
pdf(paste0(Output_path,"0.25D_210m_agg_boxplots_AI_20250215.pdf"),
    height=6,
    width=10)
print(g)
dev.off()

# Make Panel d =======================================
# Boxplots grouped by PFT
df <- data.frame(Alpha = values(raster_all$Theta),
                 VPD50qt = values(raster_all$SM_50th),
                 SM50qt = values(raster_all$VPD_50th),
                 VPD50abs = values(raster_all$SM_50_abs),
                 SM50abs = values(raster_all$VPD_50_abs),
                 PFT = values(raster_all$LC))

df <- na.omit(df)
df <- df %>%
  mutate(PFT = replace(PFT,PFT==41,"DF"),
         PFT = replace(PFT,PFT==42,"EF"),
         PFT = replace(PFT,PFT==43,"MF"),
         PFT = replace(PFT,PFT==51,"Dwarf scrub"),
         PFT = replace(PFT,PFT==52,"SHB"),
         PFT = replace(PFT,PFT==71,"GRA"),
         PFT = replace(PFT,PFT==72,"Sedge/Herbaceous"),
         PFT = replace(PFT,PFT==81,"PAS"),
         PFT = replace(PFT,PFT==82,"CRO"))
df$PFT <- factor(df$PFT,levels = c("CRO","SHB","GRA","PAS","EF","DF","MF"))

# Initialize a list to store boxplots
BP_ls <- list()
for(i in 1:5){
  var_name <- names(df)[i]
  g <- ggplot(data=df,aes(x=PFT,y=.data[[var_name]],fill = PFT))+
    geom_boxplot(outlier.size = 0.2,
                 outlier.alpha = 0.3)+
    my_theme+
    theme(aspect.ratio = 1/0.85,
          axis.text.x = element_blank())+
    labs(x="")
  #scale_y_discrete(guide = guide_axis(n.dodge = 2))
  if(i == 1){
    g <- g+
      labs(y = bquote(alpha~ (degree)))
  }else if(i == 2){
    g <- g+
      labs(y = bquote(VPD[50]),fill="")
  }else if(i == 3){
    g <- g+
      labs(y = bquote(SM[50]))
  }else if(i == 4){
    g <- g+
      labs(y = bquote(VPD[50] (kPa)))
  }else if(i == 5){
    g <- g+
      labs(y = bquote(SM[50] (m^3/m^3)))+
      theme(legend.position = "bottom")
  }
  BP_ls[[i]] <- g
  print(i)
}
# Combine boxplots
g <- plot_grid(plotlist = BP_ls,nrow=2,
               align = "hv",
               axis = "tblr")
pdf(paste0(Output_path,"0.25D_210m_agg_boxplots_PFT_20250211.pdf"),
    height=10,
    width=8.5)
print(g)
dev.off()

# Make Panel d 2025.2.15 version =======================================
# Half violin plots grouped by PFT, conditional on aridity
df <- data.frame(Alpha = values(raster_all$Theta),
                 VPD50abs = values(raster_all$SM_50_abs),
                 SM50abs = values(raster_all$VPD_50_abs),
                 PFT = values(raster_all$LC),
                 AI = values(raster_all$AI))

df <- na.omit(df)
df <- df %>%
  mutate(PFT = replace(PFT,PFT==41,"DF"),
         PFT = replace(PFT,PFT==42,"EF"),
         PFT = replace(PFT,PFT==43,"MF"),
         PFT = replace(PFT,PFT==51,"Dwarf scrub"),
         PFT = replace(PFT,PFT==52,"SHB"),
         PFT = replace(PFT,PFT==71,"GRA"),
         PFT = replace(PFT,PFT==72,"Sedge/Herbaceous"),
         PFT = replace(PFT,PFT==81,"PAS"),
         PFT = replace(PFT,PFT==82,"CRO"))
df$PFT <- factor(df$PFT,levels = c("CRO","SHB","GRA","PAS","EF","DF","MF"))

# Semihumid and humid are considered humid
# Semiarid and arid are considered arid

# Get dryland type classification
df$Type[df$AI < 0.2] <- "1"
df$Type[df$AI >= 0.2 & df$AI < 0.5] <- "2"
df$Type[df$AI >= 0.5 & df$AI < 0.65] <- "3"
df$Type[df$AI >= 0.65] <- "4"

df$Type <- as.factor(df$Type)
levels(df$Type) <- c("Arid","Semiarid","Semihumid","Humid")

df <- df %>%
  mutate(AI_type = case_when(
    Type == "Humid"| Type == "Semihumid" ~ "Humid",
    Type == "Arid"|Type == "Semiarid" ~ "Arid",
    TRUE ~ "Other"))

# Make half violin plots
VP_ls <- list()
for(i in 1:3){
  var_name <- names(df)[i]
  g <- ggplot()+
    # Top side is for Humid
    geom_half_violin(data=df[df$AI_type == "Humid",],
                     aes(x=PFT,y=.data[[var_name]]),
                     fill="#EFF3FF",
                     side="r",
                     draw_quantiles = 0.5,
                     size=0.3)+
    # Bottom side is for arid
    geom_half_violin(data=df[df$AI_type == "Arid",],
                     aes(x=PFT,y=.data[[var_name]]),
                     fill="#6BAED6",
                     side="l",
                     draw_quantiles = 0.5,
                     size=0.3)+
    my_theme+
    theme(aspect.ratio = 1/0.85)+
    labs(x="")+
    coord_flip()
  
  if(i == 1){
    g <- g+
      labs(y = bquote(alpha~ (degree)))
  }else if(i == 2){
    g <- g+
      labs(y = bquote(VPD[50] (kPa)))
  }else if(i == 3){
    g <- g+
      labs(y = bquote(SM[50] (m^3/m^3)))
  }  
  VP_ls[[i]] <- g
}

g <- plot_grid(plotlist = VP_ls,nrow=1,
               align = "hv",
               axis = "tblr")

pdf(paste0(Output_path,"0.25D_210m_agg_violin_PFT_20250215.pdf"),
    height=5,
    width=9.5)
print(g)
dev.off()


