# Author: Zhaozhe Chen (chen.8926@osu.edu)
# Date: 2024.3.10

# This code is to check RF variables correlations for 0.25D and 210m 
# Only keep complete observations

#########
# Global
#########

library(raster)
library(ggplot2)

dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/gdal/3.2.1/lib/libgdal.so")
dyn.load("/apps/R/gnu/9.1/4.0.2/site/lib/geos/3.8.1/lib/libgeos_c.so", local=FALSE)

library(raster)
library(sf)

# Manually change 0.25D or 210m
Resolution <- "0.25D"

# Whole raster for CONUS
raster_all <- readRDS(paste0("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/",Resolution,"_raster_combined/Combined_raster.rds"))

# Variables to test
var_ls <- c("Theta",
            "T_gs_mean","p_gs_mean","vpd_gs_mean","SM_gs_mean","AI","LAI",
            "CH","CH_sd","RD","LC","LC_shannon","DEM","DEM_sd",
            "Ks","n_parameter","T_Sand",
            "g1","P50","gpmax","C")

Output_path <- paste0("/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CONUS_RF_Final_4/Figures/",Resolution,"/")

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

# Only keep required variables
raster_all <- raster_all[[var_ls]]
df <- as.data.frame(raster_all)
# Only keep complete observations
df <- na.omit(df)
corr_matrix <- cor(df,use="pairwise.complete.obs",method="pearson")
g <- plot_CM(corr_matrix)

png(paste0(Output_path,"Correlation_matrix.png"),
    width=12,
    height=12,
    res=600,
    unit='in')
print(g)
dev.off()



