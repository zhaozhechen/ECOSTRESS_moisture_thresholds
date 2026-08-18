# Author: Zhaozhe Chen (zhaozhe.chen@wisc.edu)
# Date: 2026.8.18

# This code is to show the fraction of future projected days outside
# the pixel-specific historical ranges of VPD and SM

# ---- Global --------
library(ggplot2)
library(cowplot)
library(RColorBrewer)

# Path to input table
Input_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/02_Results/Tables/Future_outside_historical_range_by_model.csv"
# Path to output figures
Output_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/02_Results/Figures"

# Source plotting functions
source("D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes/01_Codes/Plotting_functions.R")

# Colors used for VPD and SM in the manuscript
my_color <- brewer.pal(6,"Set2")

VPD_color <- my_color[2]
SM_color <- my_color[1]

# ------ Main -------
dir.create(Output_path,recursive=TRUE,showWarnings=FALSE)

# Read model-level statistics
df <- read.csv(Input_path)
# Only keep end-century projections
df <- df[df$Time=="End",]
df$Scenario <- ifelse(df$SSP=="ssp245","SSP2-4.5","SSP5-8.5")
df$Scenario <- factor(df$Scenario,levels=c("SSP2-4.5","SSP5-8.5"))
df$Group <- factor(df$Scenario,levels=c("SSP2-4.5","SSP5-8.5"))

# Panel A: fraction outside either historical range
df_either <- df[df$Metric=="Either_outside",]
df_either$Percent <- df_either$Mean_fraction*100

df_either_mean <- aggregate(Percent~Group,data=df_either,FUN=mean)
df_either_mean$Within <- 100-df_either_mean$Percent

df_bar <- rbind(
  data.frame(Group=df_either_mean$Group,
             Range="Outside historical range",
             Percent=df_either_mean$Percent),
  data.frame(Group=df_either_mean$Group,
             Range="Within historical range",
             Percent=df_either_mean$Within)
)
df_bar$Range <- factor(df_bar$Range,
                       levels=c("Within historical range",
                                "Outside historical range"))

g1 <- ggplot()+
  geom_col(data=df_bar,
           aes(x=Group,y=Percent,fill=Range),
           width=0.65,color="black",linewidth=0.6)+
  geom_point(data=df_either,
             aes(x=Group,y=Percent),
             position=position_jitter(width=0.13,height=0),
             size=1.5,alpha=0.65)+
  geom_text(data=df_either_mean,
            aes(x=Group,y=Percent/2,
                label=paste0(round(Percent,1),"%")),
            size=4.5)+
  scale_fill_manual(values=c("Within historical range"="grey85",
                             "Outside historical range"=my_color[3]),
                    guide=guide_legend(override.aes=list(color="black",
                                                        linewidth=0.6)))+
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20),
                     expand=c(0,0))+
  my_theme+
  theme(text=element_text(size=12),
        axis.title=element_text(size=13),
        axis.text=element_text(size=11),
        axis.text.x=element_text(size=10),
        legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=11),
        plot.title=element_text(size=14))+
  labs(x="End-century scenario",y="Projected days (%)")+
  ggtitle("a  Outside either historical range")

# Panel B: separate contributions of VPD and SM
df_variable <- df[df$Metric %in% c("VPD_outside","SM_outside"),]
df_variable$Variable <- ifelse(df_variable$Metric=="VPD_outside","VPD","Soil moisture")
df_variable$Variable <- factor(df_variable$Variable,levels=c("VPD","Soil moisture"))
df_variable$Percent <- df_variable$Mean_fraction*100

g2 <- ggplot(df_variable,
             aes(x=Group,y=Percent,fill=Variable))+
  geom_boxplot(aes(group=interaction(Group,Variable)),
               position=position_dodge(width=0.65),
               width=0.5,color="black",linewidth=0.6,
               outlier.shape=21,outlier.color="black",
               outlier.fill="white",outlier.size=1.5)+
  scale_fill_manual(values=c("VPD"=VPD_color,
                             "Soil moisture"=SM_color),
                    guide=guide_legend(override.aes=list(color="black",
                                                        linewidth=0.6)))+
  scale_y_continuous(limits=c(0,60),breaks=seq(0,60,10),
                     expand=c(0,0))+
  my_theme+
  theme(text=element_text(size=12),
        axis.title=element_text(size=13),
        axis.text=element_text(size=11),
        axis.text.x=element_text(size=10),
        legend.position="top",
        legend.title=element_blank(),
        legend.text=element_text(size=11),
        plot.title=element_text(size=14))+
  labs(x="End-century scenario",y="Days outside range (%)")+
  ggtitle("b  Outside individual historical ranges")

# Combine the plots
g_all <- plot_grid(g1,g2,nrow=1,rel_widths=c(1,1.15))
print_g(g_all,"Future_outside_historical_range",9,4.5)
