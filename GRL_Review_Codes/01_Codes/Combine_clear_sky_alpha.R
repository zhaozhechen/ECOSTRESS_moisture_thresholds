# Date: 2026.9.13
# Combine the daily ESI comparison (a,b) and AMF alpha comparison (c,d).
# Alpha uncertainty remains in the results table, not in the figure.

Review_path <- "D:/Research/ECOSTRESS/Github repo/ECOSTRESS_moisture_thresholds/GRL_Review_Codes"
source(file.path(Review_path,"01_Codes/Clear_Sky_comparison.R"))
clear_sky_panels <- g_ls

# Reuse the alpha plotting code and calculated summary without repeating the fits.
summary_df <- read.csv(file.path(Review_path,"02_Results/Tables/AMF_alpha_netradiation_summary.csv"))
Site_ls <- c("US-A32","US-CF3")
alpha_code <- readLines(file.path(Review_path,"01_Codes/AMF_alpha_netradiation.R"))
plot_start <- grep("^# Direct comparison of alpha",alpha_code)
stopifnot(length(plot_start)==1)
eval(parse(text=alpha_code[plot_start:length(alpha_code)]))
alpha_panels <- g_ls

combined <- plot_grid(plotlist=c(clear_sky_panels,alpha_panels),ncol=2,
                      labels=c("a","b","c","d"),label_size=16,
                      align="hv",axis="lrbt")
print_g(combined,"Clear_Sky_comparison_with_alpha",8,8)
