# This code includes functions for plotting

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

print_g <- function(g,title,w,h){
  pdf(paste0(Output_path,"/",title,".pdf"),
      width=w,height=h)
  print(g)
  dev.off()
  png(paste0(Output_path,"/",title,".png"),
      width=w,height=h,units = "in",
      res=600)
  print(g)
  dev.off()
}