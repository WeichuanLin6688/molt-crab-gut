theme_linweichuan <- function(){
  pacman::p_load(ggplot2,showtext)
  theme_bw()+
  theme(panel.grid = element_blank(),panel.background = element_blank(),
        plot.title = element_text(color='black',size=8,hjust = 0.5),
        axis.title = element_text(color='black',size=8,hjust = 0.5), 
        axis.text = element_text(color='black',size=8,hjust = 0.5),
        axis.title.x=element_text(colour='black', size=8,hjust = 0.5), 
        axis.title.y=element_text(colour='black', size=8,hjust = 0.5),
        legend.text=element_text(colour='black', size=8, hjust = 0),
        legend.title =element_text(colour='black', size=8, hjust = 0),
        legend.background = element_rect(fill = 'transparent'), #c(0.9,0.5)
        axis.ticks.x = element_line(linewidth = 0.3),
        axis.ticks.y = element_line(linewidth = 0.3))
}
