#fig S4----
library(ggplot2);library(ggpmisc)
##Mixed gut----
site1<- c('Foregut','Midgut','Hindgut')
site2 <- c('Shannon','Richness','Stability')
data <- read.csv('~/Desktop/molt github/Fig. S4/alpha stability.csv')
####Mixed gut
te <- NULL
for (i in site1) {
	data1 <- subset(data,tissue==i)
	data2 <- subset(data,tissue=='Mixed gut')
	for (j in site2) {
		data5 <- data.frame(x=data1[,colnames(data1) %in% j],y=data2[,colnames(data2) %in% j],
												time=data1$time,index=j,tissue=i)
		te <- rbind(te,data5)
	}
}
dat1 <- aggregate(. ~ tissue*time*index, data = te, mean)
dat1$time <- factor(dat1$time,levels  = c('CK','0','3','6','9','12','24','48','72'))
dat1$tissue <-factor(dat1$tissue,levels = c('Foregut','Midgut','Hindgut'))
dat1$index <-factor(dat1$index,levels = c('Richness','Shannon','Stability'),
										labels = c('ASV richness','Shannon index','Community stability'))
dat1 <- subset(dat1,!time=='CK')
p1 <- ggplot(dat1, aes(y, x,color=tissue)) +
	geom_point(aes(color=tissue),size=2,stroke=0)+#shape=tissue,
	scale_color_manual(values = c('#8931B1','#0E70BC','#D72B23'))+#'#FD992D','#D72B23','#2D469C'
	geom_smooth(method = 'lm',#color='gray2',
							formula=y ~ x,linewidth= 0.5,se = F,show.legend = F)+
	stat_poly_eq(
		aes(label = paste(after_stat(adj.rr.label), ..p.value.label.., sep = "*\", \"*")),#..eq.label..,
		formula = y ~ x, parse = TRUE,  small.r = F,small.p = T,
		vstep = 0.1, label.x.npc = 0.05,#label.y.npc = 0.9,
		rr.digits = 2,p.digits = 3,size = 2.8) +
	facet_wrap(~index,scales = 'free',shrink = T)+
	#facet_grid(index~tissue,scales = 'free')
	theme_linweichuan()+labs(x=NULL,y=NULL)+#guides(linetype=guide_axis())
	theme(strip.background = element_blank(),strip.text = element_text(size=8),
				legend.direction = 'horizontal', legend.position = 'none')
p1

##Whole gut----
library(reshape2)
site1<- c('Foregut','Midgut','Hindgut')
site2 <- c('Shannon','Richness','Stability','tissue','time')
dat <- read.csv('~/Desktop/molt github/Fig. S4/alpha stability.csv')
data <- dat[,colnames(dat) %in% site2]
dat2 <- NULL
for (i in site1) {
	data1 <- subset(data,tissue==i)
	df <-melt(data1)
	te1 <- aggregate(df$value, by  =list(df$time,df$variable),FUN='mean')
	data2 <- subset(data,tissue=='Whole gut')
	df <-melt(data2)
	te2 <- aggregate(df$value, by  =list(df$time,df$variable),FUN='mean')
	mg <- merge(te1,te2,by = c('Group.1','Group.2'))
	mg$tissue <- i
	colnames(mg) <- c('time','index','x','y',"tissue")
	dat2 <- rbind(dat2,mg)
}
dat2$time <- factor(dat2$time,levels  = c('CK','0','3','6','9','12','24'))
dat2$tissue <-factor(dat2$tissue,levels = c('Foregut','Midgut','Hindgut'))
dat2$index <-factor(dat2$index,levels = c('Richness','Shannon','Stability'),
										labels = c('ASV richness','Shannon index','Community stability'))
dat2 <- subset(dat2,!time=='CK')
p2 <- ggplot(dat2, aes(y, x,color=tissue)) +
	geom_point(aes(color=tissue),size=2,stroke=0)+
	scale_color_manual(values = c('#8931B1','#0E70BC','#D72B23'))+#'#FD992D','#D72B23','#2D469C'
	geom_smooth(method = 'lm',#color='gray2',
							formula=y ~ x,linewidth= 0.5,se = F,show.legend = F)+
	stat_poly_eq(
		aes(label = paste(..adj.rr.label.., ..p.value.label.., sep = "*\", \"*")),#..eq.label..,
		formula = y ~ x, parse = TRUE,  small.r = F,small.p = T,
		vstep = 0.1, #label.x.npc = 0.05,#label.y.npc = 0.9,
		rr.digits = 2,p.digits = 3,size = 2.8) +
	facet_wrap(~index,scales = 'free',shrink = T)+
	theme_linweichuan()+labs(x=NULL,y=NULL)+
	theme(strip.background = element_blank(),strip.text = element_text(size=8),
				legend.direction = 'horizontal', legend.position = 'none')
p2
#
p3 <- ggarrange(p1,p2,nrow = 2,labels = c('A','B'),align = 'hv',legend = 'bottom',common.legend = T,
								font.label = list(size = 10,face = "plain"))
p3 
#ggsave('~/Desktop/palpha lm exCK.pdf',p3,width = 5.68,height = 3.88,dpi = 600,units = 'in')

#Finally, these plots are arranged into fig.S4 in AI software and further modify legends.
