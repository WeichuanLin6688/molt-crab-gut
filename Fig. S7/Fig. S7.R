#fig. S7----
#plot lm cor network attrubition----
###Mixed gut----
source('~/Desktop/molt github/molt plot/theme_linweichuan.R')
site1<- c('Foregut','Midgut','Hindgut')
sub_graph_stat <- read.csv('~/Desktop/molt github/Fig. S6/subnetwork.csv')
group <- read.delim('~/Desktop/molt github/molt plot/groupALL.txt')
data <- merge(sub_graph_stat,group[,c(1,3,4)],by = 'SampleID')
site2 <- colnames(data)[4:9]
mg <- NULL
for (i in site1) {
	data1 <- subset(data,tissue==i)
	data2 <- subset(data,tissue=='Mixed gut')
	for (j in site2) {
		data5 <- data.frame(x=data1[,colnames(data1) %in% j],y=data2[,colnames(data2) %in% j],
												time=data1$time,index=j,tissue=i)
		mg <- rbind(mg,data5)
	}
}
mg$time <- factor(mg$time,levels  = c('CK','0','3','6','9','12','24','48','72'))
mg$tissue <-factor(mg$tissue,levels = c('Foregut','Midgut','Hindgut'))
mg$index <-factor(mg$index,levels = site2,labels = gsub("\\.", " ",site2))
#Use the aggregation function to find the average according to the tissue*time*index variable
te <- aggregate(. ~ tissue*time*index, data = mg, mean)
te1 <- subset(te,!time=='CK')

p1 <- ggplot(te1, aes(y, x,color=tissue)) +
	geom_point(aes(color=tissue),size=2,stroke=0)+#shape=tissue,
	scale_color_manual(values = c('#8931B1','#0E70BC','#D72B23'))+#'#FD992D','#D72B23','#2D469C'
	geom_smooth(method = 'lm',#color='gray2',
							formula=y ~ x,linewidth= 0.5,se = F,show.legend = F)+
	stat_cor(method = "pearson",label.x.npc = 0.05,size=2.8,
					 r.digits = 2,p.digits = 3,r.accuracy = 0.01,p.accuracy = 0.001)+
	facet_wrap(~index,scales = 'free',shrink = T)+
	theme_linweichuan()+labs(x=NULL,y=NULL)+
	theme(strip.background = element_blank(),strip.text = element_text(size=8),
				legend.direction = 'horizontal', legend.position = 'none')
p1
#ggsave('~/Desktop/p1.pdf',p1,width = 8.00,height = 2.80,dpi = 300,units = 'in')

###whole gut----
site1<- c('Foregut','Midgut','Hindgut')
sub_graph_stat <- read.csv('~/Desktop/molt github/Fig. S6/subnetwork.csv')
group <- read.delim('~/Desktop/molt github/molt plot/groupALL.txt')
data <- merge(sub_graph_stat,group[,c(1,3,4)],by = 'SampleID')
site2 <- colnames(data)[c(4:9,11)]
mg <- NULL
for (i in site1) {
	data1 <- subset(data,tissue==i)
	df <-melt(data1[,c(4:9,11)])
	te1 <- aggregate(df$value, by  =list(df$time,df$variable),FUN='mean')
	data2 <- subset(data,tissue=='Whole gut')
	df <-melt(data2[,c(4:9,11)])
	te2 <- aggregate(df$value, by  =list(df$time,df$variable),FUN='mean')
	te <- merge(te1,te2,by = c('Group.1','Group.2'))
	te$tissue <- i
	colnames(te) <- c('time','index','x','y',"tissue")
	mg <- rbind(mg,te)
}
mg$time <- factor(mg$time,levels  = c('CK','0','3','6','9','12','24','48','72'))
mg$tissue <-factor(mg$tissue,levels = c('Foregut','Midgut','Hindgut'))
mg$index <-factor(mg$index,levels = site2,labels = gsub("\\.", " ",site2))
mg1 <- subset(mg,!time=='CK')
p2 <- ggplot(mg1, aes(y, x,color=tissue)) +
	geom_point(aes(color=tissue),size=2,stroke=0)+#shape=tissue,
	scale_color_manual(values = c('#8931B1','#0E70BC','#D72B23'))+#'#FD992D','#D72B23','#2D469C'
	geom_smooth(method = 'lm',#color='gray2',
							formula=y ~ x,linewidth= 0.5,se = F,show.legend = F)+
	stat_cor(method = "pearson",label.x.npc = 0.05,size=2.8,
					 r.digits = 2,p.digits = 3,r.accuracy = 0.01,p.accuracy = 0.001)+
	facet_wrap(~index,scales = 'free',shrink = T,nrow = 2)+
	theme_linweichuan()+labs(x=NULL,y=NULL)+#guides(linetype=guide_axis())
	theme(strip.background = element_blank(),strip.text = element_text(size=8),
				legend.direction = 'horizontal', legend.position = 'none')
p2

p3 <- ggarrange(p1,p2,nrow = 2,labels = c('A','B'),align = 'hv',
								font.label = list(size = 10,face = "plain"))
p3
#ggsave('~/Desktop/p3.pdf',p3,width = 6.18,height = 7.88,dpi = 600,units = 'in')

#Finally, these plots are arranged into fig.S7 in AI software and further modify legends.
