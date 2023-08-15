#fig. 4 and fig. S11----
#we get fig. 4 and fig. S11 with the same code
###mc---
library(dplyr);library(ggplot2);library(vegan)
source('~/Desktop/molt plot/theme_linweichuan.R')
otu <- read.delim("~/Desktop/molt plot/ASV_34848.txt",check.names = F,row.names = 1,header = T)
group<-read.csv("~/Desktop/molt plot/groupMOLT mc.csv")
ann <- read.delim('~/Desktop/molt plot/taxonomy_qc mc hc.txt',row.names = 1)
#exclude weight
df <- read.delim('~/Desktop/randomforest/mc/Mc importance_otu0325.txt',check.names = F,row.names = 1)
df1 <- read.csv('~/Desktop/randomforest/mc/mc imp weight csv',check.names = F,row.names = 1)
n11 <- c('ASV_177','ASV_32','ASV_29','ASV_4743','ASV_14','ASV_1','ASV_444','ASV_233','ASV_132','ASV_7')

###hc---
otu <- read.delim("~/Desktop/molt plot/ASV_34848.txt",check.names = F,row.names = 1,header = T)
group<-read.csv("~/Desktop/randomforest/groupMOLT hc.csv")
ann <- read.delim('~/Desktop/molt plot/taxonomy_qc mc hc.txt',row.names = 1)
#exclude weight
df <- read.delim('~/Desktop/randomforest/hc/hc importance_otu0325.txt',check.names = F,row.names = 1)
df1 <- read.csv('~/Desktop/randomforest/hc/hc imp.weight.csv',check.names = F,row.names = 1)
n11 <- c('ASV_31','ASV_69','ASV_1','ASV_444','ASV_62')

#divide three category
diff <- setdiff(rownames(df),rownames(df1))
n <- 1:3
n
for (i in n) {
	if (i==1) {
		n1 <- n11
	}
	if (i==2) {
		n1 <- rownames(otu[!rownames(otu) %in% diff,])
	}
	if (i==3) {
		n1 <- diff[!diff %in% n11]
	}
	otu1 <- otu[rownames(otu) %in% n1,colnames(otu) %in% group$SampleID] 
	group2 <- group[,c('SampleID','time')]
	otu1 <- t(otu1)
	otu2 <- cbind(otu1,group2[,c('time')])
	otu3 <-melt(otu2,id.vars = 'time',variable.name = 'ASV',value.name = 'value')
	colnames(otu3) <- c('SampleID','ASV','abundance')
	otu3 <- merge(otu3,group2[,c('SampleID','time')],by ='SampleID')
	otu3$time <- factor(otu3$time,levels = c("CK","0.000001","3","6","9","12","24",'48','72'))#
	otu3$abundance <- as.numeric(factor(otu3$abundance))# to avoid NA
	####
	result <- NULL
	m <- unique(otu3$SampleID) 
	for (j in 1:length(m)) {
		otu4 <- dplyr::filter(otu3,SampleID==m[j])
		otu4 <- otu4[which(otu4[['abundance']] > 0),]
		output <- stability_onerep(otu4, 'abundance')
		s <- data.frame(SampleID=m[j],stability=output)
		s1 <- merge(s,group2,by = 'SampleID')
		result <- rbind(result,s1)
	}
	result$stability<- (result$stability-min(result$stability))/(max(result$stability)-min(result$stability))
	#calculate signficance
	comparison<-with(result,kruskal(stability,time,group=TRUE, p.adj="BH"))
	cp <- comparison$means
	cp$time <- rownames(cp)
	sig <- comparison$groups
	sig$time <-rownames(sig)
	sig <- sig[-1]
	df3 <- merge(cp[,c(1:4,10)],sig,by = 'time')
	#plot#
	df3$time <- factor(df3$time,levels = c('CK','0','3','6','9','12','24','48','72'))
	result$time <- factor(result$time,levels = c('CK','0','3','6','9','12','24','48','72'))
	assign(paste0('p',i),
				 ggplot(result,aes(x=time,y=stability))+
				 	geom_boxplot(outlier.shape = NA,color='black',position=position_dodge(0.8),
				 							 linewidth=0.3,fatten = 1.5,fill='white')+
				 	#geom_errorbar(aes(ymin=value-std/1.1,ymax=value+std/1.1),linewidth=0.3,width=0.3)+
				 	geom_jitter(aes(color=time),size=1.5,alpha=0.8,stroke=0,width = 0.2)+
				 	scale_color_manual(values = c('#E64B35FF','#F39B7FFF','#91D1C2FF','#4DBBD5FF',
				 																'#00A087FF','#8491B4FF','#3C5488FF',"#984EA3",'#000000'))+
				 	geom_text(data=df3,aes(label=groups,y=stability+std+0.1),size=2.8)+
				 	labs(x = NULL,y= 'Community stability')+
				 	theme_linweichuan()+
				 	theme(legend.position = 'none')
	)
	cat("now is",i,'\n')
	otu1 <- otu[rownames(otu) %in% n1,colnames(otu) %in% group$SampleID] 
	richness <- estimateR(t(otu1))[1, ]
	shannon <- diversity(t(otu1), index = 'shannon', base = 2)
	dat <- data.frame(cbind(richness,shannon))
	dat$SampleID <- rownames(dat)
	group2 <- group[,-c(3:5)]
	dat1 <- merge(dat,group2,by='SampleID')
	colnames(dat1)[2:3] <- c('ASV richness','Shannon index')
	for (d in colnames(dat1[2:3])) {
		dat2 <- dat1[,colnames(dat1)%in% c(d,'time')]
		colnames(dat2)[1] <- "mean"
		comparison<-with(dat2,kruskal(mean,time,group=TRUE, p.adj="BH"))
		cp <- comparison$means
		cp$time <- rownames(cp)
		sig <- comparison$groups
		sig$time <-rownames(sig)
		sig <- sig[-1]
		df3 <- merge(cp[,c(1:4,10)],sig,by = 'time')
		df3$time <- factor(df3$time,levels = c('CK','0','3','6','9','12','24','48','72'))
		dat2$time <- factor(dat2$time,levels = c('CK','0','3','6','9','12','24','48','72'))
		assign(paste0('d',d,i),
					 ggplot(dat2,aes(x=time,y=mean))+
					 	#stat_boxplot(geom = "errorbar", width = 0.15,linetype=2,position=position_dodge(0.8)) + 
					 	geom_boxplot(outlier.shape = NA,color='black',position=position_dodge(0.8),
					 							 linewidth=0.3,fatten = 1.5,fill='white')+
					 	geom_jitter(aes(color=time),size=1.5,alpha=0.8,stroke=0,width = 0.2)+
					 	scale_color_manual(values = c('#E64B35FF','#F39B7FFF','#91D1C2FF','#4DBBD5FF',
					 																'#00A087FF','#8491B4FF','#3C5488FF',"#984EA3",'#000000'))+
					 	geom_text(data=df3,aes(label=groups,y=mean+std+0.1),size=2.8)+
					 	labs(x = NULL,y= d)+
					 	theme_linweichuan()+
					 	theme(legend.position = 'none')
		)
	}
	##community composition
	otu2 <- merge(otu1,ann[,c('PC','Family')],by = 'row.names')
	otu3 <- otu2[,-c(1,57)]
	t<-unique(otu3$PC)
	nd <- NULL
	for (j in 1:length(t)) {
		t1 <- dplyr::filter(otu3,PC==t[j])
		t2 <- data.frame(t(colSums(t1[1:ncol(otu3)-1])))
		rownames(t2) <- t[j]
		nd <- rbind(nd,t2)
	}
	nd <- nd[order(rowSums(nd),decreasing = T),]
	rownames(nd)
	if ('Unassigned' %in% rownames(nd)) {
		a <- rownames(nd[!rownames(nd)%in% 'Unassigned',])
	}else{
		a <- rownames(nd)
	}
	
	if (nrow(nd)>=10) {
		nd10<- nd %>%filter(rownames(nd) %in% a[1:nrow(nd)])
	}else{
		nd10<- nd %>%filter(rownames(nd) %in% a[1:10])
	}
	##
	nd11<- rbind(nd10,colSums(nd)-colSums(nd10))
	if (sum(nd11[nrow(nd11),])>0) {
		rownames(nd11)[nrow(nd11)] <- 'Others'
	}else{
		nd11 <- nd11[-nrow(nd11),]
	}
	nd11$Taxonomy <- rownames(nd11)
	colnames(nd11) <- stringr::str_replace(colnames(nd11), "X", "")
	a1  <-reshape2::melt(nd11,id.vars = 'Taxonomy',variable.name = 'SampleID',value.name = 'value')
	df  <- merge(a1,group[,c('SampleID',"time")],by='SampleID')
	####geom_area
	te1 = aggregate(df$value/348.48,                      
									by  =list(df$Taxonomy,df$time),FUN='mean')
	te2 = aggregate(df$value/348.48,                      
									by  =list(df$Taxonomy,df$time),FUN='sd')
	te <- cbind(te1,te2[3])
	names(te) = c('Taxonomy','group','mean','sd')
	te$group <- factor(te$group,levels = c('CK',"0",'3','6','9','12','24','48','72'))
	color_mapping <- c("Betaproteobacteria" = "#1F78B4","Gammaproteobacteria" = "#B2DF8A", 
										 "Alphaproteobacteria"="#FB9A99","Firmicutes"="#FF7F00",
										 "Bacteroidetes"= "#6A3D9A","Actinobacteria"="#FFFF99",
										 "Chlamydiae"="#A6CEE3","Deltaproteobacteria"="#4C97FF",
										 "Fusobacteriia"='#E31A1C','Verrucomicrobia'='#33A02C','Others'="gray",
										 'Unassigned'="gray",'Epsilonproteobacteria'="#FFFF99")
	te$Taxonomy<-factor(te$Taxonomy,levels = rownames(nd11))
	assign(paste0('b',i),
				 ggplot(te,aes(group,mean,fill=Taxonomy))+
				 	geom_col(position = position_stack(reverse = F),width = 0.6)+
				 	scale_fill_manual(values = color_mapping)+
				 	labs(x=NULL,y="Relative abundance (%)",title=NULL)+
				 	scale_y_continuous(limits = c(0,100.1),breaks = seq(0,100.1,25))+
				 	theme_linweichuan()+
				 	theme(strip.background = element_blank(),strip.text.y.right = element_text(angle = 0,size=8),
				 				legend.position = 'right',legend.key.size = unit(0.4,'cm'))+
				 	guides(fill = guide_legend(ncol = 1,reverse = FALSE))
	)
}

te2 = aggregate(subset(te,!group=='CK')$mean,                      
								by  =list(subset(te,!group=='CK')$group),FUN='sum')
range(te2$x)

b4 <- b3+scale_y_continuous(limits = c(0,20),breaks = seq(0,20,5))+
	theme(legend.position = 'none',axis.ticks.x = element_blank(),
				axis.title.y = element_blank(),
				axis.text.x = element_blank())
b4
library(grid)
b5<-ggplotGrob(b4)
b6 <-b3+annotation_custom(b5,xmin=1,xmax=9,ymin=25,ymax=95)
b6
#pb1 <- b2+	theme(legend.position = 'right')
#ggsave('pd1.pdf',pb1,width = 4.77,height = 2.77,dpi = 600,units = 'in')
p <- ggarrange(p1,p2,p3,b1,b2,b6,
							 `dASV richness1`,`dASV richness2`,`dASV richness3`,
							 nrow = 3,ncol=3,legend = 'none',align = 'hv',
							 labels = c('A','D','G','B','E','H','C','F','I','D'),
							 font.label = list(size = 10, face = "plain"))
p
#ggsave('~/Desktop/p.pdf',p,width = 7.80,height = 5.45,dpi = 600,units = 'in')

#Finally, these plots are arranged into fig.4 and fig. S11 in AI software and further modify legends.




