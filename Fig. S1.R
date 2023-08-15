#fig. s1 top 10 phyla and classes-----
library(picante);	library(tidyverse)
#when load picante at the same time by default loading vegan
source('~/Desktop/molt plot/alpha/alpha.R')
source('~/Desktop/molt plot/theme_linweichuan.R')
group <- read.delim('~/Desktop/molt plot/groupALL.txt')
site <- unique(group$tissue)
site
data  <- NULL
for (i in site) {
	group1 <- subset(group,tissue==i)
	if (i== 'Foregut'|i== 'Midgut'|i== 'Hindgut') {
		otu<- read.delim('~/Desktop/molt plot/ASV_34848.txt',row.names = 1,header = T,check.names = F)
		otu <- otu[,colnames(otu) %in% group1$SampleID]
		ann <- read.delim('~/Desktop/molt plot/taxonomy_qc mc hc.txt',header = T,check.names = F)
	}
	if (i== 'Mixed gut') {
		otu<- read.delim('~/Desktop/molt plot/ASV_mix gut35399.txt',row.names = 1,header = T,check.names = F)
		otu <- otu[,colnames(otu) %in% group1$SampleID]
		ann <- read.delim('~/Desktop/molt plot/taxonomy_mix gut.txt',header = T,check.names = F)
	}	
	if (i=='Whole gut') {
		otu<- read.delim('~/Desktop/molt plot/ASV_whole gut.txt',row.names = 1,header = T,check.names = F)
		otu <- otu[,colnames(otu) %in% group1$SampleID]
		ann <- read.delim('~/Desktop/molt plot/taxonomy_whole gut.txt',header = T,check.names = F)
	}
	num <- mean(colSums(otu))/100
	otu$ASV_ID <- rownames(otu)
	otu1 <- merge(otu,ann[,c('ASV_ID','PC')],by = 'ASV_ID')
	otu1[which(otu1$ASV_ID == ''),'PC'] <- 'Others'
	rownames(otu1) <- otu1$ASV_ID
	otu2 <- otu1[-1]
	t<-table(otu2$PC)
	d <- otu2
	nd <- NULL
	for (j in 1:length(t)) {
		t1 <- dplyr::filter(otu2,PC==names(t)[j])
		t2 <- data.frame(t(colSums(t1[1:ncol(otu)-1])))
		rownames(t2) <- names(t)[j]
		nd <- rbind(nd,t2)
	}
	nd <- nd[order(rowSums(nd/num),decreasing = T),]
	colnames(nd) <- stringr::str_replace(colnames(nd), "X", "")
	a <- c("Betaproteobacteria","Gammaproteobacteria","Alphaproteobacteria","Unassigned","Firmicutes",         
				 "Bacteroidetes","Epsilonproteobacteria","Chlamydiae","Deltaproteobacteria","Fusobacteriia")
	nd10<- nd %>%filter(rownames(nd) %in% a)
	nd11<- rbind(nd10,num*100-colSums(nd10))
	nd11 <- nd11/num
	colSums(nd11)
	rownames(nd11)[nrow(nd11)] <- 'Others'
	nd11$Taxonomy <- rownames(nd11)
	a1  <-reshape2::melt(nd11,id.vars = 'Taxonomy',variable.name = 'SampleID',value.name = 'value')
	df  <- merge(a1,group,by='SampleID')
	####geom_bar
	te1 = aggregate(df$value,                      
									by  =list(df$Taxonomy,df$time,df$tissue),FUN='mean')
	te2 = aggregate(df$value,                      
									by  =list(df$Taxonomy,df$time,df$tissue),FUN='sd')
	te <- cbind(te1,te2[4])
	names(te) = c('Taxonomy','time','tissue','mean','sd')
	data <- rbind(data,te)
}
unique(data$Taxonomy)
colors <-c("#1F78B4","#B2DF8A","#FB9A99","#FF7F00",
					 "#6A3D9A","#FFFF99","#A6CEE3","#4C97FF",
					 '#E31A1C','#33A02C',"gray")
data$Taxonomy<-factor(data$Taxonomy,levels = c(a,'Others'))
data$time <- factor(data$time,levels = c('CK','0','3','6','9','12','24','48','72'))
data$tissue <- factor(data$tissue,levels = c('Foregut','Midgut','Hindgut','Mixed gut','Whole gut') )
bar <-ggplot(data,aes(time,mean,fill=Taxonomy))+
	geom_col(position = position_stack(reverse = F),width = 0.6)+
	facet_wrap(~tissue,ncol=3,nrow=2,scales = 'free_x')+
	labs(x=NULL,y="Relative abundance (%)",title=NULL)+
	scale_fill_manual(values = colors)+theme_bw()+
	scale_y_continuous(expand=c(0,2.5),limits = c(0,100.1),breaks = seq(0,100.1,25))+
	theme_linweichuan()+
	theme(strip.background = element_blank(),strip.text = element_text(size=8,hjust = 0.5),
				legend.title=element_blank(),legend.position = "right",legend.key.size=unit(0.4,'cm'),
				legend.text = element_text(size=8,hjust = 0,color = "black"))
bar
#ggsave('~/Desktop/bar.pdf',bar,width = 6.65,height = 4.25,dpi = 600,units = 'in')
