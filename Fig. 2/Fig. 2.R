#fig. 2a gut part ggtree------------------------------
library(vegan);library(ggtree);library(ggplot2)
otu <- read.delim("~/Desktop/molt plot/ASV_table all.txt",row.names = 1,header = T,check.names = F)
group<-read.delim("~/Desktop/molt plot/groupALL.txt",header = T)
colnames(group)[1]<-"label"
pclust<- vegdist(t(otu),method = 'bray')
tree <- hclust(dist(pclust), method = 'complete')
tree <- hclust(dist(pclust))%>%ggtree(layout="circular")+geom_text(aes(label=node),hjust=1)+geom_tiplab()
tree$data
tree <-full_join(tree$data,group,by="label")
tree$tissue <- factor(tree$tissue,levels = c('Foregut','Midgut','Hindgut','Mixed gut','Whole gut'))
tree$time <- factor(tree$time,levels = c('CK','0','3','6','9','12','24','48','72'))
color_mapping <- c('Foregut'='#8931B1','Midgut'='#0E70BC','Hindgut'='#D72B23',
									 'Mixed gut'="#E18727",'Whole gut'="#20854E",
									 'CK'='#E64B35','0'='#F39B7F','3'='#91D1C2','6'='#4DBBD5',
									 '9'='#00A087','12'='#8491B4','24'='#3C5488','48'="#984EA3",'72'='black')

ptree<-ggtree(tree,layout = "circular",branch.length = "none",
							linetype=1,size=0.3,ladderize = F)+#circular
	geom_tree(aes(color=time),size=0.3,show.legend=T)+
	geom_tippoint(aes(color=tissue),size=2,stroke=0)+
	scale_color_manual(values = color_mapping)+
	xlim_tree(0)+theme_tree()+
	guides(colour = guide_legend(order = 1,nrow = 2))+
	theme(legend.title = element_blank(), legend.position = 'none',
				legend.background = element_blank(),legend.key.size = unit(0.4,'in'),
				legend.text = element_text(size=8,color="black"),
				plot.margin = margin(           # adjust the margin
					t = 0.01,                       # top 
					r = 0.01,                       # right
					b = 0.01,                       # bottom
					l = 0.01,                       # left'
					unit = "in"                   # unit
				))
ptree
#ggsave('~/Desktop/ptree1.pdf',ptree,width = 3.12,height = 3.00,dpi = 600,units = 'in')

#fig. 2b top 10 phyla and classes-----
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
#write.csv(data,'~/Desktop/data.csv')
#ggsave('~/Desktop/bar.pdf',bar,width = 6.65,height = 4.25,dpi = 600,units = 'in')

#fig. 2c,d alpha diversity----
library(picante)
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
		otu1 <- otu[,colnames(otu) %in% group1$SampleID]
		tree <- read.tree('~/Desktop/molt plot/ASV_qc mc hc.tree')}
	if (i== 'Mixed gut') {
		otu<- read.delim('~/Desktop/molt plot/ASV_mix gut35399.txt',row.names = 1,header = T,check.names = F)
		otu1 <- otu[,colnames(otu) %in% group1$SampleID]
		tree <- read.tree('~/Desktop/molt plot/ASV_mix.tree')}	
	if (i=='Whole gut') {
		otu<- read.delim('~/Desktop/molt plot/ASV_whole gut.txt',row.names = 1,header = T,check.names = F)
		otu1 <- otu[,colnames(otu) %in% group1$SampleID]
		tree <- read.tree('~/Desktop/molt plot/ASV_whole gut.tree')}
	alpha_all <- alpha(t(otu1), tree, base = 2)
	alpha_all$SampleID <- rownames(alpha_all)
	alpha1 <- merge(alpha_all,group1,by='SampleID')
	data <- rbind(data,alpha1)
}
#write.csv(data,'~/Desktop/Fig. 2/result.csv',row.names = F)
####calculate significance
library(reshape2);library(agricolae);library(ggh4x)
data <- read.csv('~/Desktop/molt plot/alpha/result.csv')
dat1 <- melt(data,id.vars = c('SampleID', 'group', 'tissue','time'),
						 variable.name = 'name',value.name = 'value') 
site <- unique(data$tissue)
site
df <- NULL;res <- NULL
for (i in site) {
	dat2 <- subset(dat1,tissue==i)
	n <- unique(dat2$name) 
	for (j in n) {
		dat3 <- subset(dat2,name==j)
		comparison<-with(dat3,kruskal(value,time,group=TRUE, p.adj="BH"))
		te <- comparison$means
		te$time <- rownames(te)
		sig <- comparison$groups
		sig$time<-rownames(sig)
		sig <- sig[-1]
		te1 <- merge(te[,c(1:4,10)],sig,by='time')
		te1$tissue <- i
		te1$name <- j
		df <- rbind(df,te1)
	}
}
df1 <- dplyr::filter(df,name=='Shannon'|name=='Richness')
#write.csv(df1,'~/Desktop/palpha.csv',row.names = F)

#fig. 2e calculate community stability-----
pacman::p_load(reshape2,dplyr,ggplot2,stats,agricolae)
stability_onerep <- function(df, x){
	assertthat::assert_that(assertthat::has_name(df, x))
	assertthat::assert_that(is.numeric(df[[x]]))
	sync_var <- df[[x]]
	stability <- mean(sync_var, na.rm = TRUE)/stats::sd(sync_var, na.rm = TRUE)
	return(stability)
}
group <- read.delim('~/Desktop/molt plot/groupALL.txt')
site <- unique(group$tissue);site
data <- NULL
for (i in site) {
	group1 <- subset(group,tissue==i)
	if (i== 'Foregut'|i== 'Midgut'|i== 'Hindgut') {
		otu<- read.delim('~/Desktop/molt plot/ASV_34848.txt',row.names = 1,header = T,check.names = F)
	}
	if (i== 'Mixed gut') {
		otu<- read.delim('~/Desktop/molt plot/ASV_mix gut.txt',row.names = 1,header = T,check.names = F)
	}	
	if (i=='Whole gut') {
		otu<- read.delim('~/Desktop/molt plot/ASV_whole gut.txt',row.names = 1,header = T,check.names = F)
	}
	otu <- otu[,colnames(otu) %in% group1$SampleID]
	group2 <- group1[,c('SampleID','time')]
	otu1 <- t(otu)
	otu2 <- cbind(otu1,group2[,c('time')])
	otu3 <-melt(otu2,id.vars = 'time',variable.name = 'ASV',value.name = 'value')
	colnames(otu3) <- c('SampleID','ASV','abundance')
	otu3 <- merge(otu3,group2[,c('SampleID','time')],by ='SampleID')
	otu3$time <- factor(otu3$time,levels = c("CK","0 h","3 h","6 h","9 h","12 h","24 h"))#
	otu3$abundance <- as.numeric(factor(otu3$abundance))# to avoid NA
	####
	result <- NULL
	n <- unique(otu3$SampleID) 
	for (j in 1:length(n)) {
		otu4 <- dplyr::filter(otu3,SampleID==n[j])
		otu4 <- otu4[which(otu4[['abundance']] > 0),]
		output <- stability_onerep(otu4, 'abundance')
		s <- data.frame(SampleID=n[j],stability=output)
		s1 <- merge(s,group2,by = 'SampleID')
		result <- rbind(result,s1)
	}
	comparison<-with(result,kruskal(stability,time,group=TRUE, p.adj="BH"))
	cp <- comparison$means
	cp$time <- rownames(cp)
	sig <- comparison$groups
	sig$time <-rownames(sig)
	sig <- sig[-1]
	df3 <- merge(cp[,c(1:4,10)],sig,by = 'time')
	df3$tissue <- i
	data <- rbind(data,df3)
}
colnames(data) <- c('time','mean','rank','sd','N','sign','tissue')
#write.csv(data,'~/Desktop/molt plot/pstability.csv',row.names = F)

#fig. 2c-e plot alpha and stability----
#combine shannon, richness, and stability into palpha.csv
library(ggh4x)
df1 <- read.csv('~/Desktop/molt plot/alpha/palpha.csv')
colnames(df1) <- c('time','mean','rank','sd','N','sign','tissue','name')
df1$time <- factor(df1$time,levels  = c('CK','0','3','6','9','12','24','48','72'))
df1$name <- factor(df1$name,levels  = c('Richness','Shannon','Stability'),
									 labels = c('ASV richness','Shannon index','Community stability'))
df1$line <- 1
df1$tissue <- factor(df1$tissue,levels = unique(group$tissue))
pal <- c('#E64B35FF','#F39B7FFF','#91D1C2FF','#4DBBD5FF','#00A087FF','#8491B4FF','#3C5488FF',"#984EA3",'black')
p <- ggplot(df1, aes(x = time,y=mean,color=time))+
	geom_line(aes(group=line),linewidth=0.3) +
	geom_errorbar(aes(ymax = (mean+sd/1.05),ymin = (mean-sd/1.05)),
								width = 0.25,color='black',linewidth=0.3,
								position = position_dodge(0.8)) +
	geom_point(aes(color=time),shape=19,size=3,stroke=0) + # 21 is filled circle
	geom_text(aes(label=sign,y = (mean+sd)*1.15),color='black',size=2.8)+
	scale_color_manual(values = pal)+labs(x=NULL,y=NULL)+
	facet_grid(name~tissue,scales = "free",shrink = T)+#,space='free_x'
	facetted_pos_scales(
		y = list(name== 'ASV richness' ~
						 	scale_y_continuous(expand = c(0,20),breaks = seq(0,800,200),limits = c(0,800)),
						 name== 'Shannon index' ~
						 	scale_y_continuous(expand = c(0,0.15),breaks = seq(0,6,1.5),limits = c(0,6)),
						 name== 'Community stability' ~
						 	scale_y_continuous(expand = c(0,0.01),breaks = seq(0.1,0.3,0.05),limits = c(0.1,0.3))))+
	theme_linweichuan()+guides(color='none')+
	theme(strip.background = element_blank(),strip.text = element_text(size=8))
p
#ggsave('~/Desktop/molt plot/p.pdf',p,width = 8.50,height = 4.60,dpi = 600,units = 'in')

#Finally, these plots are arranged into Fig.2 in AI software and further modify legends.
