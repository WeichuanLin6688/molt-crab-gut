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

#fig. 2b comm Rstne-----
source('~/Desktop/molt plot/theme_linweichuan.R')
##foregut Rtsne-----
library(Rtsne)
otu<- read.delim('~/Desktop/molt plot/ASV_34848.txt',row.names = 1,header = T,check.names = F)
group<-read.delim("~/Desktop/molt plot/groupALL.txt",header = T)
otu1 <- otu[,colnames(otu1) %in% group$SampleID]
otu2 <- otu1[,colnames(otu1) %in% subset(group,tissue=='Foregut')$SampleID]#
otu3 <- data.frame(t(otu2))
bray_dis <- vegdist(otu3, method = 'bray')
groups1 <- filter(groups,tissue=='Foregut')
set.seed(123)
tsne_out <- Rtsne(otu3,
									pca=F,perplexity=(nrow(otu3)-1)/3,theta=0.45)#you can adjust theta value
tsne_res <- data.frame(tsne_out$Y,groups1)
colnames(tsne_res) <- c("tSNE1", "tSNE2",'SampleID', "tissue",'time','group')
#write.csv(tsne_res,'qc_rtsne.csv',row.names = F)
##foregut Rtsne-----
tsne_res <- read.csv('qc_rtsne.csv')
tsne_res$time <- factor(tsne_res$time,levels  = c('CK','0h','3h','6h','9h','12h','24h','48h','72h'),
												labels=c('CK','0','3','6','9','12','24','48','72'))
qc <-  ggplot(tsne_res,aes(tSNE1,tSNE2,color=time)) + 
	geom_point(aes(color = time),shape=19,size=2,stroke=0) +
	scale_color_manual(values = c('#E64B35FF','#F39B7FFF','#91D1C2FF',
																'#4DBBD5FF','#00A087FF','#8491B4FF',
																'#3C5488FF',"#984EA3",'black'))+
	geom_vline(xintercept = 0, color = 'gray', size = 0.5, linetype="dashed") +
	geom_hline(yintercept = 0, color = 'gray', size = 0.5, linetype="dashed") +
	labs(title='Foregut')+
	guides(color=guide_legend(size=2,nrow=3,ncol=3))+
	theme_linweichuan()+theme(legend.position = 'bottom',legend.title = element_blank(),
														legend.key.size = unit(0.2,'in'))
qc
##midgut Rtsne-----
tsne_res <- read.csv('mc_rtsne.csv')
tsne_res$time <- factor(tsne_res$time,levels  = c('CK','0h','3h','6h','9h','12h','24h','48h','72h'),
												labels=c('CK','0','3','6','9','12','24','48','72'))
mc <- ggplot(tsne_res,aes(tSNE1,tSNE2,color=time)) + 
	geom_point(aes(color = time),shape=19,size=2,stroke=0) +
	scale_color_manual(values = c('#E64B35FF','#F39B7FFF','#91D1C2FF',
																'#4DBBD5FF','#00A087FF','#8491B4FF',
																'#3C5488FF',"#984EA3",'black'))+
	geom_vline(xintercept = 0, color = 'gray', size = 0.5, linetype="dashed") +
	geom_hline(yintercept = 0, color = 'gray', size = 0.5, linetype="dashed") +
	labs(title='Midgut')+
	guides(color=guide_legend(size=2,nrow=3,ncol=3))+
	theme_linweichuan()+theme(legend.position = 'bottom',legend.title = element_blank(),
														legend.key.size = unit(0.2,'in'))
mc
##hindgut Rtsne-----
tsne_res <- read.csv('hc_rtsne.csv')
tsne_res$time <- factor(tsne_res$time,levels  = c('CK','0h','3h','6h','9h','12h','24h','48h','72h'),
												labels=c('CK','0','3','6','9','12','24','48','72'))
hc <- ggplot(tsne_res,aes(tSNE1,tSNE2,color=time)) + 
	geom_point(aes(color = time),shape=19,size=2,stroke=0) +
	scale_color_manual(values = c('#E64B35FF','#F39B7FFF','#91D1C2FF',
																'#4DBBD5FF','#00A087FF','#8491B4FF',
																'#3C5488FF',"#984EA3",'black'))+
	geom_vline(xintercept = 0, color = 'gray', size = 0.5, linetype="dashed") +
	geom_hline(yintercept = 0, color = 'gray', size = 0.5, linetype="dashed") +
	labs(title='Hindgut')+
	guides(color=guide_legend(size=2,nrow=3,ncol=3))+
	theme_linweichuan()+theme(legend.position = 'bottom',legend.title = element_blank(),
														legend.key.size = unit(0.2,'in'))
hc
##Mixed gut Rtsne-----
tsne_res <- read.csv('mix rtsne.csv')
tsne_res$time <- factor(tsne_res$time,levels  = c('CK','0','3','6','9','12','24','48','72'))
mix <- ggplot(tsne_res,aes(tSNE1,tSNE2,color=time)) + 
	geom_point(aes(color = time),shape=19,size=2,stroke=0) +
	scale_color_manual(values = c('#E64B35FF','#F39B7FFF','#91D1C2FF',
																'#4DBBD5FF','#00A087FF','#8491B4FF',
																'#3C5488FF',"#984EA3",'black'))+
	labs(title='Mixed gut')+
	geom_vline(xintercept = 0, color = 'gray', size = 0.5, linetype="dashed") +
	geom_hline(yintercept = 0, color = 'gray', size = 0.5, linetype="dashed") +
	guides(color=guide_legend(size=2,nrow=3,ncol=3))+
	theme_linweichuan()+theme(legend.position = 'bottom',legend.title = element_blank(),
														legend.key.size = unit(0.2,'in'))
mix
##whole gut Rtsne-----
tsne_res <- read.csv('pre_rtsne.csv')
tsne_res$time <- factor(tsne_res$time,levels  = c('CK','0','3','6','9','12','24','48','72'))
pre <- ggplot(tsne_res,aes(tSNE1,tSNE2,color=time)) + 
	geom_jitter(aes(color = time),shape=19,size=2,stroke=0,height = 2) +
	scale_color_manual(values = c('#E64B35FF','#F39B7FFF','#91D1C2FF',
																'#4DBBD5FF','#00A087FF','#8491B4FF',
																'#3C5488FF'))+
	geom_vline(xintercept = 0, color = 'gray', size = 0.5, linetype="dashed") +
	geom_hline(yintercept = 0, color = 'gray', size = 0.5, linetype="dashed") +
	labs(title='Whole gut')+
	guides(color=guide_legend(size=2,nrow=3,ncol=3))+
	theme_linweichuan()+theme(legend.position = 'bottom',legend.title = element_blank(),
														legend.key.size = unit(0.2,'in'))
pre
#
library(ggpubr)
# Extract the legend from the first plot
plot.list <- list(qc,mc,hc,mix,pre)
legend_grob <- ggplotGrob(plot.list[[1]])$grobs[[which(ggplotGrob(plot.list[[1]])$layout$name == "guide-box")]]
# Remove legends from all plots
plot.list.nolegend <- lapply(plot.list, function(x) x + theme(legend.position = "none"))
p <- ggarrange(plot.list.nolegend[[1]], plot.list.nolegend[[2]], plot.list.nolegend[[3]],
							 plot.list.nolegend[[4]], plot.list.nolegend[[5]], legend_grob,
							 ncol = 3, nrow = 2,	align = 'hv',#labels = c('A','B','C','D','E'),
							 font.label = list(size = 10, face = "plain"))
p
#ggsave('~/Desktop/fig. 2b.pdf',p,width = 5.60,height = 3.65,dpi = 600,units = 'in')


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

#Finally, these plots are arranged into fig.2 in AI software and further modify legends.
