# DESeq2------------------------------
library(DESeq2)
group <- read.delim('~/Desktop/molt github/molt plot/groupMOLT.txt')
site <- c("Foregut","Midgut","Hindgut")
site
data  <- NULL
#generate a compared DESeq2 data in each gut group
for (i in site) {
	group1 <- subset(group,tissue==i)
	otu<- read.delim('~/Desktop/molt github/molt plot/ASV_34848.txt',row.names = 1,header = T,check.names = F)
	ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_qc mc hc.txt')
	otu1 <- otu[,colnames(otu) %in% group1$SampleID]
	n <- unique(group1$time)
	n1 <- n[!n==0]
	for (j in n1) {
		group2 <- subset(group1,time==j|time==0)
		com <- otu1[,as.character(group2$SampleID)]
		com <- com[which(rowSums(com)>0),] 
		group_com <- group2[,c("SampleID",'time')] 
		group_com$time <- factor(group_com$time)
		dds <- DESeqDataSetFromMatrix(countData = com, colData = group_com, design = ~time)
		dds <- DESeq(dds, parallel = FALSE)	#parallel = TRUE 将启用多线程模式
		suppressMessages(dds)
		res <- results(dds, contrast = c('time',j,'0'), pAdjustMethod = 'BH', alpha = 0.05)
		summary(res)
		resdata <- as.data.frame((counts(dds)),by="row.names",sort=FALSE)
		as.data.frame(resdata) -> mut_exp
		deseq_res <- as.data.frame(res[order(res$padj), ])
		#输出
		deseq_res$ASV_ID <- rownames(deseq_res)
		deseq_res$tissue<- i
		deseq_res$time <- j
		##合并两列
		deseq_res <- merge(deseq_res,ann[,c("ASV_ID",'PC','Family','Genus')],by="ASV_ID")
		deseq_res=deseq_res[order(deseq_res$padj,decreasing = F),]
		#例如这里根据 |log2FC| >= 1 & FDR p-value < 0.05 定义“差异”
		deseq_res[which(deseq_res$padj %in% NA),'sig'] <- 'Unchanged'
		deseq_res[which(deseq_res$log2FoldChange >= 1 & deseq_res$padj < 0.05),'sig'] <- 'Enriched'
		deseq_res[which(deseq_res$log2FoldChange <= -1 & deseq_res$padj < 0.05),'sig'] <- 'Depleted'
		deseq_res[which(abs(deseq_res$log2FoldChange) < 1 | deseq_res$padj >= 0.05),'sig'] <- 'Unchanged'
		data <- rbind(data,deseq_res)
	}
}
#write.csv(data,"~/Desktop/molt github/Fig. S8/deseq_res.csv",row.names = F)

#fig. S8 GHI plot volcano------------------------------
library(ggplot2)
source('~/Desktop/molt github/molt plot/theme_linweichuan.R')
group <- read.delim('~/Desktop/molt github/molt plot/groupMOLT.txt')
plotdata <- read.csv("~/Desktop/molt github/Fig. S8/deseq_res.csv")
plotdata$sig <- factor(plotdata$sig,levels = c('Enriched','Unchanged','Depleted'))
plotdata$tissue <- factor(plotdata$tissue,levels = c('Foregut','Midgut','Hindgut'))
plotdata$time  <- factor(plotdata$time,levels =c('CK','3','6','9','12','24','48','72'))
site <- unique(group$tissue)
t <- unique(group$time)
ann <- NULL
for (i in 1:length(site)) {
	plotdata1 <- subset(plotdata,tissue == site[i])
	if (site[i]=='Foregut') {x=15;y=10^-6}
	if (site[i]=='Midgut')  {x=15;y=10^-12}
	if (site[i]=='Hindgut') {x=15;y=10^-20}
	for (j in 1:length(t)) {
		plotdata2 <- subset(plotdata1,time == t[j])
		if (nrow(plotdata2)==0) {
			res <- data.frame(E=0,
												U=nrow(subset(plotdata2, sig == 'Unchanged')),
												D=0,
												log2FoldChange=x,padj=y,
												tissue= site[i],time=t[j])
		}else{
			plotdata3 <- subset(plotdata2,!log2FoldChange==Inf&!log2FoldChange==-Inf)
			res <- data.frame(E=nrow(subset(plotdata2,sig == 'Enriched')),
												U=nrow(subset(plotdata2,sig == 'Unchanged')),
												D=nrow(subset(plotdata2,sig == 'Depleted')),
												log2FoldChange=x,padj=y,
												tissue= site[i],time=t[j])
			ann <- rbind(ann,res)
		}
	}
}
ann$label <- paste0('E = ', ann$E, '\nU = ', ann$U, '\nD = ', ann$D)
ann$tissue <- factor(ann$tissue,levels = c('Foregut','Midgut','Hindgut'))
ann$time  <- factor(ann$time,levels =c('CK','3','6','9','12','24','48','72'))
p <- ggplot(plotdata, aes(x =log2FoldChange, y =-log(padj,10))) +
	geom_vline(xintercept = c(-1, 1), color = 'black', linewidth = 0.3,linetype=2)+
	geom_hline(yintercept = -log(0.05, 10), color = 'black', linewidth= 0.3,linetype=2)+
	geom_point(aes(shape=sig,fill=sig,color = sig),size=2.5,stroke=0.3,alpha=0.6) +
	scale_shape_manual(values = c(24,1,25))+
	scale_fill_manual(values  = c("#E41A1C",'gray',"#377EB8")) +
	scale_color_manual(values  = c("#E41A1C",'gray',"#377EB8")) +
	labs(x=bquote(Log[2]~Foldchange),y=bquote(-Log[10]~(padj)),title =NULL)+ 
	geom_text(data=ann,aes(label=label),size=2.2)+
	facet_grid(tissue~time,scales = "free",shrink = T)+#,space='free_x'
	theme_linweichuan()+guides(fill=guide_legend(reverse = T),
														 shape=guide_legend(reverse = T),
														 color=guide_legend(reverse = T))+
	theme(legend.title = element_blank(),axis.title.x =  element_text(size=8),
				axis.title.y = element_text(size=8),legend.position = 'none',
				axis.text = element_text(size=8),
				strip.background = element_blank(),
				strip.text.y.right = element_text(size=8,angle = 0,hjust = 0),
				strip.text.x.top  = element_text(size=8,hjust = 0.5))
p
#ggsave('~/Desktop/vol.pdf',p,width =11.50 ,height =5.00 ,dpi = 600,units = 'in')
head(plotdata)

#fig. S8 ABC------------------------------
##Disappearing species ------------------------------
pacman::p_load(tidyverse,ggplot2,ggsci)
otu <- read.delim('~/Desktop/molt github/molt plot/ASV_34848.txt',row.names = 1,check.names = F)
group <- read.delim('~/Desktop/molt github/molt plot/groupMOLT.txt')
otu <- otu[,colnames(otu) %in% group$SampleID]
group$tissue <- factor(group$tissue,levels = c('Foregut','Midgut','Hindgut'))
n <- c('Foregut','Midgut','Hindgut')
ctr <- c('CK','3','6','9','12','24','48','72')
ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_qc mc hc.txt')
disappear <- NULL
for (i in 1:length(n)) {
	group1 <- filter(group,tissue==n[i])
	otu1 <- otu[,colnames(otu) %in% group1$SampleID]
	otu1 <-otu1[which(rowSums(otu1)>0),] 
	for (j in 1:length(ctr)) {
		control1 <- subset(group1,time== '0')
		control <- otu[,colnames(otu) %in% control1$SampleID]
		control <- control[which(rowSums(control)>0),] 
		c2 <- filter(group1,time== ctr[j])
		c3 <- otu[,colnames(otu) %in% c2$SampleID]
		c4 <- c3[which(rowSums(c3)>0),] 
		d1 <- setdiff(rownames(c4),rownames(control))
		d2 <- otu[rownames(otu) %in% as.character(d1),colnames(otu) %in% c2$SampleID]#找0h以外中不同于0h的物种
		d2 <- d2[which(rowSums(d2)>0),]
		cat(n[i],ctr[j],nrow(d2),'\n')
		d2$mean <-rowMeans(d2)/348.48
		d2$ASV_ID <- rownames(d2)
		d3 <- merge(d2,ann[,c("ASV_ID",'PC','Family')],by = 'ASV_ID')
		a1 <- d3[,c("ASV_ID",'mean','PC','Family')]
		te = aggregate(a1$mean,                      
									 by  =list(a1$PC),FUN='sum')
		colnames(te) <- c('PC','abundance')
		te$time <- ctr[j];te$tissue <- n[i]
		disappear <- rbind(disappear,te)
	}
	#write.csv(disappear,"~/Desktop/molt github/Fig. S8/disapper.csv",row.names = T)
}

color_mapping <- c("Betaproteobacteria" = "#1F78B4","Gammaproteobacteria" = "#B2DF8A", 
									 "Alphaproteobacteria"="#FB9A99","Firmicutes"="#FF7F00",
									 "Bacteroidetes"= "#6A3D9A","Actinobacteria"="#FFFF99",
									 "Chlamydiae"="#A6CEE3","Chlamydiae"="#4C97FF",
									 "Deltaproteobacteria"='#E31A1C','Verrucomicrobia'='#33A02C','Others'="black",
									 'Unassigned'='gray')

for (i in 1:length(n)) {
	df <- filter(disappear,tissue==n[i])
	df1 <- data.frame(PC='Acidobacteria',abundance=0,time='0',tissue=n[i])
	df2 <- rbind(df,df1)
	df2$time <- factor(df2$time,levels=c('CK','0','3','6','9','12','24','48','72'))
	assign(paste0('p',n[i]),
				 ggplot(df2,aes(time,abundance,fill=PC))+
				 	geom_col(position = position_stack(reverse = F),width = 0.6)+
				 	labs(x=NULL,y="Relative abundance (%)",title=NULL)+
				 	scale_fill_manual('Phylum/Class',values = color_mapping)+
				 	scale_y_continuous(expand = c(0,0.3),limits = c(0,15),breaks = seq(0,15,3))+
				 	theme_linweichuan()+theme(legend.key.size = unit(0.4,'cm'))+
				 	guides(fill = guide_legend(ncol = 1, nrow = 10, byrow = F))+
				 	annotate('text', label = n[i],                    
				 					 x = 1.5, y = 9.5/10*15,size = 2.8, parse = TRUE)
	)
}
head(df)

#fig. S8 EFG------------------------------
# Enriched species ------------------------------
pacman::p_load(tidyverse,ggplot2,ggsci)
otu <- read.delim('~/Desktop/molt github/molt plot/ASV_34848.txt',row.names = 1,check.names = F)
group <- read.delim('~/Desktop/molt github/molt plot/groupMOLT.txt')
otu <- otu[,colnames(otu) %in% group$SampleID]
ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_qc mc hc.txt',header = T,check.names = F)
n <- c('Foregut','Midgut','Hindgut')
ctr <- c('CK','3','6','9','12','24','48','72')
deseq <- read.csv('~/Desktop/molt github/Fig. S8/deseq_res.csv')
enriched <- NULL
for (i in 1:length(n)) {
	group1 <- subset(group,tissue==n[i])
	deseq1 <- subset(deseq,tissue==n[i])
	for (j in 1:length(ctr)) {
		deseq2 <- subset(deseq1,time==ctr[j]&sig=='Enriched')
		if (nrow(deseq2)==0) {
			te <- data.frame(abundance=0,PC='Alphaproteobacteria',
											 tissue= n[i],time=ctr[j])
			enriched <- rbind(enriched,te)
		}else{
			otu1 <- otu[as.character(deseq2$ASV_ID),group1$SampleID]
			otu1 <-otu1[which(rowSums(otu1)>0),] 
			otu1$ASV_ID <- rownames(otu1)
			otu2 <- merge(otu1,ann[,colnames(ann) %in% c('ASV_ID','PC')],by = 'ASV_ID')
			otu2[which(otu2$ASV_ID == ''),'PC'] <- 'Others'
			rownames(otu2) <- otu2$ASV_ID
			otu2 <- otu2[-1]
			t<-table(otu2$PC)
			nd <- NULL
			for (k in 1:length(t)) {
				t1 <- subset(otu2,PC==names(t)[k])
				t2 <- data.frame(t(colSums(t1[1:ncol(otu2)-1])))
				rownames(t2) <- names(t)[k]
				nd <- rbind(nd,t2)
			}
			nd <- nd[order(rowSums(nd),decreasing = T),]
			nd$PC <- rownames(nd)
			a1 <- melt(nd,id.vars = c('PC'),value.name = 'value',variable.name = 'SampleID')
			head(a1)
			a1$SampleID <- gsub('X', '',a1$SampleID)
			
			df  <- merge(a1,group1,by='SampleID')
			te = aggregate(df$value/348.48,                      
										 by  =list(df$PC),FUN='mean')
			colnames(te) <- c('PC','abundance')
			te$time <- ctr[j];te$tissue <- n[i]
			enriched <- rbind(enriched,te)
		}
		
	}
	print(n[i])
}

color_mapping <- c("Betaproteobacteria" = "#1F78B4","Gammaproteobacteria" = "#B2DF8A", 
									 "Alphaproteobacteria"="#FB9A99","Firmicutes"="#FF7F00",
									 "Bacteroidetes"= "#6A3D9A","Actinobacteria"="#FFFF99",
									 "Chlamydiae"="#A6CEE3","Chlamydiae"="#4C97FF",
									 "Deltaproteobacteria"='#E31A1C','Verrucomicrobia'='#33A02C','Others'="gray",
									 'Unassigned'='gray')

for (i in 1:length(n)) {
	df <- subset(enriched,tissue==n[i])
	df1 <- data.frame(PC='Alphaproteobacteria',abundance=0,time='0',tissue=n[i])
	df2 <- rbind(df,df1)
	df2$time <- factor(df2$time,levels = c('CK','0','3','6','9','12','24','48','72'))
	assign(paste0('p',n[i],'1'),
				 ggplot(df2,aes(time,abundance,fill=PC))+
				 	geom_col(position = position_stack(reverse = F),width = 0.6)+
				 	labs(x=NULL,y="Relative abundance (%)",title=NULL)+
				 	scale_fill_manual('Phylum/Class',values = color_mapping)+
				 	scale_y_continuous(expand = c(0,1),limits = c(0,30),breaks = seq(0,30,6))+
				 	theme_linweichuan()+theme(legend.key.size = unit(0.4,'cm'))+
				 	guides(fill = guide_legend(ncol = 1, nrow = 10, byrow = F))+
				 	annotate('text', label = n[i],                    
				 					 x = 1.5, y = 9.5/10*30,size = 2.8, parse = TRUE)
	)
}
head(df)
p2 <- ggpubr::ggarrange(pForegut,pMidgut,pHindgut,
												pForegut1,pMidgut1,pHindgut1,
												ncol = 3,nrow = 2,
												labels = c('A','B','C','D','E','F'),font.label = list(size = 10,face = "plain"),
												common.legend = T,align = 'hv',legend = 'right' )
p2

p3 <- ggarrange(p2,p,ncol = 1,nrow = 2, heights = c(0.6,0.6),
								labels = c('','GHI'),font.label = list(size = 10,face = "plain"))
p3#928*800
#ggsave('~/Desktop/pFig S8.pdf',p3,width = 9.40,height = 8.00,dpi = 600,units = 'in')

#Finally, these plots are arranged into fig. S8 in AI software and further modify legends.

