#fig. S3
#calculate disimmilarity---=
library(vegan)
source("~/Desktop/molt github/Fig. S3/pairwise.anosim.R")
group <- read.delim('~/Desktop/molt github/molt plot/groupALL.txt')
#读入物种数据
site <- unique(group$tissue)
data  <- NULL
for (i in site) {
	group1 <- subset(group,tissue==i)
	if (i== 'Foregut'|i== 'Midgut'|i== 'Hindgut') {
		otu<- read.delim('~/Desktop/molt github/molt plot/ASV_34848.txt',row.names = 1,header = T,check.names = F)
		otu1 <- otu[,colnames(otu) %in% group1$SampleID]}
	if (i== 'Mixed gut') {
		otu<- read.delim('~/Desktop/molt github/molt plot/ASV_mix gut35399.txt',row.names = 1,header = T,check.names = F)
		otu1 <- otu[,colnames(otu) %in% group1$SampleID]}	
	if (i=='Whole gut') {
		otu<- read.delim('~/Desktop/molt github/molt plot/ASV_whole gut.txt',row.names = 1,header = T,check.names = F)
		otu1 <- otu[,colnames(otu) %in% group1$SampleID]}
	#calculate disimmilarity
	g <- unique(group1$group)
	dis <- as.matrix(vegdist(t(otu1),method = 'bray'))
	comm<- matrix(0, length(g),length(g))
	colnames(comm) <- g;rownames(comm) <- g
	for (m in 1:nrow(comm)){
		group2 <- subset(group1,group== g[m])
		for (n in 1:ncol(comm)){
			group3 <- subset(group1,group== g[m])
			num <- mean(dis[rownames(dis) %in% group2$SampleID,colnames(dis) %in% group3$SampleID])
			comm[m,n] <- num
		}
	}
	diag(comm) <- 0
	#calculate significance
	otu2 <- data.frame(t(otu1))
	pairwise.anosim_result <- pairwise.anosim(otu2, group1$group, sim.method="bray", p.adjust.m= "BH")
	###export files
  write.csv(comm,paste0('~/Desktop/molt github/Fig. S3/',i,'_dis bray.csv'),row.names = T)
	write.csv(pairwise.anosim_result,paste0('~/Desktop/molt github/Fig. S3/',i,'_pairwise.anosim_result.csv'),row.names = F)
}

##plot dis----
library(corrplot);library(tidyverse)
group <- read.delim('~/Desktop/molt github/molt plot/groupALL.txt')
col <- colorRampPalette(c("#FFFFFF","#EE9988",'#BC3C29'))
site <- unique(group$tissue)
par(mfrow = c(2, 3))#Basic drawing patch
#Read the species data
for (i in site) {
	group1 <- subset(group,tissue==i)
	pattern <- print(i) # The characters you want to find
	files <- list.files(path ='~/Desktop/molt github/Fig. S3/',pattern = pattern)
	comm <- read.csv(paste0('~/Desktop/molt github/Fig. S3/',files[1]),row.names = 1)
	df <- read.csv(paste0('~/Desktop/molt github/Fig. S3/',files[2]))
	df1 <-df %>% separate(pairs,c('group1','group2'), " vs ")#divide into two columns
	ng <- unique(group1$group)
	head(df1)
	df2 <- df1[,c(1:2,5)]
	library(igraph)
	# Create a diagram from the adjacency table
	graph <- graph_from_data_frame(df2, directed = FALSE)
	# Convert the graph into an adjacent matrix
	adj <- as_adjacency_matrix(graph, attr = "p.adjusted", sparse = FALSE)
	# Check the adjacency matrix
	#print(adj_matrix)
	rownames(comm) <- unique(group1$time);colnames(comm) <- unique(group1$time)
	rownames(adj) <- unique(group1$time);colnames(adj) <- unique(group1$time)
	comm <- as.matrix(comm);adj <- as.matrix(adj)
	par(mar = c(2, 2, 2, 2)) 
	corrplot(comm,method = "square", type = "lower",
					 p.mat =adj,insig = 'blank', col.lim = c(0,1),
					 pch.cex = 0.6, pch.col = "black",
					 colorRampPalette(c("#4477AA", "#77AADD", "gray98", "#EE9988", "#BB4444"))(10),
					 number.cex = 0.8,cl.pos = "b",bg = "white",
					 cl.cex = 0.7,cl.offset = 0.3,cl.length = 6,cl.ratio = 0.14,
					 tl.offset = 0.6,tl.cex = 0.8,diag = F,tl.col = "black", #文本标签为黑色字体
					 tl.srt = 0 )# 上方文本标签为45度倾斜)
	cat('Now is in ',i,'\n')
}

#Finally, these plots are arranged into fig.S3 in AI software and further modify legends.
