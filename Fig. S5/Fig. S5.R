#fig. S5----
#construct network----
library(igraph);library(WGCNA);library(tidyverse);library(psych);library(multtest)
setwd('~/Desktop/molt github/Fig. S5/')
if(!dir.exists("all adj1")){dir.create("all adj1")}
group <- read.delim('~/Desktop/molt github/molt plot/groupALL.txt')
site <- unique(group$tissue)
#generate a adj matrix in each gut group
for (i in site) {
	group1 <- subset(group,tissue==i)
	if (i== 'Foregut'|i== 'Midgut'|i== 'Hindgut') {
		df <- read.delim('~/Desktop/molt github/molt plot/ASV_34848.txt',row.names = 1,header = T,check.names = F)
		ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_qc mc hc.txt')
	}
	if (i== 'Mixed gut') {
		df <- read.delim('~/Desktop/molt github/molt plot/ASV_mix gut35399.txt',row.names = 1,header = T,check.names = F)
		ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_mix gut.txt')
	}
	if (i=='Whole gut') {
		df <- read.delim('~/Desktop/molt github/molt plot/ASV_whole gut.txt',row.names = 1,header = T,check.names = F)
		ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_whole gut.txt')
	}
	###
	otu <- df[,colnames(df) %in% group1$SampleID]
	otu <- otu[which(rowSums(otu)>0),]
	myfun<-function(x){sum(x>0.001)}
	otu<-otu[as.logical(apply(otu/colSums(otu),1,myfun)),]
	occor<-WGCNA::corAndPvalue(t(otu),method = 'pearson')
	mtadj<-multtest::mt.rawp2adjp(unlist(occor$p),proc='BH')
	adpcor<-mtadj$adjp[order(mtadj$index),2]
	occor.p<-matrix(adpcor,dim(t(otu)/colSums(otu))[2])
	## R value
	occor.r<-occor$cor
	diag(occor.r) <- 0
	occor.r[occor.p>0.001|abs(occor.r)< 0.75] = 0  ##|abs(occor.r)< 0.6
	occor.r[is.na(occor.r)]=0
	diag(occor.r) <- 0
	#write.csv(occor.r,file=paste("~/Desktop/molt github/Fig. S5/all adj1/",i,".csv",sep=""),row.names = T)
	print(i)
}


##plot each netwrok2----
setwd('~/Desktop/molt github/Fig. S5/all adj1/');
library(igraph);library(tidygraph);library(ggraph)
path <- c("~/Desktop/molt github/Fig. S5/all adj1/")
group <- read.delim('~/Desktop/molt github/molt plot/groupALL.txt')
col_g <- "#999999"
file_names <- c("Foregut.csv","Midgut.csv","Hindgut.csv", "Mixed gut.csv","Whole gut.csv")
plots_list <- list()
for(i in 1:length(file_names)){
	name<-gsub(".csv","",file_names[i])
	adj<-assign(name,read.csv(paste0(path,file_names[i]),row.names = 1,header = TRUE, stringsAsFactors = FALSE))
	if (name== 'Foregut'|name== 'Midgut'|name== 'Hindgut') {
		otu <- read.delim('~/Desktop/molt github/molt plot/ASV_34848.txt',row.names = 1,header = T,check.names = F)
		ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_qc mc hc.txt')
	}
	if (name== 'Mixed gut') {
		otu <- read.delim('~/Desktop/molt github/molt plot/ASV_mix gut35399.txt',row.names = 1,header = T,check.names = F)
		ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_mix gut.txt')
	}
	if (name=='Whole gut') {
		otu <- read.delim('~/Desktop/molt github/molt plot/ASV_whole gut.txt',row.names = 1,header = T,check.names = F)
		ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_whole gut.txt')
	}
	head(otu)
	group1 <- subset(group,tissue==name)
	otu1 <- otu[,colnames(otu) %in% as.character(group1$SampleID)]
	n <- unique(group1$group)
	for (j in 1:length(n)) {
		group2 <- subset(group1,group==n[j])
		sample_i <- otu1[, colnames(otu1) %in% group2$SampleID, drop = FALSE]
		select_node <- rownames(sample_i[which(rowSums(sample_i)>0),])
		adj1 <- adj[rownames(adj) %in% select_node,colnames(adj) %in%select_node]
		set.seed(123)
		g <- graph_from_adjacency_matrix(as.matrix(abs(adj1)), mode = 'undirected', weight = T)
		g <- simplify(g)
		g <- delete_vertices(g, names(degree(g)[degree(g) == 0])) 
		g1 <- g
		####
		edge <- data.frame(as_edgelist(g1)) 
		edge_list <- data.frame(
			source = edge[[1]],
			target = edge[[2]])
		####
		node_list <- data.frame(
			ASV_ID = names(V(g1)),
			label = names(V(g1))
		)
		node_list1 <- merge(node_list,ann[,c('ASV_ID','PC','Family')],by = 'ASV_ID')
		ifelse(nrow(node_list1) < length(unique(node_list$ASV_ID)), 
					 dif <- setdiff(as.character(unique(node_list$ASV_ID)),
					 							 as.character(node_list1$ASV_ID)),
					 dif <- 0
		)
		if (length(dif)>0){
			for (k in 1:length(dif)) {
				node_list1 <- rbind(node_list1,
														data.frame(ASV_ID= dif[k],label= dif[k],
																			 PC='Others',Family='Others'))}
		}else{
			node_list1 <- node_list1
		}
		color_mapping <- c("Betaproteobacteria" = "#FB9A99","Gammaproteobacteria" = "#B2DF8A", 
											 "Alphaproteobacteria"= "#1F78B4","Firmicutes"="#FF7F00",
											 "Bacteroidetes"= "#6A3D9A","Actinobacteria"="#FFFF99",
											 "Fusobacteria"="#A6CEE3","Chlamydiae"="#4C97FF",
											 "Deltaproteobacteria"='#E31A1C','Tenericutes'='#33A02C','Unassigned'="gray",
											 'Candidatus_Saccharibacteria'="#E55C6D",'Epsilonproteobacteria'="#8DA0CB")
		g2 <- tbl_graph(nodes = node_list1, edges = edge_list,directed = F)
		g2 <- delete_vertices(g2, names(degree(g2)[degree(g2) == 0])) 
		g2 <- simplify(g2)
		V(g2)$size <- degree(g2)*0.2
		p <- ggraph(g2,layout = 'with_fr')+#layout = 'linear', circular = TRUE) +
			#geom_edge_fan(edge_color = "#8DA0CB",edge_width = 0.2) + 
			geom_node_point(aes(size=size,color=PC), shape = 19)+
			scale_color_manual(values = color_mapping) +
			scale_size(range =c(0.1,2))+#,breaks=seq(1,25,5),limits=c(1,25)
			labs(title = paste0('Nodes = ',length(node_list1$ASV_ID),', ',
													'edges = ',nrow(data.frame(as_edgelist(g2)))),
					 x=unique(subset(group1,group==n[j])$time))+
			theme_graph(title_family = "sans")+#title_family = "Helvetica"
			theme(legend.position = 'none',
						axis.title.x = element_text(size=8,hjust = 0),
						plot.title = element_text(face ="plain",size=6,hjust = 0.5),
						plot.margin = margin(           # 调整边距
							t = 0.2,                       # 上边距
							r = 0.1,                       # 右边距
							b = 0.1,                       # 下边距
							l = 0.1,                       # 左边距
							unit = "in"                   # 单位（例如，“pt”，“cm”，“inch”等）
						))
		p
		plots_list[[n[j]]] <- p
		print(n[j])
	}
}
library(gridExtra)
p1 <- ggarrange(plotlist = plots_list, ncol = 9,nrow=5,common.legend = T,legend = 'right')
p1
#ggsave('~/Desktop/p1.pdf',p1,width = 12.00,height = 8.10,dpi = 600,units = 'in')


#Finally, these plots are arranged into fig.S5 in AI software and further modify legends.



