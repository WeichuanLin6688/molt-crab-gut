#fig. S6----
#caclulate sub every sample network----
source("~/Desktop/molt github/Fig. S6/info.centrality.R");
source("~/Desktop/molt github/Fig. S6/nc.R")##
path <- c("~/Desktop/molt github/Fig. S6/all adj1/")
file_names<- list.files("~/Desktop/molt github/Fig. S6/all adj1/")
group <- read.delim('~/Desktop/molt github/molt plot/groupALL.txt')
sub_graph_stat <- NULL
for (i in 1:length(file_names)) {
	name<-gsub(".csv","",file_names[i])
	adj<-assign(name,read.csv(paste0(path,file_names[i]),row.names = 1,header = TRUE, stringsAsFactors = FALSE))
	if (name== 'Foregut'|name== 'Midgut'|name== 'Hindgut') {
		otu <- read.delim('~/Desktop/molt github/molt plot/ASV_34848.txt',row.names = 1,header = T,check.names = F)
	}
	if (name== 'Mixed gut') {
		otu <- read.delim('~/Desktop/molt github/molt plot/ASV_mix gut35399.txt',row.names = 1,header = T,check.names = F)
	}
	if (name=='Whole gut') {
		otu <- read.delim('~/Desktop/molt github/molt plot/ASV_whole gut.txt',row.names = 1,header = T,check.names = F)
	}
	head(otu)
	group1 <- subset(group,tissue==name)
	otu1 <- otu[,colnames(otu) %in% as.character(group1$SampleID)]
	n <- colnames(otu1)
	sub_graph <- list()
	for (j in 1:length(n)) {
		sample_i <- otu1[, colnames(otu1) %in% n[j], drop = FALSE]
		colnames(sample_i)[1] <- 'mean'
		select_node <- rownames(dplyr::filter(sample_i,mean>0))
		adj1 <- adj[rownames(adj) %in% select_node,colnames(adj) %in%select_node]
		remove_node <- sample(1:nrow(adj1), round(nrow(adj1)*0.5,0))
		adj_remove <- adj1[-remove_node,-remove_node]
		network_stability <- nc(adj_remove)
		fit <- try(graph_from_adjacency_matrix(as.matrix(abs(adj1)), mode = 'undirected', weight = T),T)
		if("try-error"%in% class(fit)){
			if (j > 56) {j <- 56}
			next
		}else{
			g <- graph_from_adjacency_matrix(as.matrix(abs(adj1)), mode = 'undirected', weighted = TRUE)
			g <- simplify(g)
			g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
			sub_graph[[j]] <- g
			node_num <- vcount(sub_graph[[j]])
			edge_num <- ecount(sub_graph[[j]])
			average_degree  <- mean(degree(sub_graph[[j]]))
			clustering_coefficient <-transitivity(sub_graph[[j]])
			####
			fit <- try(modularity(sub_graph[[j]], membership(cluster_fast_greedy(sub_graph[[j]]))),T) 
			if("try-error"%in% class(fit)){
				fc <- 0
			} else	{
				fc <-modularity(sub_graph[[j]], membership(cluster_fast_greedy(sub_graph[[j]])))
			}
			####
			fit <- try(max(data.frame(info.centrality.vertex(sub_graph[[j]]))),T) 
			if("try-error"%in% class(fit)){
				vulnerability <-  0
			} else	{
				vulnerability <-  max(data.frame(info.centrality.vertex(sub_graph[[j]])))
			}
			####
			fit <- try(average.path.length(sub_graph[[j]]),T) 
			if("try-error"%in% class(fit)){
				average_path_length = 0
			} else	{
				average_path_length <- average.path.length(sub_graph[[j]])
			}
			####
			sub_graph_stat1 <- data.frame(SampleID=n[j],
																		node=node_num,edge=edge_num,
																		average.degree=average_degree,
																		cc=clustering_coefficient,
																		apl=average_path_length,
																		modularity=fc,
																		network.stability=network_stability,
																		vul=vulnerability)
			sub_graph_stat <- rbind(sub_graph_stat,sub_graph_stat1)
			cat("正在处理任务", name, "，当前为：", n[j], "\n")
			print(sub_graph_stat1)
		}
	}
}
colnames(sub_graph_stat) <- c('SampleID','node num','edge num',
															'Average degree','Clustering coefficient',
															'Average path length','Modularity','Network stability',
															'Vulnerability')
#write.csv(sub_graph_stat,'~/Desktop/molt github/Fig. S6/subnetwork.csv',row.names = F)


##each network attributions significance----
library(agricolae);library(reshape2);library(ggplot2)
sub_graph_stat <- read.csv('~/Desktop/molt github/Fig. S6/subnetwork.csv')
group <- read.delim('~/Desktop/molt github/molt plot/groupALL.txt')
sub_graph_stat1 <- merge(sub_graph_stat,group[,c(1,3,4)],by = 'SampleID')
sub <- melt(sub_graph_stat1[4:11],id.vars = c("tissue",'time'),variable.name = 'name',value.name = 'value')
res <- NULL
n <- unique(sub$name)
n
t <- unique(group$tissue)
for (i in n) {
	sub1 <- subset(sub,name==i)	
	for (j in t) {
		sub2 <- subset(sub1,tissue==j)
		comparison<-with(sub2,kruskal(value,time,group=TRUE, p.adj="BH"))
		mark <- data.frame(comparison$groups)
		mark$time = row.names(mark)
		means <- comparison$means
		means$time = row.names(means)
		means$name =i 
		means$tissue =j
		sign <- merge(means,mark[,-1],by='time')
		res <- rbind(res,sign[,-c(5:9)])
	}
}
res 
head(res)
source('~/Desktop/molt github/molt plot/theme_linweichuan.R')
res$time <- factor(res$time,levels =c('CK','0','3','6','9','12','24','48','72'))
res$tissue <- factor(res$tissue,levels = c("Foregut","Midgut","Hindgut",'Mixed gut','Whole gut'))
res$name <-  gsub("\\.", " ",res$name)

sub$time <- factor(sub$time,levels =c('CK','0','3','6','9','12','24','48','72'))
sub$name <-  gsub("\\.", " ",sub$name)

n <- unique(res$tissue)
for (i in 1:length(n)) {
	df <- subset(res,tissue==n[i])
	df1 <- subset(sub,tissue==n[i])
	#df1$line <- 1
	assign(paste0('p',i),
				 ggplot(df1, aes(x = time, y = value,color=time,fill=time)) +
				 	geom_boxplot(outlier.shape = NA,color='black',position=position_dodge(0.8),
				 							 linewidth=0.3,fatten = 1.5,fill='white')+
				 	#geom_line(aes(group=line))+
				 	geom_jitter(aes(color=time),size=1.5,alpha=0.8,stroke=0,width=0.2)+
				 	scale_fill_manual(values = c('#E64B35FF','#F39B7FFF','#91D1C2FF','#4DBBD5FF',
				 															 '#00A087FF','#8491B4FF','#3C5488FF',"#984EA3",'black'))+
				 	scale_color_manual(values = c('#E64B35FF','#F39B7FFF','#91D1C2FF','#4DBBD5FF',
				 																'#00A087FF','#8491B4FF','#3C5488FF',"#984EA3",'black'))+
				 	geom_text(data=df,aes(x=time,y=round( Q75*1.2,2),label=groups),
				 						position=position_dodge(0.8),color="black",size = 2.8,fontface="plain")+
				 	labs(y = NULL,x=NULL)+
				 	facet_wrap(~name,scales = 'free',nrow = 6,ncol = 1)+
				 	theme_linweichuan()+
				 	guides(color=guide_legend(size=3,nrow=1,ncol=10),fill=guide_legend(size=3,nrow=1,ncol=10))+
				 	theme(legend.position = 'bottom',#legend.direction = 'horizontal',
				 				legend.title = element_blank(),
				 				strip.background = element_blank(),strip.text = element_text(size=8))
	)
}
library(ggpubr)

p <- ggarrange(p1,p2,p3,p4,p5,align = 'hv',nrow=1,ncol = 5,common.legend = T,legend = 'bottom')
p
#ggsave('~/Desktop/p.pdf',p,width = 10.00,height = 8.9,dpi = 600,units = 'in')

#Finally, these plots are arranged into fig.S6 in AI software and further modify legends.

