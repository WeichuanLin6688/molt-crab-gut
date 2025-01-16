#fig. 2----
#fig. 2A----
##calculate the metabolism function significance-----
library(dplyr);library(agricolae);library(reshape2);
source('~/Desktop/molt github/molt plot/theme_linweichuan.R')
ann <- read.delim('~/Desktop/molt github/Fig. 2/KO1-5.txt',header = T)#Read the KEGG classification database
m <- filter(ann,PathwayL1=='Metabolism')
group <- read.delim('~/Desktop/molt github/molt plot/groupMOLT.txt')
site <- c('Foregut','Midgut','Hindgut')
data <- NULL;result <- NULL
for (i in site) {
	group1 <- subset(group,tissue==i)
  df <-  read.delim('~/Desktop/molt github/Fig. 2/KEGG.Pathway.raw.txt',check.names = F,row.names = 1)
	df1 <- df[,colnames(df) %in% group1$SampleID]
	df21 <- df1[unique(m$Pathway),]
	df21 <-df21[complete.cases(df21),]
	df2 <- df21
	colSums(df2)
	df3 <- melt(t(df2),value.name = 'value')
	colnames(df3) <- c('SampleID','name','value')
	n <- unique(df3$name)
	df3$SampleID <- gsub('X', '',df3$SampleID)
	df4 <- merge(df3,group,by = 'SampleID')
	colnames(df4) <- c('SampleID','name','value','group','tissue','time')
	df5 <- df4
	for (j in 1:length(n)) {
		df6 <- filter(df5, name==n[j])
		if (sum(df6$value)==0) {
			j=j+1
			next
		}else{
			comparison<-with(df6,kruskal(value,time,group=TRUE, p.adj="BH"))
			c <- comparison$means
			c
			p <- comparison$groups
			p
			if (length(unique(p$groups))==1) {
				j=j+1
				next
			}else{
				c$time <- rownames(c)
				c$funcs <-n[j]
				c$tissue <- i
				sig <- comparison$groups
				sig$time<-rownames(sig)
				sig <- sig[-1]
				names(sig) = c('sign','time')
				df7 <- merge(c,sig,by = 'time')
				result <- rbind(result,df7)
			}
		}
	}
	print(i)
}
##add over anonate----
a <- result
ko1 <- read.delim('~/Desktop/molt github/Fig. 2/KO1-5.txt')
plotdata <- NULL
for (i in 1:nrow(a)) {
	a1 <- a[i,]
	name <- unique(filter(ko1,Pathway==a1$funcs)$PathwayL2)
	a1$pathway <- name
	plotdata <- rbind(plotdata,a1)
	print(i)
}
unique(plotdata$pathway)
res <- filter(plotdata,!pathway=="Not included in regular maps")
res1 <- res
head(res1)
res1$time <- factor(res1$time,levels  = c('CK','0','3','6','9','12','24','48','72'))
res1$funcs <-gsub("\\[.*\\]","",res1$funcs)#delete [] context
result <- NULL
t <- c('Foregut','Midgut','Hindgut')
for (i in 1:length(t)) {
	res2 <- filter(res1,tissue==t[i])
	f <- unique(res1$funcs)
	for (j in 1:length(f)) {
		res3 <- filter(res2,funcs==f[j])
		if (length(unique(res3$sign))==1) {
			j=j+1
			res3$foldchange <- 1
			result <- rbind(result,res3)
			#next
		}else{
			res3$foldchange<- res3$value/min(filter(res3,!value==0)$value)#filter(res3,time=='0')$value#min(res3$value)
			result <- rbind(result,res3)
		}
	}
}
#write.csv(result,'~/Desktop/molt github/Fig. 2/plot_metabolism.csv',row.names = F)

##plot picrust----
result <- read.csv('~/Desktop/molt github/Fig. 2/plot_metabolism.csv')
log(max(filter(result,!foldchange==Inf)$foldchange))
result$tissue <- factor(result$tissue,levels = c('Foregut','Midgut','Hindgut'))
result1 <- result[complete.cases(result),]
result1$pathway <- factor(result1$pathway,levels=rev(sort(unique(result1$pathway))))
pal <- c("#59BD92","#EA963A","#2660A6","#E55C6D","#F7CB3C","#8D4093","#1F78B4",'#E31A1C',"#000000",'#AF5930',"#8DA0CB")
pal <- pal[1:length(unique(result1$pathway))]
#forcats::fct_reorder(funcs,pathway) #order funcs by pathway
p <- ggplot(result1,aes(y=time,x=forcats::fct_reorder2(funcs,funcs,pathway)))+
	geom_tile(color='transparent',fill='transparent',linewidth=0.3)+
	scale_size(range = c(0.5,3.5),breaks = c(0.5,1,2,4,6,8))+
	facet_wrap(~tissue,nrow = 5,strip.position = 'top')+
	geom_point(aes(size=log2(foldchange),color=pathway),alpha=1,stroke=0)+#4 cross
	scale_color_manual(values = rev(pal))+
	scale_y_discrete(limits=c('72','48','24','12','9','6','3','0','CK'))+
	#scale_x_reverse()+
	labs(x='Metabolism',y=NULL)+
	guides(size=guide_legend(title=NULL,#"Log2(foldchange)",
													 override.aes = list(nrow=3,ncol=2),order = 1),
				 color=guide_legend(title=NULL,reverse=T,override.aes = list(nrow=2,ncol=5),order = 2))+
	theme(strip.background = element_blank(),strip.text = element_blank(),
				panel.grid.major = element_line(linewidth=0.2,color='gray'),
				#panel.border = element_rect(linewidth=0.3,color='black',fill='transparent'),
				panel.border = element_blank(),
				panel.background  = element_blank(),
				axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
				axis.title.x = element_text(angle=0,hjust=0.5,vjust = 0.5,size=10,color = 'black'),
				axis.title.y = element_text(angle=90,hjust=0.5,vjust = 0.5,size=10,color = 'black'),
				axis.text.x = element_blank(),legend.spacing = unit(0, "cm"),
				#axis.text.x = element_text(angle=75,hjust=0.98,vjust = 1,size=9,color = 'black'),
				axis.text.y = element_text(angle=0,hjust=0.5,vjust = 0.5,size=8,color = 'black'),
				legend.text = element_text(size=6,color = 'black'),legend.title= element_text(size=6),
				legend.key.size = unit(0.2,'cm'),legend.background = element_blank(),
				legend.position = 'bottom')+
	annotate("rect", xmin = 0.0, xmax = 9.5,  ymin=-Inf, ymax=Inf,fill="#59BD92", alpha = 0.15)+
	annotate("rect", xmin = 9.5, xmax = 29.5, ymin=-Inf, ymax=Inf,fill="#EA963A", alpha = 0.15)+
	annotate("rect", xmin = 29.5, xmax =44.5, ymin=-Inf, ymax=Inf,fill="#2660A6", alpha = 0.15)+
	annotate("rect", xmin = 44.5, xmax = 52.5, ymin=-Inf, ymax=Inf,fill="#E55C6D", alpha = 0.15)+
	annotate("rect", xmin = 52.5, xmax = 63.5, ymin=-Inf, ymax=Inf,fill="#F7CB3C", alpha = 0.15)+
	annotate("rect", xmin = 63.5, xmax = 79.5, ymin=-Inf, ymax=Inf,fill="#8D4093", alpha = 0.15)+
	annotate("rect", xmin = 79.5, xmax = 91.5, ymin=-Inf, ymax=Inf,fill="#1F78B4", alpha = 0.15)+
	annotate("rect", xmin = 91.5, xmax = 100.5, ymin=-Inf, ymax=Inf,fill='#E31A1C', alpha = 0.15)+
	annotate("rect", xmin = 100.5, xmax = 117.5, ymin=-Inf, ymax=Inf,fill="#999999", alpha = 0.15)+
	annotate("rect", xmin = 117.5, xmax = 119.5, ymin=-Inf, ymax=Inf,fill='#AF5930', alpha = 0.15)+
	annotate("rect", xmin = 119.5, xmax = 141, ymin=-Inf, ymax=Inf,fill="#8DA0CB", alpha = 0.15)
p
#ggsave('~/Desktop/pPICRUST2 metabolism.pdf',p,width = 10.50,height = 6.5,dpi = 600,units = 'in')


#fig. 2B----
df <- read.delim('~/Desktop/molt github/Fig. 2/pathway.txt',row.names = 1)
group <- read.csv('~/Desktop/molt github/Fig. 2/group meta.csv',row.names = 1)
df2 <- cbind(t(df),group)
df3 <- melt(df2,id.vars = 'time',variable.name = 'name',value.name = 'value')
n <- unique(df3$name)
result <- NULL
for (i in 1: length(n)) {
	df4 <- subset(df3, name==n[i])
	if (sum(df4$value)==0) {
		i=i+1
		next
	}else{
		comparison<-with(df4,kruskal(value,time,group=TRUE))
		c <- comparison$means
		p <- comparison$groups
		length(unique(p$groups))>1
		if (length(unique(p$groups))==1) {
			i=i+1
			next
		}else{
			c$time <- rownames(c)
			c$funcs <-n[i]
			sig <- comparison$groups
			sig$time<-rownames(sig)
			sig <- sig[-1]
			names(sig) = c('sign','time')
			df5 <- merge(c,sig,by = 'time')
			result <- rbind(result,df5)}
	}
	print(i)
}
length(unique(result$funcs))
colnames(result)[11] <- "Pathway"
data_separated <- separate(result, Pathway, into = c("PathwayID", "Pathway"), sep = " ", extra = "merge")
#write.csv(data_separated ,'~/Desktop/molt github/Fig. 2/pathway significance.csv',row.names = F)

#plot sign functions by heatmap geom_tile----
m <- read.csv('~/Desktop/molt github/Fig. 2/KO1-5.csv')
m$Pathway <- gsub("\\[.*\\]","",m$Pathway)
m$Pathway <- trimws(m$Pathway, "right")
data <- read.csv('~/Desktop/molt github/Fig. 2/pathway significance.csv',header = TRUE, stringsAsFactors = FALSE)
f <- unique(data$Pathway)
result <- NULL
for (j in 1:length(f)) {
	res <- subset(data,Pathway==f[j])
	if (nrow(subset(m,Pathway==f[j]))==0) {
		res$Pathway2 <- 'Others'
	}else{
		res$Pathway2 <- unique(subset(m,Pathway==f[j])$PathwayL1)[1]
	}
	res$foldchange<- scale(res$value)
	result <- rbind(result,res)
}
result1 <- subset(result,!Pathway2=='Human Diseases')
result1 <- subset(result1,!Pathway2=='Others')
unique(result1$Pathway)
color_mapping <- c(c='#D41E26',bc='#AF5930',b="#FCAE6A",ab="#ADD9E8",a="#317BB3")
p3 <- ggplot(result1,aes(y=time,x=forcats::fct_reorder2(Pathway,Pathway,Pathway2)))+
	geom_tile(aes(fill=sign),width=2, height=1,alpha=1)+
	#scale_fill_gradient2(low = "#1F78B4",mid = 'white',high = "red",space = "Lab",
	#										 midpoint = 0,limits=c(-1.5,1.5),breaks=c(-1.5,-0.5,0,0.5,1.5))+
	scale_fill_manual('Significance',values = color_mapping)+
	scale_y_discrete(limits=c('48','12','0','CK'),expand = c(0,0.5))+
	labs(x=NULL,y=NULL)+#theme_bw()+
	#guides(fill=guide_colorbar(barwidth=5,barheight=0.5))+
	theme(panel.grid.major =  element_line(size=0.5),
				panel.background  = element_blank(),
				axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
				axis.title.x = element_text(angle=0,hjust=0.5,vjust = 0.5,size=8,color = 'black'),
				axis.title.y = element_text(angle=90,hjust=0.5,vjust = -1,size=8,color = 'black'),
				axis.text.y = element_text(angle=0,hjust=1,vjust = 1,size=8,color = 'black'),
				axis.text.x = element_blank(),
				legend.text = element_text(size=8,color = 'black'),
				legend.title= element_text(size=8,vjust = 0.5,hjust = 0),
				legend.key.size = unit(0.36,'cm'),legend.background = element_blank(),
				legend.position = 'none',
				plot.margin = margin(           # 调整边距
					t = 0.1,                       # 上边距
					r = 0.1,                       # 右边距
					b = 0.5,                       # 下边距
					l = 0.1,                       # 左边距
					unit = "in"                   # 单位（例如，“pt”，“cm”，“inch”等）
				)
	)+
	annotate("segment", y = 0.25, yend = 0.25, x = 243.5, xend = 273,colour = "black",linewidth=0.6)+
	#annotate("text", label='Cellular Processes',y = 0.25, size=2.4,x = 274,colour = "black",angle=45)+
	annotate("segment", y = 0.4, yend = 0.4, x = 216.5, xend = 243.5,colour = "black",linewidth=0.6)+
	annotate("segment", y = 0.3, yend = 0.3, x = 199.5, xend = 216.5,colour = "black",linewidth=0.6)+
	annotate("segment", y = 0.4, yend = 0.4, x = 62.5, xend = 199.5,colour = "black",linewidth=0.6)+
	annotate("segment", y = 0.3, yend = 0.3, x = 0, xend = 62.5,colour = "black",linewidth=0.6)
p3

#fig. 2C functional stability----
source('~/Desktop/molt github/molt plot/theme_linweichuan.R')
stability_onerep <- function(df, x){
	assertthat::assert_that(assertthat::has_name(df, x))
	assertthat::assert_that(is.numeric(df[[x]]))
	sync_var <- df[[x]]
	stability <- mean(sync_var, na.rm = TRUE)/stats::sd(sync_var, na.rm = TRUE)
	return(stability)
}
otu <- read.delim('~/Desktop/molt github/Fig. 2/enzyme.txt',row.names = 1,check.names = F)
group <- read.csv('~/Desktop/molt github/Fig. 2/group meta.csv')#[-12,]
otu1 <- otu[,colnames(otu) %in% group$SampleID] 
otu2 <- t(otu1[which(rowSums(otu1)>0),])
otu3 <-melt(otu2)
colnames(otu3) <- c("SampleID",'ASV','abundance')
otu3$SampleID <- factor(otu3$SampleID)#
otu3$abundance <- as.numeric(factor(otu3$abundance))# to avoid NA
####
result <- NULL
m <- unique(otu3$SampleID) 
for (j in 1:length(m)) {
	otu4 <- dplyr::filter(otu3,SampleID==m[j])
	otu4 <- otu4[which(otu4[['abundance']] > 0),]
	output <- stability_onerep(otu4, 'abundance')
	s <- data.frame(SampleID=m[j],stability=output)
	s1 <- merge(s,group,by = 'SampleID')
	result <- rbind(result,s1)
}
comparison <-with(result,kruskal(stability,time,group=TRUE))
cp <- comparison$means
cp$time <- rownames(cp)
sig <- comparison$groups
sig$time <-rownames(sig)
sig <- sig[-1]
df3 <- merge(cp[,c(1:4,10)],sig,by = 'time')
df3
df3$time <- factor(df3$time,levels = c('CK','0','3','6','9','12','24','48','72'))
result$time <- factor(result$time,levels = c('CK','0','3','6','9','12','24','48','72'))
p1 <- ggplot(result,aes(x=time,y=stability-0.05))+
	geom_boxplot(outlier.shape = NA,color='black',position=position_dodge(0.8),
							 linewidth=0.3,fatten = 1.5,fill='white')+
	geom_jitter(aes(color=time),size=1.5,alpha=0.8,stroke=0,width = 0.2)+
	scale_color_manual(values = c('#E64B35FF','#F39B7FFF','#91D1C2FF','#4DBBD5FF',
																'#00A087FF','#8491B4FF','#3C5488FF',"#984EA3",'#000000'))+
	geom_text(data=df3,aes(label=groups,y=stability-std-0.1),size=3.2)+
	labs(x = NULL,y= 'Functional stability')+
	theme_linweichuan()+
	theme(legend.position = 'none')
p1
#write.csv(result,'~/Desktop/molt github/Fig. 2/functional stability.csv',row.names = F)

#fig. 2D functional stability----
library(multifunc);library(agricolae)
otu <- read.delim('~/Desktop/molt github/Fig. 2/enzyme.txt',row.names = 1,check.names = F)
group <- read.csv('~/Desktop/molt github/Fig. 2/group meta.csv',row.names = 1)
a <- subset(read.csv('~/Desktop/molt github/Fig. 2/alpha stability.csv',row.names = 1),tissue=='Midgut')
otu1 <- otu[, colnames(otu)%in%rownames(group)]
otu1 <- otu1[which(rowSums(otu1)>0),]
otu1 <- otu1/colSums(otu1)*100
otu2<- merge(t(otu1),a,by='row.names') 
#head(otu2)
mf <-cbind(otu2[,colnames(otu2) %in% colnames(a)],getStdAndMeanFunctions(otu2, rownames(otu1)))
rownames(mf) <- rownames(group)
#write.csv(mf,'~/Desktop/molt github/Fig. 2/meta sign multifunctions.csv')
c<-with(mf,kruskal(meanFunction,time,group = T))
c
cp <- c$means
cp$time <- rownames(cp)
sig <- c$groups
sig$time <-rownames(sig)
sig <- sig[-1]
df3 <- merge(cp[,c(1:4,10)],sig,by = 'time')
df3
df3$time <- factor(df3$time,levels = c('CK','0','3','6','9','12','24','48','72'))
mf$time <- factor(mf$time,levels = c('CK','0','3','6','9','12','24','48','72'))
p2 <- ggplot(mf,aes(x=time,y=meanFunction))+
	geom_boxplot(outlier.shape = NA,color='black',position=position_dodge(0.8),
							 linewidth=0.3,fatten = 1.5,fill='white')+
	geom_jitter(aes(color=time),size=1.5,stroke=0,alpha=0.8,width = 0.2)+
	scale_color_manual(values = c('#E64B35FF','#F39B7FFF','#91D1C2FF','#4DBBD5FF',
																'#00A087FF','#8491B4FF','#3C5488FF',"#984EA3",'#000000'))+
	geom_text(data=df3,aes(label=groups,y=meanFunction-std-0.1),size=2.8)+
	labs(x = NULL,y= 'Multifunctionality')+
	theme_linweichuan()+theme(legend.position = 'none')
p2

#fig. 2E rfPermute----
library(rfPermute)
mf <- read.csv('~/Desktop/molt github/Fig. 2/meta sign multifunctions.csv',row.names = 1)
otu <- read.delim("~/Desktop/molt github/Fig. 2/pathway.txt",check.names = F,row.names = 1,header = T)
otu <- otu/colSums(otu)*100
group <- read.csv("~/Desktop/molt github/molt plot/groupMOLT mc.csv")
time <- c('CK','0','12','48')
otu1 <- otu[,colnames(otu)%in%group$SampleID] 
otu1 <- otu1[which(rowSums(otu1)>0),]
#match SampleID
matched_indices <- match(rownames(t(otu1)), rownames(mf))
matched_data2 <- mf[matched_indices,]
train1 <- data.frame(cbind(t(otu1),matched_data2[,c("meanFunction"),drop=F]))
colnames(train1)[ncol(train1)] <-"meanFunction"
train1$SampleID <- rownames(train1)
train2 <- merge(train1,group,by='SampleID')
set.seed(123)
otu_rfP <- rfPermute(meanFunction~., data =train2,ntree = 2000,importance = TRUE)#num.rep = 5000
print(otu_rfP)
importance_otu.scale <- data.frame(importance(otu_rfP, scale = TRUE), check.names = FALSE)
importance_otu.scale <- importance_otu.scale[order(importance_otu.scale$'%IncMSE', decreasing = TRUE), ]
importance_otu.scale$ASV_ID <- rownames(importance_otu.scale)
importance_otu.scale$ASV_ID <- factor(importance_otu.scale$ASV_ID, levels = importance_otu.scale$ASV_ID)
head(importance_otu.scale)
imp.select <- subset(importance_otu.scale, IncNodePurity.pval<0.05&`%IncMSE.pval`<0.05)#
#write.csv(imp.select,'~/Desktop/rfPermute4.csv',row.names = T)
nrow(subset(importance_otu.scale, IncNodePurity.pval<0.05&`%IncMSE.pval`<0.05))
nrow(subset(importance_otu.scale, `%IncMSE.pval`<0.05))

##plot function rfPermute----
source('~/Desktop/molt github/molt plot/theme_linweichuan.R')
library(reshape2)
df1 <- read.csv('~/Desktop/molt github/Fig. 2/rfPermute4.csv',check.names = F)
df1 <- df1[order(df1$`%IncMSE`,decreasing = T),]
df1$name <- 1
n <- unique(df1$Pathway)
df4 <- NULL
for (i in n) {
	res <- subset(df1,Pathway==i)
	wrapped_string <- strwrap(unique(res$Pathway), width = 40)
	res$Pathway <- paste(wrapped_string, collapse = "\n")
	df4 <- rbind(df4,res)
}
df4$Pathway
h1 <- ggplot(df4, aes(x=factor(name),y=reorder(Pathway,`%IncMSE`))) + 
	geom_tile(fill='transparent')+
	geom_text(aes(label=round(`%IncMSE`,2)),size=2.8)+
	labs(x=NULL,y=NULL,title = 'Importance')+
	theme_linweichuan()+
	theme(panel.grid = element_blank(),
				panel.border = element_blank(),
				legend.position = 'none',
				axis.text.x = element_blank(),   
				axis.text.y = element_blank(), 
				axis.ticks.x = element_blank(),	axis.ticks.y = element_blank(),
				plot.title = element_text(colour='black', size=8,angle = 0,hjust = 0.5),
				plot.margin = margin(           # 调整边距
					t = 0.1,                       # 上边距
					r = 0,                       # 右边距
					b = 0.1,                       # 下边距
					l = 0,                       # 左边距
					unit = "in"                   # 单位（例如，“pt”，“cm”，“inch”等）
				))
h1
#fig. 2F----
m <- read.csv('~/Desktop/molt github/Fig. 2/KO1-5.csv')
m$Pathway <- gsub("\\[.*\\]","",m$Pathway)
m$Pathway <- trimws(m$Pathway, "right")
df2 <- read.csv('~/Desktop/molt github/Fig. 2/pathway significance.csv',check.names = F)
n <- unique(df1$Pathway)
result <- NULL
for (i in n) {
	if (nrow(subset(df2,Pathway==i))==0) {
		next
	}else{
		res <- subset(df2,Pathway==i)
		res$Pathway2 <- unique(subset(m,Pathway==i)$PathwayL1)[1]
		res$foldchange<- res$value/min(subset(res,!value==0)$value)
		wrapped_string <- strwrap(unique(res$Pathway), width = 40)
		res$Pathway <- paste(wrapped_string, collapse = "\n")
		result <- rbind(result,res)
	}
}
result[which(result$Pathway2 %in% NA),'Pathway2'] <-  "Others"
color_mapping <- c(c='#D41E26',bc='#AF5930',b="#FCAE6A",ab="#ADD9E8",a="#317BB3")
#forcats::fct_reorder(funcs,pathway) #order funcs by pathway
result$Pathway <- factor(result$Pathway,levels = rev(df4$Pathway))
p <- ggplot(result,aes(x=time,y = Pathway))+
	#geom_tile(color='transparent',fill='transparent',linewidth=0.3)+
	scale_size(range = c(1.5,6),breaks = c(0.1,0.5,1,2))+
	geom_point(aes(size= log2(foldchange),color=sign),alpha=1,stroke=0)+#4 cross
	scale_color_manual(values = color_mapping)+
	scale_x_discrete(limits=rev(c('48','12','0','CK')))+
	labs(x=NULL,y=NULL,title = 'Abundance')+
	guides(size=guide_legend(title="Log2(Foldchange)",override.aes = list(nrow=2,ncol=3),order = 1),
				 color=guide_legend(title="Significance",reverse=F,override.aes = list(size=3,nrow=2,ncol=3),order = 2))+
	theme(panel.grid.major = element_line(linewidth=0.2,color='gray'),
				panel.border = element_blank(),
				panel.background  = element_blank(),text = element_text(family = 'Arial'),
				plot.title = element_text(angle=0,hjust=0.5,vjust = 0.5,size=8,color = 'black'),
				axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
				axis.title.x = element_text(angle=0,hjust=0.5,vjust = 0.5,size=8,color = 'black'),
				axis.title.y = element_text(angle=90,hjust=0.5,vjust = 0.5,size=8,color = 'black'),
				axis.text.x = element_text(angle=0,hjust=0.5,vjust = 0.5,size=8,color = 'black'),
				axis.text.y = element_text(angle=0,hjust=0.5,vjust = 0.5,size=8,color = 'black'),
				legend.text = element_text(size=8,color = 'black',hjust=0.5),legend.title= element_text(size=8),
				legend.key= element_blank(),legend.key.size = unit(0.2,'in'),
				legend.background = element_blank(),
				legend.position = 'none',#'bottom',
				legend.direction = "horizontal",
				plot.margin = margin(           # 调整边距
					t = 0.1,                       # 上边距
					r = 0.0,                       # 右边距
					b = 0.1,                       # 下边距
					l = 0.0,                       # 左边距
					unit = "in"                   # 单位（例如，“pt”，“cm”，“inch”等）
				))
p
#ggsave('~/Desktop/p.pdf',p,width = 9.12,height = 5.40,dpi = 600,units = 'in')

#fig. 2G----
##correction between community and function-------
source('~/Desktop/molt github/molt plot/theme_linweichuan.R')
library(psych);library(reshape2);library(dplyr);library(RColorBrewer)
df <- read.csv('~/Desktop/molt github/Fig. 2/rfPermute4.csv',row.names = 1,check.names = F)
func <- read.delim('~/Desktop/molt github/Fig. 2/pathway1.txt',row.names = 2,check.names = F)[-1]
func1 <- func[df$Pathway,]
group <- read.csv('~/Desktop/molt github/Fig. 2/group meta.csv',row.names = 1)
comm <- read.delim('~/Desktop/molt github/molt plot/ASV_34848.txt',row.names = 1,check.names = F)
#load community data
rf <- read.delim('~/Desktop/molt github/Fig. 1 and Fig. S11-S12/Mc importance_otu.txt',check.names = F,row.names = 1)
rf1 <- read.csv('~/Desktop/molt github/Fig. 1 and Fig. S11-S12/mc imp weight.csv',check.names = F,row.names = 1)
diff <- setdiff(rownames(rf),rownames(rf1))
otu <-comm[rownames(comm)%in%diff,colnames(comm)%in%rownames(group)] 
identical(colnames(func1),colnames(otu))
otu1 <- cbind(t(otu),t(func1))
otu2 <- scale(otu1[-12,])
occor = corr.test(otu2,use="pairwise",method="pearson",adjust="BH",alpha=0.05)
occor.r = occor$r[35:60,1:34] # 取相关性矩阵R值
occor.p = occor$p[35:60,1:34]  # 取相关性矩阵p值
occor.r[occor.p>=0.05]= 0
cor <- melt(occor.r)
cor1 <- cor[!cor$value==0,]
head(cor1)
#match order
df1 <- read.csv('~/Desktop/molt github/Fig. 2/rfPermute4.csv',row.names = 1,check.names = F)
df1 <- df1[order(df1$`%IncMSE`,decreasing = T),]
df1$name <- 1
n <- unique(df1$Pathway)
df4 <- NULL
for (i in n) {
	res <- subset(df1,Pathway==i)
	wrapped_string <- strwrap(unique(res$Pathway), width = 40)
	res$Pathway <- paste(wrapped_string, collapse = "\n")
	df4 <- rbind(df4,res)
}
#[df1order(df1$`%IncMSE`)
result <- NULL
for (i in 1:nrow(cor1)) {
	cor2 <- cor1[i,]
	wrapped_string <- strwrap(cor2$Var1, width = 40)
	cor2$Var1 <- paste(wrapped_string, collapse = "\n")
	# 输出结果
	result <- rbind(result,cor2)
}
result[which(result$value>0),'Correlation'] <- 'Positive'
result[which(result$value<0),'Correlation'] <- 'Negative'
result$Var1 <- factor(result$Var1,levels = rev(df4$Pathway))
result$Var2 <- factor(result$Var2,levels = c('ASV_444','ASV_127','ASV_1','ASV_177','ASV_445'))
result$Correlation <- factor(result$Correlation,levels = c('Positive','Negative'))
pcor <- ggplot(result,aes(Var2,Var1))+
	geom_tile(fill = 'transparent')+
	geom_point(aes(shape=Correlation),stroke=0.5,size=3.5)+#color=value,fill=value,size=abs(value),
	scale_shape_manual(values = c(2,6))+
	labs(x=NULL,y=NULL,title = "Pearson's correlation")+	scale_size(range = c(4,4),breaks = c(-1,0,1))+
	guides(color=guide_legend(title="Pearson's R",reverse = F,override.aes = list(order = 1)),
				 fill='none')+
	theme(
		panel.grid.major = element_blank(),
		panel.border = element_blank(),
		panel.background  = element_blank(),text = element_text(family = 'Arial'),
		plot.title = element_text(angle=0,hjust=0.5,vjust = 0.5,size=8,color = 'black'),
		axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
		axis.title.x = element_text(angle=0,hjust=0.5,vjust = 0.5,size=8,color = 'black'),
		axis.title.y = element_text(angle=90,hjust=0.5,vjust = 0.5,size=8,color = 'black'),
		axis.text.x = element_text(angle=45,hjust=0.8,vjust = 0.8,size=8,color = 'black'),
		axis.text.y = element_blank(),
		legend.text = element_text(size=8,color = 'black',hjust=0.5),legend.title= element_text(size=8),
		legend.key= element_blank(),legend.key.size = unit(0.2,'in'),
		legend.background = element_blank(),
		legend.box.background = element_rect(color = "gray"),
		legend.position = 'none',#c(0.8,0.5),
		plot.margin = margin(           # 调整边距
			t = 0.1,                       # 上边距
			r = 0.0,                       # 右边距
			b = 0.1,                       # 下边距
			l = 0.0,                       # 左边距
			unit = "in"                   # 单位（例如，“pt”，“cm”，“inch”等）
		))
pcor

library(aplot)
s1 <- h1 %>% 
	insert_right(p, width=4)%>%
	insert_right(pcor,width=3)
s1
#ggsave('~/Desktop/s1.pdf',s1,width = 5.75,height = 6.15,dpi = 600,units = 'in')


#fig. 2H and 2I -----
#devtools::install_github("cozygene/FEAST")
library(FEAST)
comm <- read.delim('~/Desktop/molt github/Fig. 2/ASV_table_norm feast.txt',row.names = 1,check.names = F)
comm1 <- comm[rownames(comm)%in%'ASV_444',]
comm2 <- data.frame(t(comm1))
comm2$ASV_5800<- 1;## one speices is not allowed, so generate a simulated species
comm2$ASV_5800<- as.integer(comm2$ASV_5800);
metadata <- Load_metadata('~/Desktop/molt github/Fig. 2/9hgroup_mc ASV_444.txt')
comm3 <- as.matrix(comm2[rownames(comm2) %in% rownames(metadata),])
comm3[rownames(subset(metadata,!SourceSink=='Sink')),][,2] <- 0
FEAST_output <- FEAST(C = comm3, metadata = metadata, different_sources_flag =0,COVERAGE = 22507,
											dir_path = '~/Desktop/',outfile=NULL,EM_iterations = 1000)
prop <- round(FEAST_output$data_prop,10);prop
unknown <- prop[nrow(prop),];unknown$source <- "Unknown"
prop1 <- cbind(prop[-nrow(prop),],subset(metadata,!SourceSink=='Sink')[,1])
prop2 = aggregate(prop1$pred_emnoise_all,by  =list(prop1$`subset(metadata, !SourceSink == "Sink")[, 1]`),FUN='sum')
prop3 = aggregate(prop1$pred_em_all,by  =list(prop1$`subset(metadata, !SourceSink == "Sink")[, 1]`),FUN='sum')
prop4 = cbind(prop2,prop3[,2])
colnames(prop4) <- c('source','pred_emnoise_all','pred_em_all')
prop5 = rbind(prop4,unknown)
#write.csv(prop5,'~/Desktop/molt github/Fig. 2/mc ASV444 source.csv',row.names = T)

##fig. 2H-----
library(dplyr);library(tidyverse);library(tidytext)
ds <- read.csv('~/Desktop/molt github/Fig. 2/mc ASV444 source.csv')
head(ds);source('~/Desktop/molt github/molt plot/theme_linweichuan.R')
ds$tissue <- 'Midgut'
ds1 <- ds %>%
	mutate(source = as.factor(source),
				 name = reorder_within(source, pred_emnoise_all,tissue),
				 name = fct_relevel(name, paste0("Unknown___", unique(source))))
pdmc <- ggplot(ds1,aes(y=name,x=pred_emnoise_all*100))+
	geom_col(width = 0.7,fill='#1F78B4')+
	scale_x_continuous(expand = c(0,2.5),limits = c(0,100),breaks=seq(0,100,25))+
	geom_text(data=filter(ds1,pred_emnoise_all>0.0),
						aes(label=round(pred_emnoise_all*100,2)),
						hjust=-0.1,size=2.8)+labs(y='Source',x='Contribution (%)')+
	scale_y_reordered() +
	theme_linweichuan()+#coord_trans(y="reverse")+
	theme(strip.background = element_blank(),strip.text = element_blank(),
				axis.text = element_text(size=8),axis.title.x = element_text(size=8),
				axis.title.y = element_text(size=8))
pdmc 

###plot ASV444 abundance----
library(reshape2)
comm <- read.delim('~/Desktop/molt github/Fig. 2/ASV_table_norm feast.txt',row.names = 1,check.names = F)
colSums(comm)
group <- subset(read.delim('~/Desktop/molt github/molt plot/groupALL.txt'),tissue=='Midgut')
comm1 <- comm[rownames(comm)%in%'ASV_444',colnames(comm) %in% group$SampleID]
comm2 <- melt(comm1,variable.name = 'SampleID',value.name = 'value')
df <- merge(comm2,group[,c('SampleID','time')],by='SampleID')
te1 = aggregate(df$value/225.07,                      
								by  =list(df$time),FUN='mean')
te2 = aggregate(df$value/225.07,                      
								by  =list(df$time),FUN='sd')
te <- cbind(te1,te2[2])
colnames(te) = c('time','mean','sd')
te$time <- factor(te$time,levels = rev(c('CK','0','3','6','9','12','24','48','72')))
p444 <- ggplot(te,aes(mean,time))+
	geom_col(width = 0.6,fill="#B2DF8A")+
	labs(y=NULL,x="Relative abundance (%)",title='ASV444 g_Pseudomonas')+
	scale_x_continuous(expand = c(0,0.005),limits = c(0,0.1),breaks = seq(0,0.1,0.05))+
	geom_text(data=subset(te,mean>0),
						aes(label=round(mean,3)),hjust=-0.1,size=2.2)+
	theme_linweichuan()+
	theme(axis.text = element_text(color='black',size=6,hjust = 0.5),
				axis.title.x = element_text(color='black',size=6,hjust = 0.5),
				plot.title = element_text(color='black',size=6,hjust = 0.5),
				legend.text = element_text(color='black',size=5,hjust = 0),
				legend.title = element_blank(),
				legend.key.size = unit(0.1,'cm'),	legend.position = 'none',
				plot.background = element_rect(fill = 'transparent',color = 'transparent'))
p444

library(grid)
pa1<-ggplotGrob(p444)
p444 <-pdmc+annotation_custom(pa1,xmin=30,xmax=100,ymin=0.5,ymax=7.5)
p444 
#
#ggsave('~/Desktop/pASV 444.pdf',p,width = 3.80,height = 2.75,dpi = 600,units = 'in')

##fig. 2I-----
#devtools::install_github("cozygene/FEAST")
library(FEAST)
comm <- read.delim('~/Desktop/molt github/Fig. 2/ASV_table_norm feast.txt',row.names = 1,check.names = F)
comm1 <- comm[rownames(comm)%in%'ASV_127',]
comm2 <- data.frame(t(comm1))
comm2$ASV_5800<- 1;## one speices is not allowed, so generate a simulated species
comm2$ASV_5800<- as.integer(comm2$ASV_5800);
metadata <- Load_metadata('~/Desktop/molt github/Fig. 2/48hgroup_mc.txt')
comm3 <- as.matrix(comm2[rownames(comm2) %in% rownames(metadata),])
comm3[rownames(subset(metadata,!SourceSink=='Sink')),][,2] <- 0
set.seed(123)
FEAST_output <- FEAST(C = comm3, metadata = metadata, different_sources_flag =0,COVERAGE = 22507,
											dir_path = '~/Desktop/',outfile=NULL,EM_iterations = 1000)
prop <- round(FEAST_output$data_prop,10);prop
unknown <- prop[nrow(prop),];unknown$source <- "Unknown"
prop1 <- cbind(prop[-nrow(prop),],subset(metadata,!SourceSink=='Sink')[,1])
prop2 = aggregate(prop1$pred_emnoise_all,by  =list(prop1$`subset(metadata, !SourceSink == "Sink")[, 1]`),FUN='sum')
prop3 = aggregate(prop1$pred_em_all,by  =list(prop1$`subset(metadata, !SourceSink == "Sink")[, 1]`),FUN='sum')
prop4 = cbind(prop2,prop3[,2])
colnames(prop4) <- c('source','pred_emnoise_all','pred_em_all')
prop5 = rbind(prop4,unknown)
#write.csv(prop5,'~/Desktop/molt github/Fig. 2/mc ASV127 source.csv',row.names = T)

source('~/Desktop/molt github/molt plot/theme_linweichuan.R')
library(dplyr);library(tidyverse);library(tidytext)
ds <- read.csv('~/Desktop/molt github/Fig. 2/mc ASV127 source.csv')
head(ds);
ds$tissue <- 'Midgut'
ds1 <- ds %>%
	mutate(source = as.factor(source),
				 name = reorder_within(source, pred_emnoise_all,tissue),
				 name = fct_relevel(name, paste0("Unknown___", unique(source))))
pdmc <- ggplot(ds1,aes(y=name,x=pred_emnoise_all*100))+
	geom_col(width = 0.7,fill='#1F78B4')+
	scale_x_continuous(expand = c(0,2.5),limits = c(0,100),breaks=seq(0,100,25))+
	geom_text(data=filter(ds1,pred_emnoise_all>0.0),
						aes(label=round(pred_emnoise_all*100,2)),
						hjust=-0.1,size=2.8)+labs(y='Source',x='Contribution (%)',
																			title = 'ASV127 s_Maritalea porphyrae')+
	scale_y_reordered() +
	theme_linweichuan()+
	theme(strip.background = element_blank(),strip.text = element_blank(),
				axis.text = element_text(size=8),axis.title.x = element_text(size=8),
				axis.title.y = element_text(size=8))
pdmc
###plot ASV127 abundance----
comm <- read.delim('~/Desktop/molt github/Fig. 2/ASV_table_norm feast.txt',row.names = 1,check.names = F)
colSums(comm)
group <- subset(read.delim('~/Desktop/molt github/molt plot/groupALL.txt'),tissue=='Midgut')
comm1 <- comm[rownames(comm)%in%'ASV_127',colnames(comm) %in% group$SampleID]
comm2 <- melt(comm1,variable.name = 'SampleID',value.name = 'value')
df <- merge(comm2,group[,c('SampleID','time')],by='SampleID')
te1 = aggregate(df$value/225.07,                      
								by  =list(df$time),FUN='mean')
te2 = aggregate(df$value/225.07,                      
								by  =list(df$time),FUN='sd')
te <- cbind(te1,te2[2])
colnames(te) = c('time','mean','sd')
te$time <- factor(te$time,levels = rev(c('CK','0','3','6','9','12','24','48','72')))
p127 <- ggplot(te,aes(mean,time))+
	geom_col(width = 0.6,fill="#FB9A99")+
	labs(y=NULL,x="Relative abundance (%)",title='Temporal change')+
	scale_x_continuous(expand = c(0,0.005),limits = c(0,0.1),breaks = seq(0,0.1,0.05))+
	geom_text(data=filter(te,mean>0),
						aes(label=round(mean,3)),hjust=-0.1,size=2.2)+
	theme_linweichuan()+
	theme(axis.text = element_text(color='black',size=6,hjust = 0.5),
				axis.title.x = element_text(color='black',size=6,hjust = 0.5),
				plot.title = element_text(color='black',size=6,hjust = 0.5),
				legend.text = element_text(color='black',size=5,hjust = 0),
				legend.title = element_blank(),
				legend.key.size = unit(0.1,'cm'),	legend.position = 'none', #c(0.7,0.6),
				plot.background = element_rect(fill = 'transparent',color = 'transparent'))
p127
library(grid)
pa1 <- ggplotGrob(p127)
p127 <-pdmc+annotation_custom(pa1,xmin=30,xmax=100,ymin=0.5,ymax=7.5)
p127
#ggsave('~/Desktop/pASV 127.pdf',p,width = 3.80,height = 2.75,dpi = 600,units = 'in')

#Finally, these plots are arranged into fig. 2 in AI software and further modify legends.


