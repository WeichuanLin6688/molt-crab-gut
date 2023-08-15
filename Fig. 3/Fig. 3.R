#Fig. 3----
setwd("~/Desktop/Fig. 3")
#RandomForest----
library(randomForest);library(ggsci)
otu <- read.delim("~/Desktop/molt plot/ASV_34848.txt",check.names = F,row.names = 1,header = T)
group<-read.csv("~/Desktop/molt plot/groupMOLT mc.csv",row.names = 1)
group<- subset(group,!time=='CK')
group$weight <- group$wt-group$wo
otu1 <-otu[,colnames(otu)%in%rownames(group)] 
otu1 <-(otu1[which(rowSums(otu1)>0),])
train <- cbind(t(otu1),group[,c('time1')])
colnames(train)[ncol(train)] <- 'time'
train <- data.frame(train)
train$time <- factor(train$time)
set.seed(123)
otu_forest <- randomForest(time~.,data =train,ntree = 2000,importance = TRUE,proximity=TRUE)# nPerm = 5000,
print(otu_forest);plot(otu_forest)
imp <- data.frame(round(importance(otu_forest, scale = TRUE),2), check.names = FALSE)
imp <- imp[order(imp$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(imp)
#write.table(imp, 'Mc importance_otu0325.txt', sep = '\t', col.names = NA, quote = FALSE)
#cross-validation
set.seed(123)
train.cv <- replicate(5, rfcv(train[-ncol(train)], train$time, cv.fold = 10,step = 1.5), simplify = FALSE)
train.cv <- data.frame(sapply(train.cv, '[[', 'error.cv'))
train.cv$otus <- rownames(train.cv)
train.cv <- reshape2::melt(train.cv, id = 'otus')
train.cv$otus <- as.numeric(as.character(train.cv$otus))
train.cv.mean <- aggregate(train.cv$value, by = list(train.cv$otus), FUN = mean)
head(train.cv.mean, 20)#check the number of OTUs corresponding to the minimum X value
plot(train.cv.mean)
min(train.cv.mean$x)
#write.csv(train.cv.mean,'mc_train.cv.mean1',row.names = F)

#fig. 3a----
train.cv.mean <- read.csv('mc_train.cv.mean1.csv')
library(splines)  
p1 <- ggplot(train.cv.mean, aes(Group.1, x)) +
	geom_smooth(se = F,	method = 'lm', formula = y~ns(x,7),colour = "black",size=0.3) +##1F78B4
	geom_point(size=1.5,shape=19,stroke=0)+
	scale_x_continuous(limits  = c(0,1000))+
	scale_y_continuous(expand = c(0,0.005),limits = c(0.08,0.24,0.04))+
	geom_vline(xintercept = 36,colour = '#000000', linetype = "dashed", size = 0.5)+
	geom_hline(yintercept = 36,colour = '#000000', linetype = "dashed", size = 0.5)+
	labs(title = NULL,x = 'Number of ASVs', y = 'Error rate of cross-validation')+
	annotate("text",x=700,y=0.15,label='Optimal = 36',color="black",family='Arial',size=2.8)+
	theme_bw(base_size = 8)+
	theme(panel.grid = element_blank(),panel.background = element_blank(),
				axis.ticks.x = element_line(linewidth = 0.3),
				axis.ticks.y = element_line(linewidth = 0.3),
				axis.text = element_text(color='black',size=8,hjust = 0.5))
p1

#fig. 3b----
#find out optimal ASVs number
imp <- read.delim('Mc importance_otu0325.txt',row.names = 1)
train1 <- train[,colnames(train) %in% c(rownames(imp)[1:36],'time')]
set.seed(123)
otu_select <- randomForest(time~.,data = train1 ,ntree = 2000,importance = TRUE,proximity=TRUE)# nPerm = 5000,
print(otu_select);plot(otu_select)
par(mar = c(4, 1,1, 1)) 
n <- MDSplot(otu_select,train$time,k = 2, palette = pal_nejm()(2))
n1 <- cbind(n$point,group)
imp_select <- data.frame(round(importance(otu_select, scale = TRUE),2), check.names = FALSE)
imp_select <- imp_select[order(imp_select$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(imp_select)
#write.csv(n1,'mc rf0606 NMDS.csv',row.names = T)
source('~/Desktop/molt plot/theme_linweichuan.R')
n1 <- read.csv('mc rf0606 NMDS.csv',row.names = 1)
n1$time1 <- factor(n1$time1)
p4 <- ggplot(n1,aes(Dim.1,Dim.2))+
	geom_point(aes(color=time1),shape=19,size=2,stroke=0)+
	scale_color_manual(values =c("#BA3E2C","#0F72B2"))+
	labs(x='Dim 1',y='Dim 2')+
	theme_linweichuan()+
	theme(legend.title = element_blank(),legend.background = element_blank(),
				legend.position = 'none',legend.key = element_blank(),)
p4

#fig. 3c----
#exclude growth associated ASVs
library(randomForest);library(rfPermute)
otu <- read.delim("~/Desktop/molt plot/ASV_34848.txt",check.names = F,row.names = 1,header = T)
group<-read.csv("~/Desktop/molt plot/groupMOLT mc.csv",row.names = 1)
group<- filter(group,!time=='CK')
group$weight <- group$wt-group$wo
otu1 <-otu[,colnames(otu)%in%rownames(group)] 
otu1 <-otu1[which(rowSums(otu1)>0),]
train1 <- data.frame(cbind(t(otu1),group[,c('weight')]))
colnames(train1)[ncol(train1)] <-"weight"
#randomfoest
set.seed(123)
otu_forest <- randomForest(weight~., data =train1,ntree = 2000,importance = TRUE,proximity=TRUE)
print(otu_forest)
#rfPermute
set.seed(123)
otu_rfP <- rfPermute(weight~., data =train1,ntree = 500,importance = TRUE)
print(otu_rfP)
importance_otu.scale <- data.frame(importance(otu_rfP, scale = TRUE), check.names = FALSE)
importance_otu.scale <- importance_otu.scale[order(importance_otu.scale$'%IncMSE', decreasing = TRUE), ]
importance_otu.scale$ASV_ID <- rownames(importance_otu.scale)
importance_otu.scale$ASV_ID <- factor(importance_otu.scale$ASV_ID, levels = importance_otu.scale$ASV_ID)
head(importance_otu.scale)
imp.select <- filter(importance_otu.scale, IncNodePurity.pval<0.05&`%IncMSE.pval`<0.05)
#write.csv(imp.select,'mc imp weight.csv',row.names = T)
#ggvenn
library(ggvenn)
df <- read.delim('mc importance_otu0325.txt',check.names = F,row.names = 1)
df1 <- read.csv('mc imp weight.csv',check.names = F,row.names = 1)
n11 <- c('ASV_177','ASV_32','ASV_29','ASV_4743','ASV_14','ASV_1','ASV_444','ASV_233','ASV_132','ASV_7')
diff <- setdiff(rownames(df),rownames(df1))
n1 <- diff[!diff %in% n11]
a1 <- list("Stage"=rownames(df),
					 "Crab weight"=df1$ASV_ID)
s1 <- ggvenn(a1,c("Stage","Crab weight"),show_percentage = F,
						 stroke_color = 'white',stroke_alpha = 0.5,stroke_size = 0.3,digits=1,
						 fill_alpha = 0.6, fill_color =c("#999999","#999999"),
						 set_name_color ="black",set_name_size = 2.8,text_size=2.8)
s1

#fig. 3d plot ASV abundance----
#fistly get the last classified name
ann <- read.delim('~/Desktop/molt plot/taxonomy_qc mc hc.txt')
head(ann)
library(dplyr)
res <- NULL
for (i in 1:nrow(ann)) {
	ann1 <- ann[i,]
	n <- length(ann1[!ann1=="Unassigned"])
	ann1$last <- paste0(tolower(substr(colnames(ann1[n]),start = 1,stop = 1)),'_',ann1[n])
	res <- rbind(res,ann1)
}
#write.csv(res,'last classify.csv',row.names = F)

library(reshape2)
df <- read.delim('Mc importance_otu0325.txt',check.names = F,row.names = 1)
df1 <- read.csv('mc imp weight.csv',check.names = F,row.names = 1)
diff <- setdiff(rownames(df),rownames(df1))
ann <- read.delim('~/Desktop/molt plot/taxonomy_qc mc hc.txt')
otu <- read.delim("~/Desktop/molt plot/ASV_34848.txt",check.names = F,row.names = 1,header = T)
group <- read.csv("~/Desktop/olt plot/groupMOLT mc.csv",row.names = 1)
group <- subset(group,!time=='CK')
ann <- read.csv('last classify.csv')
otu1 <-otu[rownames(otu)%in%rownames(df),colnames(otu)%in%rownames(group)] 
colSums(otu1)
otu2 <- cbind(t(otu1),group[1])
colnames(otu2[ncol(otu2)]) <- 'time'
otu3 <- melt(otu2,id.vars = 'time',variable.name = 'ASV_ID',value.name = 'value')
te1 = aggregate(otu3$value/348.48,                      
								by  =list(otu3$ASV_ID,otu3$time),FUN='mean')
te2 = aggregate(otu3$value/348.48,                      
								by  =list(otu3$ASV_ID,otu3$time),FUN='sd')
te <- cbind(te1,te2[3])
colnames(te) <- c('ASV_ID','time','mean','sd')
otu4 <- merge(te,ann[,colnames(ann) %in% c('ASV_ID','PC','last')],by = 'ASV_ID')
otu5 <- tidyr::unite(otu4,"last", ASV_ID, last,sep = " ",remove = F)
otu5 <- subset(otu5,!PC=='Unassigned')
f <- unique(otu5$last)
result <- NULL
for (j in 1:length(f)) {
	res3 <- subset(otu5,last==f[j])
	res3$foldchange<- res3$mean/min(subset(res3,!mean==0)$mean)
	res3$time <- as.numeric(res3$time)
	if (	mean(subset(res3,time < 47)$mean)<	mean(subset(res3,time > 47)$mean)) {
		res3$trend <- 1;
	}else{res3$trend <- 2}
	result <- rbind(result,res3)
}
result$time <-factor(result$time,levels = sort(unique(result$time)),
										 labels = c('0','3','6','9','12','24','48','72')) 
log2(range(result$foldchange))
head(result)
result <- result[order(-result[[8]],result[[6]],result[[2]]),]
result$last <- factor(result$last,levels =unique(result$last))
result$trend <- factor(result$trend,levels = c(2,1),labels = c("Early stage (0-24 h)","Late stage (48-72 h)"))

pabu <- ggplot(result,aes(time,forcats::fct_reorder2(last,trend,PC)))+
	geom_tile(color='transparent',fill='transparent')+
	geom_point(shape=19,aes(size=log2(foldchange),color=factor(trend)),alpha=0.95,stroke=0)+
	scale_size(range = c(0.01,6),breaks = c(0.5,1,2,4,6,8,10))+
	scale_color_manual(values =c("#BA3E2C","#0F72B2"))+
	scale_y_discrete(position = 'right')+
	labs(x=NULL,y=NULL)+
	guides(size=guide_legend(title="Log2(foldchange)",nrow=2,ncol=3,order = 2),#"Log2(foldchange)"
				 color=guide_legend(title = "Enriched in",size=4,order = 1))+theme_bw()+
	theme(strip.background = element_blank(),strip.text = element_blank(),
				panel.grid.major = element_line(linewidth=0.2,color='gray',alpha),
				panel.border = element_blank(),
				panel.background  = element_blank(),
				axis.ticks.x = element_blank(),axis.ticks.y = element_blank(),
				axis.title.x = element_text(angle=0,hjust=0.5,vjust = 0.5,size=8,color = 'black'),
				axis.title.y = element_text(angle=90,hjust=0.5,vjust = 0.5,size=8,color = 'black'),
				axis.text.x = element_text(angle=0,hjust=0.5,vjust = 1,size=8,color = 'black'),
				axis.text.y = element_text(angle=0,hjust=1,vjust = 0.5,size=7,color = 'black'),
				legend.text = element_text(size=8,color = 'black',hjust=0),legend.title= element_text(size=8,hjust=0.5),
				legend.key  = element_blank(),legend.background = element_blank(),
				legend.position = 'bottom',legend.direction = 'vertical')+
	#add PC annotation
	#"Gammaproteobacteria"
	annotate("rect", ymin = 0.3, ymax = 12.5, xmin=0, xmax=0.3,fill="#B2DF8A", alpha = 1)+#"#8D4093"
	#"Firmicutes"
	annotate("rect", ymin = 12.5, ymax = 14.5, xmin=0, xmax=0.3,fill="#FF7F00", alpha = 1)+#"#FFFF99"
	#"Betaproteobacteria"
	annotate("rect", ymin = 14.5, ymax = 23.5, xmin=0, xmax=0.3,fill="#1F78B4", alpha = 1)+#'#3C5488'
	#"Bacteroidetes" 
	annotate("rect", ymin = 23.5, ymax = 28.5, xmin=0, xmax=0.3,fill="#6A3D9A", alpha = 1)+#"#59BD92"
	#"Alphaproteobacteria"
	annotate("rect", ymin = 28.5, ymax = 34.5, xmin=0, xmax=0.3,fill="#FB9A99", alpha = 1)#"#FF7F00"
pabu

##
library(ggpubr)
pf1 <- ggarrange(p1,p4,s1,ncol = 1,heights = c(0.3,0.3,0.2),labels = c('A','B','C'),
								 font.label = list(size = 10, face = "plain"))
pf2 <- ggarrange(pf1,pabu,ncol = 2,widths = c(0.309,0.5),labels = c('','D'),
								 font.label = list(size = 10, face = "plain"))
pf2
#ggsave('pmc RFCM3.pdf',pf2,width = 6.80,height = 6.25,dpi = 600,units ='in')
