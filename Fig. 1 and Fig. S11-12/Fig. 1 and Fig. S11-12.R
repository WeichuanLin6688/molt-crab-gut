#fig. 1 and fig. S11-12----
#RandomForest----
library(randomForest);library(ggsci)
otu <- read.delim("~/Desktop/molt github/molt plot/ASV_34848.txt",check.names = F,row.names = 1,header = T)
group<-read.csv("~/Desktop/molt github/molt plot/groupMOLT mc.csv",row.names = 1)
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
#write.table(imp, 'Hc importance_otu0325.txt', sep = '\t', col.names = NA, quote = FALSE)
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
#write.csv(train.cv.mean,'~/Desktop/molt github/Fig. 1 and Fig. S11-12/mc_train.cv.mean.csv',row.names = F)
#write.csv(train.cv.mean,'~/Desktop/molt github/Fig. 1 and Fig. S11-12/hc_train.cv.mean.csv',row.names = F)

#fig. 1A-D----
##fig. 1A----
library(splines)  
train.cv.mean <- read.csv('~/Desktop/molt github/Fig. 1 and Fig. S11-12/mc_train.cv.mean.csv')
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

##fig. 1B----
#find out optimal ASVs number
imp <- read.delim('~/Desktop/molt github/Fig. 1 and Fig. S11-12/Mc importance_otu.txt',row.names = 1)
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
#write.csv(n1,'mc rf NMDS.csv',row.names = T)
source('~/Desktop/molt github/molt plot/theme_linweichuan.R')
n1 <- read.csv('~/Desktop/molt github/Fig. 1 and Fig. S11-12/mc rf NMDS.csv',row.names = 1)
n1$time1 <- factor(n1$time1)
p4 <- ggplot(n1,aes(Dim.1,Dim.2))+
	geom_point(aes(color=time1),shape=19,size=2,stroke=0)+
	scale_color_manual(values =c("#BA3E2C","#0F72B2"))+
	labs(x='Dim 1',y='Dim 2')+
	theme_linweichuan()+
	theme(legend.title = element_blank(),legend.background = element_blank(),
				legend.position = 'none',legend.key = element_blank(),)
p4

##fig. 1C----
#exclude growth associated ASVs
library(randomForest);library(rfPermute)
otu <- read.delim("~/Desktop/molt github/molt plot/ASV_34848.txt",check.names = F,row.names = 1,header = T)
group<-read.csv("~/Desktop/molt github/molt plot/groupMOLT mc.csv",row.names = 1)
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
#write.csv(imp.select,'~/Desktop/molt github/Fig. 1 and Fig. S11-12/mc imp weight.csv',row.names = T)
#ggvenn
library(ggvenn)
df <- read.delim('~/Desktop/molt github/Fig. 1 and Fig. S11-12/mc importance_otu.txt',check.names = F,row.names = 1)
df1 <- read.csv('~/Desktop/molt github/Fig. 1 and Fig. S11-12/mc imp weight.csv',check.names = F,row.names = 1)
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
##fig. 3D----
#fistly get the last classified name
ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_qc mc hc.txt')
head(ann)
library(dplyr)
res <- NULL
for (i in 1:nrow(ann)) {
	ann1 <- ann[i,]
	n <- length(ann1[!ann1=="Unassigned"])
	ann1$last <- paste0(tolower(substr(colnames(ann1[n]),start = 1,stop = 1)),'_',ann1[n])
	res <- rbind(res,ann1)
}
#write.csv(res,'~/Desktop/molt github/Fig. 1 and Fig. S11-12/last classify.csv',row.names = F)

library(reshape2)
df <- read.delim('~/Desktop/molt github/Fig. 1 and Fig. S11-12/Mc importance_otu.txt',check.names = F,row.names = 1)
df1 <- read.csv('~/Desktop/molt github/Fig. 1 and Fig. S11-12/mc imp weight.csv',check.names = F,row.names = 1)
diff <- setdiff(rownames(df),rownames(df1))
ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_qc mc hc.txt')
otu <- read.delim("~/Desktop/molt github/molt plot/ASV_34848.txt",check.names = F,row.names = 1,header = T)
group <- read.csv("~/Desktop/molt github/molt plot/groupMOLT mc.csv",row.names = 1)
group <- subset(group,!time=='CK')
ann <- read.csv('~/Desktop/molt github/Fig. 1 and Fig. S11-12/last classify.csv')
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

##plot merge
library(ggpubr)
pf1 <- ggarrange(p1,p4,s1,ncol = 1,heights = c(0.3,0.3,0.2),labels = c('A','B','C'),
								 font.label = list(size = 10, face = "plain"))
pf2 <- ggarrange(pf1,pabu,ncol = 2,widths = c(0.309,0.5),labels = c('','D'),
								 font.label = list(size = 10, face = "plain"))
pf2
#ggsave('Fig. 1.pdf',pf2,width = 6.80,height = 6.25,dpi = 600,units ='in')

#fig. S11A-D----
##We get the Fig. S11 about the hindgut by using the above same codes. ----
##Please replace the file that starts with 'hc' with "mc".----
##fig. S11A----
library(splines) 
source('~/Desktop/molt github/molt plot/theme_linweichuan.R')
train.cv.mean <- read.csv('~/Desktop/molt github/Fig. 1 and Fig. S11-12/hc_train.cv.mean.csv')
p1 <- ggplot(train.cv.mean, aes(Group.1, x)) +
	geom_smooth(se = F,	method = 'lm', formula = y~ns(x,7),colour = "black",size=0.3) +##1F78B4
	geom_point(size=1.5,shape=19,stroke=0)+
	scale_x_continuous(expand = c(0,70),limits  = c(0,2800),breaks = seq(0,2800,700))+
	scale_y_continuous(expand = c(0,0.01),limits = c(0.1,0.40,0.1))+
	geom_vline(xintercept = 36,colour = '#000000', linetype = "dashed", size = 0.5)+
	geom_hline(yintercept = 0.154,colour = '#000000', linetype = "dashed", size = 0.5)+
	labs(title = NULL,x = 'Number of ASVs', y = 'Error rate of cross-validation')+
	annotate("text",x=700/1000*2800,y=0.15/0.16*0.3,label='Optimal = 20',color="black",size=2.8)+
	theme_bw(base_size = 8,base_family = 'Arial')+
	theme(panel.grid = element_blank(),panel.background = element_blank(),
				axis.ticks.x = element_line(linewidth = 0.3),
				axis.ticks.y = element_line(linewidth = 0.3),
				axis.text = element_text(color='black',size=8,hjust = 0.5))#+xlim(0,50)
p1
##fig. S11B----
source('~/Desktop/molt github/molt plot/theme_linweichuan.R')
n1 <- read.csv('~/Desktop/molt github/Fig. 1 and Fig. S11-12/hc rf NMDS.csv')
n1$time1 <- factor(n1$time1)
p4 <- ggplot(n1,aes(Dim.1,Dim.2))+
	geom_point(aes(color=time1),shape=19,size=2,stroke=0)+
	scale_color_manual(values =c("#BA3E2C","#0F72B2"))+#"#E41A1C" "#377EB8"
	labs(x='Dim 1',y='Dim 2')+
	theme_linweichuan()+
	theme(legend.title = element_blank(),legend.background = element_blank(),
				legend.position = 'none',legend.key = element_blank(),)
p4

##fig. S11C----
df <- read.delim('hc importance_otu0325.txt')
df1 <- read.csv('hc imp weight.csv',check.names = F,row.names = 1)
a1 <- list("Stage"=rownames(df),
					 "Crab weight"=df1$ASV_ID)
s1 <- ggvenn(a1,c("Stage","Crab weight"),show_percentage = F,
						 stroke_color = 'white',stroke_alpha = 0.5,stroke_size = 0.3,digits=1,
						 fill_alpha = 0.6, fill_color =c("#999999","#999999"),#c("#0072B5","#E18727FF","#20854EFF"),
						 set_name_color ="black",set_name_size = 2.8,text_size=2.8)
s1

##fig. S11D----
library(reshape2)
df <- read.delim('~/Desktop/molt github/Fig. 1 and Fig. S11-12/hc importance_otu.txt',check.names = F,row.names = 1)
df1 <- read.csv('~/Desktop/molt github/Fig. 1 and Fig. S11-12/hc imp weight.csv',check.names = F,row.names = 1)
diff <- setdiff(rownames(df),rownames(df1))
ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_qc mc hc.txt')
otu <- read.delim("~/Desktop/molt github/molt plot/ASV_34848.txt",check.names = F,row.names = 1,header = T)
group <- read.csv("~/Desktop/molt github/molt plot/groupMOLT hc.csv",row.names = 1)
group <- subset(group,!time=='CK')
ann <- read.csv('~/Desktop/molt github/Fig. 1 and Fig. S11-12/last classify.csv')
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
	labs(x=NULL,y=NULL)+#	scale_y_discrete()+
	guides(size=guide_legend(title="Log2(foldchange)",nrow=2,ncol=3,order = 2),
				 color=guide_legend(title = "Enriched in",size=4,order = 1))+
	theme_bw()+
	theme(strip.background = element_blank(),strip.text = element_blank(),
				panel.grid.major = element_line(linewidth=0.2,color='gray',alpha),
				#panel.border = element_rect(linewidth=0.3,color='black',fill='transparent'),
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
	annotate("rect", ymin = 0.3, ymax = 1.5, xmin=0, xmax=0.3,fill="#B2DF8A", alpha = 1)+
	#"Firmicutes"
	annotate("rect", ymin = 1.5, ymax = 6.5, xmin=0, xmax=0.3,fill="#FF7F00", alpha = 1)+
	#"Betaproteobacteria"
	annotate("rect", ymin = 6.5, ymax = 7.5, xmin=0, xmax=0.3,fill="#1F78B4", alpha = 1)+
	#"Bacteroidetes" 
	annotate("rect", ymin = 7.5, ymax = 10.5, xmin=0, xmax=0.3,fill="#6A3D9A", alpha = 1)+
	#"Alphaproteobacteria"
	annotate("rect", ymin = 10.5, ymax = 19.5, xmin=0, xmax=0.3,fill="#FB9A99", alpha = 1)
pabu

##
pf1 <- ggarrange(p1,p4,s1,ncol = 1,heights = c(0.3,0.3,0.2),labels = c('A','B','C'),
								 font.label = list(size = 10, face = "plain"))
pf2 <- ggarrange(pf1,pabu,ncol = 2,widths = c(0.309,0.5),labels = c('','D'),
								 font.label = list(size = 10, face = "plain"))
pf2
#ggsave('~/Desktop/Fig. S11.pdf',pf2,width = 6.80,height = 6.25,dpi = 600,units ='in')

#fig. 1E-M----
otu <- read.delim("~/Desktop/molt github/molt plot/ASV_34848.txt",check.names = F,row.names = 1,header = T)
group<-read.csv("~/Desktop/molt github/molt plot/groupMOLT mc.csv")
ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_qc mc hc.txt',row.names = 1)
#exclude weight
df <- read.delim('~/Desktop/molt github/Fig. 1 and Fig. S11-12/mc importance_otu.txt',check.names = F,row.names = 1)
df1 <- read.csv('~/Desktop/molt github/Fig. 1 and Fig. S11-12/mc imp weight.csv',check.names = F,row.names = 1)
n11 <- c('ASV_177','ASV_32','ASV_29','ASV_4743','ASV_14','ASV_1','ASV_444','ASV_233','ASV_132','ASV_7')
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
	# normalized = (x-min(x))/(max(x)-min(x)) 
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
				 	geom_jitter(aes(color=time),size=1.5,alpha=0.8,stroke=0,width = 0.2)+
				 	scale_color_manual(values = c('#E64B35FF','#F39B7FFF','#91D1C2FF','#4DBBD5FF',
				 																'#00A087FF','#8491B4FF','#3C5488FF',"#984EA3",'#000000'))+
				 	geom_text(data=df3,aes(label=groups,y=stability+std+0.1),size=2.8)+
				 	labs(x = NULL,y= 'Community stability')+
				 	theme_linweichuan()+
				 	theme(legend.position = 'none')
	)
	cat("正在进行",i,'\n')
}
ggarrange(p1,p2,p3,nrow = 1)
##mc---
otu <- read.delim("~/Desktop/molt github/molt plot/ASV_34848.txt",check.names = F,row.names = 1,header = T)
group<-read.csv("~/Desktop/molt github/molt plot/groupMOLT mc.csv")
ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_qc mc hc.txt',row.names = 1)
#exclude weight
df <- read.delim('~/Desktop/molt github/Fig. 1 and Fig. S11-12/mc importance_otu.txt',check.names = F,row.names = 1)
df1 <- read.csv('~/Desktop/molt github/Fig. 1 and Fig. S11-12/mc imp weight.csv',check.names = F,row.names = 1)
n11 <- c('ASV_177','ASV_32','ASV_29','ASV_4743','ASV_14','ASV_1','ASV_444','ASV_233','ASV_132','ASV_7')#early ASV
#divide three category
diff <- setdiff(rownames(df),rownames(df1))
diff <- diff[-12]
###
n <- 1:3;n
for (i in n) {
	if (i==1) {
		n1 <- n11
	}
	if (i==2) {
		n1 <- rownames(otu[!rownames(otu) %in% diff,]) #excact late ASV
	}
	if (i==3) {
		n1 <- diff[!diff %in% n11]
	}
	otu1 <- otu[rownames(otu) %in% n1,colnames(otu) %in% group$SampleID] 
	richness <- estimateR(t(otu1))[1, ]
	dat <- data.frame(cbind(richness,shannon))
	dat$SampleID <- rownames(dat)
	group2 <- group[,-c(3:5)]
	dat1 <- merge(dat,group2,by='SampleID')
	colnames(dat1)[2] <- c('ASV richness','Shannon index')
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
				 ggplot(te,aes(group,mean,fill=Taxonomy))+#,alluvium = Taxonomy, stratum = Taxonomy
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
#range(subset(te,!group=='CK')$mean)
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
p <- ggarrange(p1,p3,p2,b1,b6,b2,
							 `dASV richness1`,`dASV richness3`,`dASV richness2`,
							 nrow = 3,ncol=3,legend = 'none',align = 'hv',
							 labels = c('A','D','G','B','E','H','C','F','I','D','',''),
							 font.label = list(size = 10, face = "plain"))
p
#ggsave('~/Desktop/pFig 1E-M.pdf',p,width = 7.80,height = 7.45,dpi = 600,units = 'in')




#fig. S12A-D----
otu <- read.delim("~/Desktop/molt github/molt plot/ASV_34848.txt",check.names = F,row.names = 1,header = T)
group<-read.csv("~/Desktop/molt github/molt plot/groupMOLT hc.csv")
ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_qc mc hc.txt',row.names = 1)
#exclude weight
df <- read.delim('~/Desktop/molt github/Fig. 1 and Fig. S11-12/hc importance_otu.txt',check.names = F,row.names = 1)
df1 <- read.csv('~/Desktop/molt github/Fig. 1 and Fig. S11-12/hc imp weight.csv',check.names = F,row.names = 1)
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
	# normalized = (x-min(x))/(max(x)-min(x)) 
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
				 	geom_jitter(aes(color=time),size=1.5,alpha=0.8,stroke=0,width = 0.2)+
				 	scale_color_manual(values = c('#E64B35FF','#F39B7FFF','#91D1C2FF','#4DBBD5FF',
				 																'#00A087FF','#8491B4FF','#3C5488FF',"#984EA3",'#000000'))+
				 	geom_text(data=df3,aes(label=groups,y=stability+std+0.1),size=2.8)+
				 	labs(x = NULL,y= 'Community stability')+
				 	theme_linweichuan()+
				 	theme(legend.position = 'none')
	)
	cat("正在进行",i,'\n')
}
ggarrange(p1,p2,p3,nrow = 1)

###hc---
otu <- read.delim("~/Desktop/molt github/molt plot/ASV_34848.txt",check.names = F,row.names = 1,header = T)
group<-read.csv("~/Desktop/molt github/molt plot/groupMOLT hc.csv")
ann <- read.delim('~/Desktop/molt github/molt plot/taxonomy_qc mc hc.txt',row.names = 1)
#exclude weight
df <- read.delim('~/Desktop/molt github/Fig. 1 and Fig. S11-12/hc importance_otu.txt',check.names = F,row.names = 1)
df1 <- read.csv('~/Desktop/molt github/Fig. 1 and Fig. S11-12/hc imp weight.csv',check.names = F,row.names = 1)
n11 <- c('ASV_31','ASV_69','ASV_1','ASV_444','ASV_62')#early ASV
#n11 <- c('ASV_177','ASV_32','ASV_29','ASV_4743','ASV_14','ASV_1','ASV_444','ASV_233','ASV_132','ASV_7')
#divide three category
diff <- setdiff(rownames(df),rownames(df1))
diff <- diff[-12]
###
n <- 1:3;n
for (i in n) {
	if (i==1) {
		n1 <- n11
	}
	if (i==2) {
		n1 <- rownames(otu[!rownames(otu) %in% diff,]) #excact late ASV
	}
	if (i==3) {
		n1 <- diff[!diff %in% n11]
	}
	otu1 <- otu[rownames(otu) %in% n1,colnames(otu) %in% group$SampleID] 
	richness <- estimateR(t(otu1))[1, ]
	dat <- data.frame(cbind(richness,shannon))
	dat$SampleID <- rownames(dat)
	group2 <- group[,-c(3:5)]
	dat1 <- merge(dat,group2,by='SampleID')
	colnames(dat1)[2] <- c('ASV richness','Shannon index')
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
				 ggplot(te,aes(group,mean,fill=Taxonomy))+#,alluvium = Taxonomy, stratum = Taxonomy
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
#range(subset(te,!group=='CK')$mean)
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
p <- ggarrange(p1,p2,p3,b1,b2,b6,
							 `dASV richness1`,`dASV richness2`,`dASV richness3`,
							 nrow = 3,ncol=3,legend = 'none',align = 'hv',
							 labels = c('A','D','G','B','E','H','C','F','I','D','',''),
							 font.label = list(size = 10, face = "plain"))
p
#ggsave('~/Desktop/pFigS12.pdf',p,width = 7.80,height = 7.45,dpi = 600,units = 'in')



