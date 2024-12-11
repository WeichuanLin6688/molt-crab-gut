#fig. S9------
source('~/Desktop/molt github/molt plot/theme_linweichuan.R')
group <- read.delim('~/Desktop/molt github/molt plot/groupMOLT.txt')
site <- unique(group$tissue)
site
data  <- NULL
#generate a bray distance matrix in each gut group
for (i in site) {
	group1 <- subset(group,tissue==i)
	otu<- read.delim('~/Desktop/molt github/molt plot/ASV_34848.txt',row.names = 1,header = T,check.names = F)
	otu1 <- otu[,colnames(otu) %in% group1$SampleID]
	dis <-  as.matrix(vegdist(t(otu1),method = 'bray'))
	n <- unique(group1$time)
	for (j in n) {
		group2 <- subset(group1,time==j)
		dis1 <- dis[rownames(dis) %in% group2$SampleID,colnames(dis) %in% group2$SampleID]
		diag(dis1) <- NA
		dis2 <- melt(dis1)
		dis3 <-dis2[complete.cases(dis2),]
		dis3$tissue <- i;dis3$time <- j
		data <- rbind(data,dis3)
	}
}
#write.csv(data,'~/Desktop/molt github/Fig. S9/dis all.csv')

#plot best model line of dissimilarity----
library(ggtrendline);library(basicTrendline);library(segmented)
model <- c("line2P", "line3P", "log2P", "exp2P", "exp3P", "power2P", "power3P")
dis <- read.csv('~/Desktop/molt github/Fig. S9/dis all.csv',row.names = 1)
t <-  unique(dis$tissue)
colnames(dis)
result <- NULL
for (i in t) {
	dis1 <- filter(dis,tissue==i & !time=='CK')
	dis1[dis1$time==0,] <- 0.000001
	x <- as.numeric(dis1$time)
	y <- as.numeric(dis1$value)
	fit_out <- NULL
	for (j in 1:length(model)) {
		fit <- try(trendline_sum(x, y, model=model[j], summary=TRUE, eDigit=5))
		if("try-error"%in% class(fit)){
			j=j+1
			next
		}else{
			r <- data.frame(AIC=fit$AIC,model=model[j])
			fit_out <- rbind(fit_out,r)}
	}
	best=filter(fit_out,AIC==min(fit_out$AIC))$model
	fit <- try(trendline_sum(x, y, model=best, summary=TRUE, eDigit=5))
	if (fit$p.value<0.05) {
		r1 <-	data.frame(tissue= i,r2=fit$adj.R.squared,p=fit$p.value,
										 slope=fit$parameter$a,model=best)
		result<- rbind(result,r1)
	}else{
		i=i+1
		next
	}
}
result

##foregut community dissimilarity----
library(segmented);library(basicTrendline);library(ggtrendline)
dis <- read.csv('~/Desktop/molt github/Fig. S9/dis all.csv',row.names = 1)
df1<- filter(dis, !time=='CK')
df2<- filter(df1, tissue=='Foregut')
y1=as.numeric(df2$value)
x1=as.numeric(df2$time)
x1[x1==0] <- 0.000001
#z=as.numeric(df3$Shannon)
fit<- glm(y1~x1)#建立对数曲线方程
r <- segmented(fit)
summary(segmented(fit))#"#3c63ad","#cc6615"
f <- trendline_sum(sqrt(x1), y1, model="log2P", summary=TRUE, eDigit=5)
a=f$parameter$a;b=f$parameter$b
f$formula
q_x = r$psi[1];q_y = a*log(q_x,exp(1)) + b
br_x = r$psi[2];br_y = a*log(br_x,exp(1)) + b
sd = r$psi[3]
c1 <-ggtrendline(x1,y1,model = "log2P",linecolor = '#0E70BC',
								 linetype = 1,linewidth = 0.5,CI.lty = 0,CI.lwd = 0.2,
								 CI.level = 0.95,CI.fill ='transparent',CI.alpha = 1,CI.color = "black",
								 summary = TRUE,show.eq = TRUE,yhat = FALSE,eq.x = 54,eq.y = 0.7/3*4,
								 Rname = 0,Pname = 0,rrp.x = 54,rrp.y = 0.8,
								 text.col = "black",eDigit = 3,eSize = 2.8,xlab = NULL,ylab = NULL)+
	geom_point(aes(x1,y1),shape=1,size=2.5,stroke=0.3,color="gray25")+
	annotate(geom = "text", x = 5, y = 1, label='Foregut',size=2.8)+
	labs(x=NULL,y=NULL,title =NULL)+#
	scale_x_continuous(limits = c(0,72),breaks = c(0,3,6,9,12,24,48,72))+
	theme_linweichuan()
c1  
##midgut community dissimilarity----
dis <- read.csv('~/Desktop/molt github/Fig. S9/dis all.csv',row.names = 1)
df1<- filter(dis, !time=='CK')
df2<- filter(df1, tissue=='Midgut')
y2=as.numeric(df2$value)
x2=as.numeric(df2$time)
x2[x2==0] <- 0.000001
fit<- glm(y2~x2)#建立对数曲线方程
r <- segmented(fit)
summary(segmented(fit))
f <- trendline_sum(x2, y2, model="exp3P", summary=TRUE, eDigit=5)
a=f$parameter$a;b=f$parameter$b;c=f$parameter$c
q_x = r$psi[1];q_y = a*exp(b*q_x)+c
br_x = r$psi[2];br_y = a*exp(b*br_x)+c
sd = r$psi[3]

c2 <-ggtrendline(x2,y2,model = "exp3P",linecolor = '#0E70BC',
								 linetype = 1,linewidth = 0.5,CI.lty = 0,CI.lwd = 0.2,
								 CI.level = 0.95,CI.fill ='transparent',CI.alpha = 1,CI.color = "black",
								 summary = TRUE,show.eq = T,yhat = F,eq.x = 54,eq.y = 0.7/3*4,
								 Rname = 1,Pname = 0,rrp.x = 54,rrp.y = 0.8,
								 text.col = "black",eDigit = 3,eSize = 2.8,xlab = NULL,ylab = NULL)+
	annotate(geom = "rect",xmin= br_x-sd,xmax= br_x+sd,ymin=-Inf,ymax=Inf,fill='#D72B23',alpha=0.2)+
	annotate(geom = "point", x = q_x, y = q_y, color = '#D72B23',fill='white',size = 3,shape=4,stroke=0.8)+
	annotate(geom = "point", x = br_x, y = br_y, color = '#D72B23',fill='white',size = 3,shape=13,stroke=0.8)+
	geom_point(aes(x2,y2),shape=1,size=2.5,stroke=0.3,color="black")+
	#geom_point(data= d,aes(value,as.numeric(time)),shape=1,size=2.5,stroke=0.3,color="gray25")+
	labs(x=NULL,y=NULL,title =NULL)+#
	annotate(geom = "text", x = 5, y = 1, label='Midgut',size=2.8)+
	scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25))+
	scale_x_continuous(limits = c(0,72),breaks = c(0,3,6,9,12,24,48,72))+
	theme_linweichuan()
c2
##hindgut community dissimilarity----
dis <- read.csv('~/Desktop/molt github/Fig. S9/dis all.csv',row.names = 1)
df1<- filter(dis, !time=='CK')
df2<- filter(df1, tissue=='Hindgut')
y3=as.numeric(df2$value)
x3=as.numeric(df2$time)
x3[x3==0] <- 0.000001
fit<- glm(y3~x3)
r <- segmented(fit)
summary(segmented(fit))#"#3c63ad","#cc6615"
f <- trendline_sum(x3, y3, model="exp3P", summary=TRUE, eDigit=5)
a=f$parameter$a;b=f$parameter$b;c=f$parameter$c
q_x = r$psi[1];q_y = a*exp(b*q_x)+c
br_x = r$psi[2];br_y = a*exp(b*br_x)+c
sd = r$psi[3]
c3 <-ggtrendline(x3,y3,model = "exp3P",linecolor = '#0E70BC',
								 linetype = 1,linewidth = 0.5,CI.lty = 0,CI.lwd = 0.2,
								 CI.level = 0.95,CI.fill ='transparent',CI.alpha = 1,CI.color = "black",
								 summary = TRUE,show.eq = TRUE,yhat = FALSE,eq.x = 54,eq.y = 0.7/3*4,
								 Rname = 0,Pname = 0,rrp.x = 54,rrp.y = 0.8,
								 text.col = "black",eDigit = 3,eSize = 2.8,xlab = NULL,ylab = NULL)+
	annotate(geom = "rect",xmin= br_x-sd,xmax= br_x+sd,ymin=-Inf,ymax=Inf,fill='#D72B23',alpha=0.2)+
	annotate(geom = "point", x = q_x, y = q_y, color = '#D72B23',fill='white',size = 3,shape=4,stroke=0.8)+
	annotate(geom = "point", x = br_x, y = br_y, color = '#D72B23',fill='white',size = 3,shape=13,stroke=0.8)+
	geom_point(aes(x2,y2),shape=1,size=2.5,stroke=0.3,color="black")+
	annotate(geom = "text", x = 5, y = 1, label='Hindgut',size=2.8)+
	labs(x=NULL,y=NULL,title =NULL)+#
	scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25))+
	scale_x_continuous(limits = c(0,72),breaks = c(0,3,6,9,12,24,48,72))+
	theme_linweichuan()
c3
library(ggpubr)
c4 <- ggarrange(c1,c2,c3, align = 'hv',labels = c("A", "B","C"),nrow = 1,ncol=3,
								font.label = list(size = 10, color = "black", face = "plain"))
c4 <- annotate_figure(c4,left = text_grob("Community disimilarity", color = "black", rot = 90,size =8))
c4

#Finally, these plots are arranged into fig. S9 in AI software and further modify legends.


