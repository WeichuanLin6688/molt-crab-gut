#fig. S13----
#functions codyn----
library(codyn)
group <- read.delim('~/Desktop/molt github/molt plot/groupCOYDN.txt',row.names = 1)
#读入物种数据
site <- unique(group$tissue)
site
res <- NULL
rc.res <- NULL
for (s in site) {
	group1 <- subset(group,tissue==s)
  df <-  read.delim('~/Desktop/molt github/Fig. S13/KEGG.Pathway.raw.txt',check.names = F,row.names = 1)
  df1 <- df[,colnames(df) %in% rownames(group1)]
	df21 <- df1
	df21 <-df21[complete.cases(df21),]
	otu <- df21
	otu1 <- t(otu[which(rowSums(otu)>0),])
	otu2 <- merge(otu1,group[,-2,drop=F],by='row.names')[-1]
	colnames(otu2)[ncol(otu2)-1] <- 'subplot'
	colnames(otu2)[ncol(otu2)] <- 'year'
	otu3<- melt(otu2,id.vars = c("subplot","year"),variable.name =  "species",value.name = "abundance")
	otu3$year <- as.integer(otu3$year)
	otu3$abundance <- as.numeric(otu3$abundance)
	otu3$subplot <- factor(otu3$subplot)
	otu3 <- otu3[order(-otu3$year,decreasing = T),]
	t <- unique(otu3$subplot)
	n <- unique(otu3$year)
	for (i  in n) {
		otu4 <- subset(otu3,year==i)
		for (j in t) {
			otu5 <- subset(otu4,subplot==j)
			for (k in t) {
				otu6 <- subset(otu3,subplot==k&!year==i)
				otu7 <- rbind(otu5,otu6)
				otu7$year <- as.integer(otu7$year)
				otu7$abundance <- as.numeric(otu7$abundance)
				otu7 <- otu7[order(-otu7$year,decreasing = T),]
				total <- turnover(otu7,time.var = "year", species.var = "species",abundance.var = "abundance",
													replicate.var = NA,metric="total")
				dis <- turnover(otu7,time.var = "year", species.var = "species",abundance.var = "abundance",
												replicate.var =NA,metric="disappearance")
				app <- turnover(otu7,time.var = "year", species.var = "species",abundance.var = "abundance",
												replicate.var = NA,metric="appearance")
				rs <- rank_shift(otu7,time.var = "year", species.var = "species",abundance.var = "abundance",
												 replicate.var = NA)
				rc <- rate_change_interval(otu7, time.var = "year", species.var = "species",abundance.var = "abundance",
																	 replicate.var = NA) 
				rc$subplot <- k
				rc$tissue <- s
				rc.res <- rbind(rc.res,rc)
				#var <- variance_ratio(otu7,time.var = "year",species.var = "species",abundance.var = "abundance",
				#											bootnumber = 1,replicate = NA,average.replicates = F)
				res1 <- cbind(total,dis,app,rs)
				res1$tissue <- s
				res1$group1 <- i
				res1$group2 <- j
				res1$group3 <- k
				res <- rbind(res,res1)
				cat("正在处理任务", s, i, "，当前为：", j, "，重复为：", k,"\n")
			}
		}}
}
#
#write.csv(res,'~/Desktop/molt github/Fig. S13/function all turnover.csv',row.names = F)
#write.csv(rc.res,'~/Desktop/molt github/Fig. S13/function rate_change_interval.csv',row.names = F)
unique(res$tissue)

##fig. S13A-----
res <- read.csv('~/Desktop/molt github/Fig. S13/function all turnover.csv')
te4 = aggregate(res$MRS,                      
								by  =list(res$year,res$tissue),FUN='mean')
te4$category <- "rank"
names(te4) <- c('time',"tissue",'Rank','category')
te4$tissue <- factor(te4$tissue,levels = site)
p1 <- ggplot(te4,aes(time,Rank/3))+ # to clearly show the fold change, we standardize value by Rank/3.
	geom_line(size=0.5,color='#3C76AF')+
	labs(x=NULL,y='Mean rank shift')+
	facet_wrap(~tissue,scale='free',nrow = 5,ncol = 1)+
	scale_y_continuous(limits = c(0,8.1),breaks = seq(0,8.1,2))+
	scale_x_continuous(breaks = c(0,6,12,24,48,72))+
	theme_linweichuan()+
	theme(legend.position = 'none',strip.background = element_blank(),
				strip.text = element_text(size=8))	
p1 
##fig. S13B-----
rc.res <- read.csv('~/Desktop/molt github/Fig. S13/function rate_change_interval.csv')
rc.res$tissue <- factor(rc.res$tissue,levels = site)
rc.res <- subset(rc.res,!tissue=="<NA>")
p2 <- ggplot(rc.res ,aes(interval,distance/1000))+
	geom_jitter(shape=1,stroke=0.3,size=1.5,height = 0.5,width = 0)+
	facet_wrap(~tissue,scales = 'free',ncol = 1)+
	geom_smooth(formula = y~x,se = F,method = 'lm',size=0.5,color= '#3C76AF')+
	ggpmisc::stat_poly_eq(
		aes(label = paste(after_stat(adj.rr.label))),#..eq.label..,
		formula = y ~ x, parse = TRUE,  small.r = F,small.p = T,label.x.npc = 6.95,label.y.npc = 0.9,
		rr.digits = 3,p.digits = 3,size = 2.8) +
	scale_x_continuous(breaks = seq(0,8,1))+
	scale_y_continuous(limits = c(0,8000),breaks = seq(0,8000,2000))+
	labs(x=NULL,y='Euclidean distance')+#'Time interval'
	theme_linweichuan()+
	theme(legend.position = 'none',strip.background = element_blank(),
				strip.text = element_text(size=8))	
p2
library(patchwork)
p3 <-p1+p2
#ggsave('~/Desktop/pPfuns.pdf',p3,width = 4.45 ,height = 4.90 ,dpi = 600,units = 'in')

#Finally, these plots are arranged into fig. S13 in AI software and further modify legends.