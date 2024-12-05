#fig. S12
#devtools::install_github("cozygene/FEAST")
# Feast emerged species ------------------------------
pacman::p_load(dplyr,tidyverse,ggplot2,ggsci,agricolae,FEAST)
otu <- read.delim('~/Desktop/molt github/Fig. S10/ASV_table_norm feast.txt',row.names = 1,check.names = F)
group <- read.delim('~/Desktop/molt github/molt plot/groupMOLT.txt')
colSums(otu)
##Foregut
df1 <- read.csv('~/Desktop/molt github/Fig. S10/Emerge/Foregut48h.csv',check.names = F,row.names = 1)
mean(colSums(df1)/225.07)
metadata <- Load_metadata(metadata_path = "~/Desktop/molt github/Fig. S10/Emerge/48hgroup_qc.txt")
otus <- Load_CountMatrix(CountMatrix_path = "~/Desktop/molt github/Fig. S10/ASV_table_norm feast.txt")
otus1 <-otus[(rownames(otus) %in% rownames(metadata)),(colnames(otus) %in% rownames(df1))]
FEAST_output <- FEAST(C = otus1, metadata = metadata, different_sources_flag =1,COVERAGE = 22507,
											dir_path = "~/Desktop/",outfile="2222",EM_iterations = 1000)  
prop1 <- round(FEAST_output$data_prop,10);prop1
##Midgut
df1 <- read.csv('~/Desktop/molt github/Fig. S10/Emerge/Midgut48h.csv',check.names = F,row.names = 1)
mean(colSums(df1)/225.07)
metadata <- Load_metadata(metadata_path = "~/Desktop/molt github/Fig. S10/Emerge/48hgroup_mc.txt")
otus <- Load_CountMatrix(CountMatrix_path = "~/Desktop/molt github/Fig. S10/ASV_table_norm feast.txt")
otus1 <-otus[(rownames(otus) %in% rownames(metadata)),(colnames(otus) %in% rownames(df1))]
FEAST_output <- FEAST(C = otus1, metadata = metadata, different_sources_flag =1,COVERAGE = 22507,
											dir_path = "~/Desktop/",outfile="2222",EM_iterations = 1000)  
prop2 <- round(FEAST_output$data_prop,10);prop2
##Hindgut
df1 <- read.csv('~/Desktop/molt github/Fig. S10/Emerge/Hindgut48h.csv',check.names = F,row.names = 1)
mean(colSums(df1)/225.07)
metadata <- Load_metadata(metadata_path = "~/Desktop/molt github/Fig. S10/Emerge/48hgroup_hc.txt")
otus <- Load_CountMatrix(CountMatrix_path = "~/Desktop/molt github/Fig. S10/ASV_table_norm feast.txt")
otus1 <-otus[(rownames(otus) %in% rownames(metadata)),(colnames(otus) %in% rownames(df1))]
FEAST_output <- FEAST(C = otus1, metadata = metadata, different_sources_flag =1,COVERAGE = 22507,
											dir_path = "~/Desktop/",outfile="2222",EM_iterations = 1000)  
prop3 <- round(FEAST_output$data_prop,10);prop3
#merge
prop1$tissue <- 'Foregut';prop1$source <- rownames(prop1)
prop2$tissue <- 'Midgut';prop2$source <- rownames(prop2)
prop3$tissue <- 'Hindgut';prop3$source <- rownames(prop3)
prop <- rbind(prop1,prop2,prop3)
#write.csv(prop,'~/Desktop/molt github/Fig. S10/Emerge/emerge source11.csv',row.names = F)

##fig. S10A plot feast emerged species----
library(dplyr);library(tidyverse);library(tidytext)
source('~/Desktop/molt github/molt plot/theme_linweichuan.R');
ds <- read.csv('~/Desktop/molt github/Fig. S10/Emerge/emerged source11.csv')
#unknown,others;add facet_wrap annotation
ds1 <- ds %>%
	mutate(source = as.factor(source),
				 name = reorder_within(source, pred_emnoise_all,tissue),
				 name = fct_relevel(name, paste0("Unknown___", unique(source))))
ds1$tissue <- factor(ds1$tissue,levels = c('Foregut','Midgut','Hindgut'))
pd <- ggplot(ds1,aes(y=name,x=pred_emnoise_all*100))+
	geom_col(width = 0.8,fill='#1F78B4')+
	scale_x_continuous(expand = c(0,2.5),limits = c(0,100),breaks=seq(0,100,25))+
	geom_text(data=filter(ds1,pred_emnoise_all>0.0),
						aes(label=round(pred_emnoise_all*100,2)),
						hjust=-0.1,size=2.8)+labs(y='Source',x='Contribution (%)')+
	scale_y_reordered() + labs(title = 'Emerged ASVs')+
	facet_wrap(~tissue,nrow = 3,scales = 'free_y',strip.position = 'right')+
	geom_text(data = data.frame(x = c(70,70,70), y = c(3,3,3), 
															lab = c( "Foregut (1.75%)",
																			 "Midgut (1.44%)",
																			 "Hindgut (4.73%)"),
															tissue= factor(c( "Foregut","Midgut","Hindgut"))),
						aes( x, y, label = lab ),vjust = 1,size=2.8,color='black')+
	theme_linweichuan()+#coord_trans(y="reverse")+
	theme(strip.background = element_blank(),strip.text = element_blank(),
				axis.text = element_text(size=8),axis.title.x = element_text(size=8),
				axis.title.y = element_text(size=8))
pd 

#Feast enriching species ------------------------------
pacman::p_load(dplyr,tidyverse,ggplot2,ggsci,agricolae,FEAST)
otu <- read.delim('~/Desktop/molt github/Fig. S10/ASV_table_norm feast.txt',row.names = 1,check.names = F)
group <- read.delim('~/Desktop/molt github/molt plot/groupMOLT.txt')
####Foregut
df <- read.csv('~/Desktop/molt github/Fig. S10/Enriched/selectForegut48h_0h.csv')
df <- filter(df,sig=='Enriched')
df1 <- read.csv('~/Desktop/molt github/Fig. S10/Emerge/Foregut48h.csv',check.names = F,row.names = 1)
df2 <- setdiff(as.character(df$ASV_ID),rownames(df1))
metadata <- Load_metadata('~/Desktop/molt github/Fig. S10/Enriched/48hgroup_qc.txt')
otu1 <- t(otu[rownames(otu)%in% as.character(df2),
							colnames(otu)%in% rownames(metadata)])
FEAST_output <- FEAST(C = otu1, metadata = metadata, different_sources_flag =1,COVERAGE = 22507,
											dir_path = '~/Desktop/',outfile='enriched',EM_iterations = 1000)     
prop1 <- round(FEAST_output$data_prop,10);prop1
####Midgut
df <- read.csv('~/Desktop/molt github/Fig. S10/Enriched/selectMidgut48h_0h.csv')
df <- filter(df,sig=='Enriched')
df1 <- read.csv('~/Desktop/molt github/Fig. S10/Emerge/Midgut48h.csv',check.names = F,row.names = 1)
df2 <- setdiff(as.character(df$ASV_ID),rownames(df1))
metadata <- Load_metadata('~/Desktop/molt github/Fig. S10/Enriched/48hgroup_mc.txt')
otu1 <- t(otu[rownames(otu)%in% as.character(df2),
							colnames(otu)%in% rownames(metadata)])
FEAST_output <- FEAST(C = otu1, metadata = metadata, different_sources_flag =1,COVERAGE = 22507,
											dir_path = '~/Desktop/',outfile='enriched',EM_iterations = 1000)     
prop2 <- round(FEAST_output$data_prop,10);prop2
####Hindgut
df <- read.csv('~/Desktop/molt github/Fig. S10/Enriched/selectHindgut48h_0h.csv')
df <- filter(df,sig=='Enriched')
df1 <- read.csv('~/Desktop/molt github/Fig. S10/Emerge/Hindgut48h.csv',check.names = F,row.names = 1)
df2 <- setdiff(as.character(df$ASV_ID),rownames(df1))
metadata <- Load_metadata('~/Desktop/molt github/Fig. S10/Enriched/48hgroup_hc.txt')
otu1 <- t(otu[rownames(otu)%in% as.character(df2),
							colnames(otu)%in% rownames(metadata)])
FEAST_output <- FEAST(C = otu1, metadata = metadata, different_sources_flag =1,COVERAGE = 22507,
											dir_path = '~/Desktop/',outfile='enriched',EM_iterations = 1000)     
prop3 <- round(FEAST_output$data_prop,10);prop3
#merge
prop1$tissue <- 'Foregut';prop1$source <- rownames(prop1)
prop2$tissue <- 'Midgut';prop2$source <- rownames(prop2)
prop3$tissue <- 'Hindgut';prop3$source <- rownames(prop3)
prop <- rbind(prop1,prop2,prop3)
#write.csv(prop,'~/Desktop/molt github/Fig. S10/Enriched/enriched source11.csv',row.names = F)

##fig. S10B plot Feast enriching species----
library(dplyr);library(tidyverse);library(tidytext)
pd <- read.csv('~/Desktop/molt github/Fig. S10/Enriched/enriched source11.csv')
source('~/Desktop/molt github/molt plot/theme_linweichuan.R')
#unknown,others;add facet_wrap annotation
pd1 <- pd %>%
	mutate(source = as.factor(source),
				 name = reorder_within(source, pred_emnoise_all,tissue),
				 name = fct_relevel(name, paste0("Unknown___", unique(source))))
pd1$tissue <- factor(pd1$tissue,levels = c('Foregut','Midgut','Hindgut'))
pe <- ggplot(pd1,aes(y=name,x=pred_emnoise_all*100))+
	geom_col(width = 0.8,fill='#1F78B4')+
	scale_x_continuous(expand = c(0,2.5),limits = c(0,100),breaks=seq(0,100,25))+
	geom_text(data=filter(pd1,pred_emnoise_all>0.0),
						aes(label=round(pred_emnoise_all*100,2)),
						hjust=-0.1,size=2.8)+labs(y='Source',x='Contribution (%)')+
	scale_y_reordered() + labs(title='Enriched ASVs')+
	facet_wrap(~tissue,nrow = 3,scales = 'free_y',strip.position = 'right')+
	geom_text(data = data.frame(x = c(70,70,70), y = c(3,3,3), 
															lab = c( "Foregut (4.80%)",
																			 "Midgut (26.32%)",
																			 "Hindgut (14.81%)"),
															tissue= factor(c( "Foregut","Midgut","Hindgut"))),
						aes( x, y, label = lab ),vjust = 1,size=2.8,color='black')+
	theme_linweichuan()+
	theme(strip.background = element_blank(),strip.text = element_blank(),
				axis.text = element_text(size=8),axis.title.x = element_text(size=8),
				axis.title.y = element_text(size=8))
pe


library(ggpubr)
p <- ggarrange(pd,pe,nrow = 1,ncol = 2,
							 labels = c('ABC','DEF'),font.label = list(size = 10,face = "plain"))
p
#ggsave('~/Desktop/psource.pdf',p,width = 6.45,height = 5.00,dpi = 600,units = 'in')

#Finally, these plots are arranged into fig. S10 in AI software and further modify legends.


