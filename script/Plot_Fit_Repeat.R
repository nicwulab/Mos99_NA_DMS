# Title     :
# Objective :
# Created by: yiquan
# Created on: 9/7/21
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
require(cowplot)


table_data <- read_csv('result/Mos99_fit.csv')

PlotFit_Rep <- function(table_data,rep1_fit,rep2_fit,title){
  textsize=7
  colorscale <- c(brewer.pal(9,"Set1"))
  R1fit = table_data[,rep1_fit]
  R2fit = table_data[,rep2_fit]
  print (paste("Pearson Cor:", cor(R1fit,R2fit,use="complete.obs"),sep=' '))
  p <- ggplot(table_data,aes(x=Rep1_norm_fitness,y=Rep2_norm_fitness,color=mutation_type)) +
	  xlim(-2.8,2)+ylim(-2.8,2)+
	 geom_point(size=0.3) +
         theme_cowplot(12) +
         scale_color_manual(values=c('gray', colorscale)) +
	 theme(plot.title=element_text(size=textsize,face="bold",hjust = 0.5),
	       axis.title=element_text(size=textsize,face="bold"),
	       axis.text=element_text(size=textsize,face="bold"),
	       legend.title=element_blank(),
	       legend.text=element_text(size=textsize,face="bold"),
	       legend.position='none') +
	 xlab(bquote(bold('Replicate 1'))) +
	 ylab(bquote(bold('Replicate 2'))) +
         ggtitle(title)
  return (p)
  }
wt_p <- PlotFit_Rep(table_data,"Rep1_norm_fitness","Rep2_norm_fitness","Mos99")
ggsave('graph/WT_fit_repeat.png',wt_p,height=4,width=4,bg='white')
