#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(gridExtra)
library(stringr)
require(cowplot)

plot_fitness_heatmap <- function(fitness_table, WTresibox, start_resi, end_resi){
  textsize <- 5
  fitness_table <- fitness_table %>%
                     filter(Pos >= start_resi & Pos <= end_resi)
  WTresibox     <- WTresibox %>%
                     filter(Pos >= start_resi & Pos <= end_resi) %>%
                     mutate(x=x-min(x)+1)
  p <-  ggplot() +
          geom_tile(data=fitness_table,aes(x=resi,y=aa,fill=Fitness)) +
          scale_fill_gradientn(colours=c("blue","blue", "white", "red"),
                limits=c(-2.8,1),
                values=rescale(c(-2.8,-1, 0, 1)),
                breaks=c(-2,-1,0,1),
                labels=c('-2','-1','0','1'),
                guide="colorbar",
                na.value="grey") +
          theme_cowplot(12) +
          theme(plot.background = element_rect(fill = "white"),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title=element_text(size=7,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1)) +
          guides(fill = guide_colorbar(title.theme=element_text(size=7,face="bold",colour='black',hjust=0.5),
                                       label.theme=element_text(size=7,face="bold",colour='black'),
                                       frame.colour="black",
                                       frame.linewidth = 1,
                                       ticks = TRUE,
                                       ticks.colour = "black",
                                       barwidth = 0.5, barheight = 6, title="fitness")) +
          geom_point(data=WTresibox, aes(x=x, y=y), color='black', size=0.2) +
          xlab("") +
          ylab("amino acid")
  }

aa_level <- rev(c('E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','_'))

fitness_table <- read_csv('result/Mos99_fit.csv') %>%
  filter(!grepl('silent',Mutation)) %>%
  mutate(resi=str_sub(Mutation,1,-2)) %>%
  mutate(aa=str_sub(Mutation,-1,-1)) %>%
  filter(aa %in% aa_level) %>%
  mutate(aa=factor(aa,levels=aa_level)) %>%
  mutate(resi=factor(resi,levels=unique(resi))) %>%
  complete(resi, aa) %>%
  mutate(Pos=str_sub(resi,2,-1)) %>%
  mutate(Pos=as.numeric(as.character(Pos))) %>%
  arrange(Pos) %>%
  mutate(Fitness=case_when(str_sub(resi,1,1)==aa ~ 0, TRUE ~ Fitness)) %>%
  mutate(Mutation=paste(resi,aa,sep='')) %>%
  select(Mutation, resi, Pos, aa, Fitness)

WTresibox  <- fitness_table %>%
  select(resi,Pos) %>%
  unique() %>%
  mutate(WT_resi=str_sub(resi,1,1)) %>%
  mutate(x=seq(82-81,465-81)) %>%
  mutate(y=match(WT_resi,aa_level)) %>%
  select(resi,WT_resi,Pos,x, y)

p1 <- plot_fitness_heatmap(fitness_table, WTresibox, 82, 177)
p2 <- plot_fitness_heatmap(fitness_table, WTresibox, 178, 273)
p3 <- plot_fitness_heatmap(fitness_table, WTresibox, 274, 369)
p4 <- plot_fitness_heatmap(fitness_table, WTresibox, 370, 465)
p <- grid.arrange(p1, p2, p3, p4, nrow=4)
ggsave('graph/Mos99_fit_heatmap.png',p,width=6.5, height=7.5, dpi=300)
