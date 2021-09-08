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
library(viridis)
library(MASS)
require(cowplot)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
  }

plot_hist <- function(t, main_title){
  textsize <- 7
  p <- ggplot(t,aes(x=fit)) +
         geom_histogram(binwidth=0.1) +
         theme_cowplot(12) +
         theme(plot.title=element_text(size=textsize,face="bold", hjust=0.5),
               axis.title=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_blank(),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='right') +
     ggtitle(main_title) +
     scale_fill_manual(values=c('black'),drop=FALSE) +
     labs(y=expression(bold('count')),x=expression(bold('fitness'))) +
     coord_cartesian(xlim=c(-2.5,1))
  return (p)
  }

plot_replicate_cor <- function(t_all, graphname){
  textsize <- 7
  t_all$density <- get_density(t_all$Rep1_norm_fitness, t_all$Rep2_norm_fitness, n = 100)
  p <- ggplot(t_all,aes(x=Rep1_norm_fitness, y=Rep2_norm_fitness, color=density)) +
         geom_point(size=0.1) +
         scale_color_viridis() +
         #geom_point(size=0.1, alpha=0.1) +
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='right') +
     scale_fill_manual(values=c('black'),drop=FALSE) +
     labs(x=expression(bold('fitness (replicate 1)')),y=expression(bold('fitness (replicate 2)')))
  ggsave(graphname, p, height=2, width=2)
  }

t_all <- read_csv('result/Mos99_fit.csv') %>%
	   mutate(fit=fitness) %>%
	   filter(!grepl('Amp', Mutation)) %>%
	   filter(!grepl('WT', Mutation)) %>%
           filter(!grepl('-', Mutation))
           #filter(Amplicon=='Amp6')
#write.table(t_all, file='test.tsv', sep="\t", quote=F)

t_non <- t_all %>%
           filter(grepl('_', Mutation)) %>%
           filter(!grepl('silent',Mutation))
t_sil <- t_all %>%     
           filter(grepl('silent',Mutation))
t_mis <- t_all %>%
           filter(!grepl('_', Mutation)) %>%
           filter(!grepl('silent',Mutation))
p_sil <- plot_hist(t_sil, 'silent')
p_non <- plot_hist(t_non, 'nonsense')
p_mis <- plot_hist(t_mis, 'missense')

#print (filter(t_sil,fit<1)$Mutation)

#p <- grid.arrange(p_mis_all, p_mis_WTplasmid0, p_non_all,p_non_WTplasmid0,nrow=4)
p <- grid.arrange(p_sil, p_non, p_mis, nrow=3)
#ggsave('graph/Mos99_mis_and_nonsense.png', p, height=6, width=1.8)
ggsave('graph/Mos99_sil_non_mis.png', p, height=3.5, width=2)

print (paste('correlation of log10 fitness between replicates: ', cor(t_all$Rep1_norm_fitness, t_all$Rep2_norm_fitness), sep=''))
print (paste('total number of missense mutants: ', length(t_mis$Mutation), sep=''))
plot_replicate_cor(t_all, 'graph/Mos99_rep_cor.png')
