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
require(cowplot)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
  }

plot_MSA_transformer <- function(t, graphname){
  textsize <- 7
  p <- ggplot(t,aes(x=msa_transformer, y=Fitness)) +
         geom_point(size=0.5, color='black',alpha=0.3, pch=16) +
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title.x=element_text(size=textsize,face="bold",vjust=-3,hjust=0.5),
               axis.title.y=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='none') +
         scale_x_continuous(breaks=c(-15,-10,-5,0,5), labels=c('-15','-10','-5','0','5')) +
         labs(x=expression(bold('MSA transformer\nprediction')),y=expression(bold('fitness')))
  ggsave(graphname, p, height=1.5, width=1.5, dpi=600)
  }

plot_ddG <- function(t, graphname){
  textsize <- 7
  p <- ggplot(t,aes(x=ddG_foldx, y=Fitness)) +
         geom_point(size=0.5, color='black',alpha=0.2, pch=16) +
         #geom_smooth(method="loess", se=FALSE, fullrange=FALSE, level=0.95, color='blue' ,size=0.3) +
         #stat_smooth(method="loess", se=FALSE, span=0.1, color='blue' ,size=0.3) +
         stat_smooth(method="loess", se=FALSE, span=0.8, degree=1, color='blue' ,size=0.3) +
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title.x=element_text(size=textsize,face="bold",vjust=-3,hjust=0.5),
               axis.title.y=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='none') +
         #scale_x_continuous(breaks=c(-15,-10,-5,0,5), labels=c('-15','-10','-5','0','5')) +
         labs(x=expression(bold('FoldX ddG\nprediction (kcal/mol)')),y=expression(bold('fitness')))
  ggsave(graphname, p, height=1.5, width=1.5, dpi=600)
  }

ceiling <- function(x){
  if (x>30){return (30)}
  else {return (x)}
  }

fit  <- read_csv('result/Mos99_fit.csv') %>%
          rename(mut=Mutation)
pred <- read_csv('data/foldx_msa_transformer.csv') %>%
          inner_join(fit, by='mut') %>%
          mutate(pos=as.numeric(str_sub(mut, 2, -2))) %>%
          mutate(ddG_foldx=mapply(ceiling, ddG_foldx))
plot_MSA_transformer(pred, 'graph/fit_vs_MSA_transformer.png')
plot_ddG(pred, 'graph/fit_vs_foldX_ddG.png')
print (paste('Pearson correlation between fitness and MSA: ',
       cor(pred$msa_transformer, pred$Fitness, method="pearson",use="complete.obs"), sep=''))
print (paste('Spearman correlation between fitness and MSA: ',
       cor(pred$msa_transformer, pred$Fitness, method="spearman", use="complete.obs"), sep=''))
print (paste('Pearson correlation between fitness and foldX_ddG: ',
       cor(pred$ddG_foldx, pred$Fitness, method="pearson",use="complete.obs"), sep=''))
print (paste('Spearman correlation between fitness and foldX_ddG: ',
       cor(pred$ddG_foldx, pred$Fitness, method="spearman", use="complete.obs"), sep=''))
