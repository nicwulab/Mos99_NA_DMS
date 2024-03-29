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
  palette  <- c(brewer.pal(9,"Set1"))
  p <- ggplot(t,aes(x=msa_transformer, y=Fitness)) +
         geom_point(size=0.3,alpha=0.3, pch=16, color='black') +
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
         scale_x_continuous(breaks=c(-15,-10,-5,0,5), labels=c('-15','-10','-5','0','5'), lim=c(-19,8)) +
         scale_y_continuous(breaks=c(-2,-1,0,1), labels=c('-2','-1','0','1'), lim=c(-2.9,1)) +
         labs(x=expression(bold('MSA transformer\nprediction')),y=expression(bold('fitness')))
  ggsave(graphname, p, height=1.5, width=1.5, dpi=600)
  }

plot_ddG <- function(t, graphname, form, w, h){
  textsize <- 7
  if (form=='tetramer'){
    t <- rename(t, ddG=ddG_tetramer)
    }
  if (form=='monomer'){
    t <- rename(t, ddG=ddG_monomer)
    }
  p <- ggplot(t,aes(x=ddG, y=Fitness)) +
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
  ggsave(graphname, p, height=h, width=w, dpi=600)
  }

ceiling_1 <- function(x){
  if (x>30){return (30)}
  else {return (x)}
  }

ceiling_2 <- function(x){
  if (x>8){return (8)}
  else {return (x)}
  }

natural_muts <- read_tsv('result/N2_mutation_freq.tsv')
natural_muts <- unique(natural_muts$mut)
fit  <- read_csv('result/Mos99_fit.csv') %>%
          rename(mut=Mutation)
pred <- read_csv('data/foldx_msa_transformer.csv') %>%
          inner_join(fit, by='mut') %>%
          mutate(pos=as.numeric(str_sub(mut, 2, -2))) %>%
          mutate(natural=ifelse(mut %in% natural_muts, 'yes', 'no'))
pred_ceil <- pred %>%
          mutate(ddG_tetramer=mapply(ceiling_1, ddG_tetramer)) %>%
          mutate(ddG_monomer=mapply(ceiling_2, ddG_monomer))
print (paste('# of natural mutations: ', length(filter(pred,natural=='yes')$mut), sep=''))
print (paste('# of unnatural mutations: ', length(filter(pred,natural=='no')$mut), sep=''))
plot_MSA_transformer(pred, 'graph/fit_vs_MSA_transformer_all.png')
plot_MSA_transformer(filter(pred,natural=='yes'), 'graph/fit_vs_MSA_transformer_natural.png')
plot_MSA_transformer(filter(pred,natural=='no'), 'graph/fit_vs_MSA_transformer_unnatural.png')
plot_ddG(pred_ceil, 'graph/fit_vs_foldX_ddG_tetramer.png', 'tetramer', 1.5, 1.5)
plot_ddG(pred_ceil, 'graph/fit_vs_foldX_ddG_monomer.png', 'monomer', 2, 2)
print (paste('Pearson correlation between fitness and MSA (Overall): ',
       cor(pred$msa_transformer, pred$Fitness, method="pearson",use="complete.obs"), sep=''))
print (paste('Spearman correlation between fitness and MSA (Overall): ',
       cor(pred$msa_transformer, pred$Fitness, method="spearman", use="complete.obs"), sep=''))
print (paste('Pearson correlation between fitness and MSA (natural): ',
       cor(filter(pred,natural=='yes')$msa_transformer, filter(pred,natural=='yes')$Fitness, method="pearson",use="complete.obs"), sep=''))
print (paste('Spearman correlation between fitness and MSA (natural): ',
       cor(filter(pred,natural=='yes')$msa_transformer, filter(pred,natural=='yes')$Fitness, method="spearman", use="complete.obs"), sep=''))
print (paste('Pearson correlation between fitness and MSA (unnatural): ',
       cor(filter(pred,natural=='no')$msa_transformer, filter(pred,natural=='no')$Fitness, method="pearson",use="complete.obs"), sep=''))
print (paste('Spearman correlation between fitness and MSA (unnatural): ',
       cor(filter(pred,natural=='no')$msa_transformer, filter(pred,natural=='no')$Fitness, method="spearman", use="complete.obs"), sep=''))
print (paste('Pearson correlation between fitness and ddG_tetramer: ',
       cor(pred$ddG_tetramer, pred$Fitness, method="pearson",use="complete.obs"), sep=''))
print (paste('Spearman correlation between fitness and ddG_tetramer: ',
       cor(pred$ddG_tetramer, pred$Fitness, method="spearman", use="complete.obs"), sep=''))
print (paste('Pearson correlation between fitness and ddG_monomer: ',
       cor(pred$ddG_monomer, pred$Fitness, method="pearson",use="complete.obs"), sep=''))
print (paste('Spearman correlation between fitness and ddG_monomer: ',
       cor(pred$ddG_monomer, pred$Fitness, method="spearman", use="complete.obs"), sep=''))
