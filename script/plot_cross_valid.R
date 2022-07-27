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
require(cowplot)

plot_cross_validate <- function(t_join, graphname){
  textsize <- 7
  p <- ggplot(t_join,aes(x=fit_dms, y=fit)) +
         geom_point(size=1) +
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
     labs(x=expression(bold('fitness (this study)')),y=expression(bold('fitness (Wang et al. 2021)')))
  ggsave(graphname, p, height=1.7, width=1.7)
  }

t_all  <- read_csv('result/Mos99_fit.csv') %>%
            mutate(fit_dms=Fitness) %>%
            rename(mut=Mutation)
t_epi  <- read_tsv('data/Mos99_single_mut_Wang_et_al.tsv') %>%
            mutate(fit=log10(fit))
t_join <- inner_join(t_all, t_epi)
print (t_join$mut)
print (paste('correlation of log10 fitness between replicates: ', cor(t_join$fit_dms, t_join$fit), sep=''))
plot_cross_validate(t_join, 'graph/DMS_cross_validate.png')
