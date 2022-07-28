#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(qualpalr)
library(gridExtra)
library(viridis)
require(cowplot)

plot_fitness_hist <- function(t, graphname){
  textsize <- 7
  palette <- qualpal(n = 3, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex
  p <- ggplot(t,aes(x=Fitness, group=stab, fill=stab)) +
         geom_histogram(binwidth=0.05) +
         scale_fill_manual(values=palette,drop=FALSE) +
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title.x=element_text(size=textsize,face="bold",vjust=-3,hjust=0.5),
               axis.title.y=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               strip.text = element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_blank(),
               legend.text=element_text(size=textsize,face="bold"),
               legend.justification='center',
               legend.position='top') +
         scale_x_continuous(breaks=c(-3,-2,-1,0,1), labels=c('-3','-2','-1','0','1'), limits=c(-3,1)) +
         labs(x=expression(bold('fitness')),y=expression(bold('# of mutations'))) +
         facet_grid(natural ~ ., scales = "free")
  ggsave(graphname, p, height=1.7, width=2, dpi=600)
  }

fit  <- read_csv('result/Mos99_fit.csv') %>%
          rename(mut=Mutation)
pred <- read_csv('data/foldx_msa_transformer.csv') %>%
          inner_join(fit, by='mut') %>%
          mutate(pos=as.numeric(str_sub(mut, 2, -2))) %>%
          mutate(stab=ifelse(ddG_foldx <= 0, 'stable', 'unstabile')) %>%
          mutate(natural=ifelse(!is.na(msa_transformer), 'observed', 'unobserved'))
obs <- filter(pred, natural=='observed') %>% nrow(.)
unobs <- filter(pred, natural=='unobserved') %>% nrow(.)
highfit_obs   <- filter(pred, Fitness >= 0) %>% filter(natural=='observed') %>% nrow(.)
highfit_unobs <- filter(pred, Fitness >= 0) %>% filter(natural=='unobserved') %>% nrow(.)
print (paste('% of observed being high fit:', round(highfit_obs/obs, 2), '(',highfit_obs,'/',obs,')', sep=''))
print (paste('% of unobserved being high fit:', round(highfit_unobs/unobs, 2), '(',highfit_unobs,'/',unobs,')', sep=''))
plot_fitness_hist(pred, 'graph/fit_dist_nat_vs_unnat.png')
