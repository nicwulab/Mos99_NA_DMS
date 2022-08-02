#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(plyr)
library(dplyr)
library(qualpalr)
library(gridExtra)
library(ggbeeswarm)
library(sinaplot)
library(ggforce)
require(cowplot)

plot_dist_vs_fit <- function(df, graphname){
  textsize <- 7
  palette <- qualpal(n = 5, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex
  p <- ggplot(df,aes(x=dist_to_SIA, y=fit, color=type)) +
         geom_point(size=0.5, pch=16, alpha=0.8) +
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_blank(),
               legend.key.height = unit(0.15, 'in'),
               legend.key.width = unit(0.05, 'in'),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='right') +
         guides(colour = guide_legend(override.aes = list(size=0.5))) +
         scale_color_manual(values=palette,drop=FALSE) +
         labs(x=expression(bold("Distance from active site (Å)")),y=expression(bold('mutational tolerance')))
  ggsave(graphname, p, height=2, width=3, dpi=600)
  }

fit  <- read_tsv('result/position_type_vs_fit.tsv')
dist <- read_tsv('result/Dist_to_active_site.tsv')
df  <- merge(x=fit, y=dist, by='pos', all=TRUE)
plot_dist_vs_fit(df, 'graph/fit_vs_dist.png')

print (paste('Pearson correlation between RSA and fitness: ', cor(df$dist, df$fit, method="pearson"), sep=''))
print (paste('Spearman correlation between RSA and fitness: ', cor(df$dist, df$fit, method="spearman"), sep=''))
