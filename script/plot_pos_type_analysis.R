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

plot_RSA_vs_fit <- function(fit_table, graphname){
  textsize <- 7
  palette  <- c(brewer.pal(9,"Set1"))
  #palette <- qualpal(n = 5, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex
  p <- ggplot(fit_table,aes(x=RSA_tetramer*100, y=fit, color=type)) +
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
         labs(x=expression(bold("RSA (%)")),y=expression(bold('mutational tolerance')))
  ggsave(graphname, p, height=2, width=3, dpi=600)
  }

plot_position_type <- function(data_table, graphname, ylab,violin){
  summary <- ddply(data_table,c("type"), summarise, mean = mean(parameter), sd = sd(parameter))
  textsize <- 7
  if (violin=='yes'){
    p <-  ggplot() +
            geom_violin(data=data_table,aes(x=type,y=parameter),size=0.3, scale='area', width=1)
    }
  else{
    p <-  ggplot()
    }
  p <- p +
          geom_boxplot(data=data_table,aes(x=type,y=parameter),width=0.2, outlier.shape=NA, size=0.3) +
          geom_sina(data=data_table,aes(x=type,y=parameter),
                    pch=16, size=0.1,method="counts", bin_limit=0.4, scale="width", maxwidth=0.5, color='gray40', alpha=0.5) +
          theme_cowplot(12) +
          theme(plot.background = element_rect(fill = "white"),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title=element_text(size=textsize,face="bold"),
                legend.title=element_blank(),
                legend.key.height = unit(0.15, 'in'),
                legend.key.width = unit(0.05, 'in'),
                legend.text=element_text(size=textsize,face="bold"),
                legend.position='right') +
          guides(colour = guide_legend(override.aes = list(size=0.5))) +
          xlab("") +
          ylab(ylab)
  ggsave(graphname, p, width=2, height=2.2, dpi=600)
  }

set.seed(26)
position_types <- c('buried','interface','active site','antigenic','exposed')
fit_table <- read_tsv('result/position_type_vs_fit.tsv') %>%
               mutate(type=factor(type, level=position_types))
plot_RSA_vs_fit(fit_table, 'graph/fit_vs_RSA.png')
print (paste('Pearson correlation between RSA and fitness: ', cor(fit_table$RSA_tetramer, fit_table$fit, method="pearson"), sep=''))
print (paste('Spearman correlation between RSA and fitness: ', cor(fit_table$RSA_tetramer, fit_table$fit, method="spearman"), sep=''))

RSA_table <- fit_table %>%
               mutate(RSA_tetramer=RSA_tetramer*100) %>%
               select(type, RSA_tetramer) %>%
               rename(parameter=RSA_tetramer)
FIT_table <- fit_table %>%
               select(type, fit) %>%
               rename(parameter=fit)
plot_position_type(RSA_table, 'graph/position_type_vs_RSA.png', expression(bold("RSA (%)")),'no')
plot_position_type(FIT_table, 'graph/position_type_vs_fit.png', expression(bold("mutational tolerance")),'yes')

for (pos_type1 in position_types){
  for (pos_type2 in position_types){
    p_value <- t.test(filter(fit_table,type==pos_type1)$fit, filter(fit_table,type==pos_type2)$fit)$p.value
    #p_value <- t.test(filter(fit_table,type==pos_type1)$RSA_tetramer, filter(fit_table,type==pos_type2)$RSA_tetramer)$p.value
    print (paste(pos_type1, 'vs', pos_type2, ':', p_value))
    }
  }
