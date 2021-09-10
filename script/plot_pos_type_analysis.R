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
library(gridExtra)
library(ggbeeswarm)
library(sinaplot)
library(ggforce)
require(cowplot)

plot_RSA_vs_fit <- function(fit_table, graphname){
  textsize <- 7
  palette  <- c(brewer.pal(9,"Set1"))
  p <- ggplot(fit_table,aes(x=RSA_tetramer, y=fit, color=type)) +
         geom_point(size=0.2) +
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
         labs(x=expression(bold("RSA ("*Å^"2"*")")),y=expression(bold('fitness')))
  ggsave(graphname, p, height=2, width=3.5, dpi=300)
  }

plot_position_type <- function(data_table, graphname, ylab){
  summary <- ddply(data_table,c("type"), summarise, mean = mean(parameter), sd = sd(parameter))
  textsize <- 7
  palette  <- c(brewer.pal(9,"Set1"))
  p <-  ggplot() +
          geom_sina(data=data_table,aes(x=type,y=parameter,color=type), size=0.1,
                    method="counts", bin_limit=0.4, scale="width", maxwidth=0.5) +
          #geom_violin(data=data_table,aes(x=type,y=parameter,fill=type), size=0.1, scale='area') +
          #geom_beeswarm(data=data_table,aes(x=type,y=parameter,color=type), size=0.1) +
          geom_errorbar(data=summary, mapping=aes(x=type, ymin=mean-sd, ymax=mean+sd),
                        size=0.3, color="black", width=0.2, position=position_dodge(1)) +
          geom_point(data=summary, mapping = aes(x=type, y=mean),
                     size=4, color="black", shape="-", position=position_dodge(1)) +
          theme_cowplot(12) +
          theme(plot.background = element_rect(fill = "white"),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title=element_text(size=textsize,face="bold"),
                #axis.line = element_line(colour = 'black', size = 0),
                #panel.border = element_rect(colour = "black", fill=NA, size=1), 
                legend.title=element_blank(),
                legend.key.height = unit(0.15, 'in'),
                legend.key.width = unit(0.05, 'in'),
                legend.text=element_text(size=textsize,face="bold"),
                legend.position='right') +
          guides(colour = guide_legend(override.aes = list(size=0.5))) +
          scale_color_manual(values=palette,drop=FALSE) +
          xlab("") +
          ylab(ylab)
  ggsave(graphname, p, width=3.5, height=3.5, dpi=300)
  }

set.seed(26)
position_types <- c('active site','antigenic site',
                    'buried','exposed','interface')
fit_table <- read_tsv('result/position_type_vs_fit.tsv') %>%
               mutate(type=factor(type, level=position_types))
plot_RSA_vs_fit(fit_table, 'graph/RSA_vs_fit.png')
RSA_table <- fit_table %>%
               select(type, RSA_tetramer) %>%
               rename(parameter=RSA_tetramer)
FIT_table <- fit_table %>%
               select(type, fit) %>%
               rename(parameter=fit)
plot_position_type(RSA_table, 'graph/position_type_vs_RSA.png', expression(bold("RSA ("*Å^"2"*")")))
plot_position_type(FIT_table, 'graph/position_type_vs_fit.png', expression(bold("fitness")))
