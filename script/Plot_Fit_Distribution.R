# Title     :
# Objective :
# Created by: yiquan
# Created on: 9/7/21
library(ggplot2)
library(readr)


table_data <- read_csv('result/Mos99_fit.csv')

r1_p <- ggplot(table_data, aes(x=Rep1_norm_fitness,fill=mutation_type)) +
    geom_histogram(binwidth=.2,colour="black", fill="white")+
    facet_grid(mutation_type ~ .)

ggsave('graph/R1_fit_distribution.png',r1_p,height=4,width=2,bg='white')

