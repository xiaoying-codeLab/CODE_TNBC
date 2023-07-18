# Settings
# Objective: To visualize the results of immune infiltration
# Data Source:
# Save data:

rm(list=ls())
options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(ggsci)
library(ggpubr)
library(patchwork)


#1、load
load('05-doEstimate/out-data/estimate_results.Rdata')
dim(scores)
head(scores)
load('04-1-doModel-lasso/out-data//basel-risk.Rdata')
rt =new_dat
dim(rt)
rt[1:4, 1:4]
table(rt$Risk_level)

# 2. 
intersect(rownames(rt), rownames(scores))
scores = scores[intersect(rownames(rt), rownames(scores)),]
rt = rt[intersect(rownames(rt), rownames(scores)),]
identical(rownames(rt), rownames(scores))
data.score = cbind(scores,rt)
colnames(data.score)
data.score$group <-  rt$Risk_level
head(data.score)


# 
data.score$purity = cos(0.61 + 0.00015 * data.score$ESTIMATEScore)

b = ggplot(dat = data.score, aes(group,purity))+
  geom_boxplot(aes(fill = group))+
  scale_fill_nejm()+
  stat_compare_means(method='kruskal.test') + labs(x='Score', title = 'purity')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        legend.position = '')
b



b1 = ggplot(dat = data.score, aes(group,StromalScore))+
  geom_boxplot(aes(fill = group))+
  scale_fill_nejm()+
  stat_compare_means(method='kruskal.test') + labs(x='Score', title = 'StromalScore')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        legend.position = '')
b1

b2 = ggplot(dat = data.score,aes(group,ImmuneScore ))+
  geom_boxplot(aes(fill = group))+
  scale_fill_nejm()+
  stat_compare_means(method='kruskal.test')+labs(x='Score', title = 'ImmuneScore')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        legend.position = '')
b2

b1+b2
ggsave('05-doEstimate/out-plot/estimate_boxplot.pdf', height = 5,width = 8)

b3 = ggplot(dat = data.score,aes(group,ESTIMATEScore ))+
  geom_boxplot(aes(fill = group))+
  scale_fill_nejm()+
  stat_compare_means(method='t.test')+labs(x='Score', title = 'ESTIMATEScore')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        legend.position = '')
b3
b1+b2+b3
ggsave('05-doEstimate/out-plot/estimate_boxplot_3.pdf', height = 5,width = 8)





