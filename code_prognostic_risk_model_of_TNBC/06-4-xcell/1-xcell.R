# Do xcell rating
rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(ggsci)
library(ggpubr)
library(patchwork)

load('05-doEstimate/out-data/exp.basel.Rdata')
library(xCell)
scores = xCellAnalysis(exp.basel, rnaseq=F)
save(scores,file = '06-xcell/score.Rdata')
load(file = '06-xcell/score.Rdata')
load('04-1-doModel-lasso/out-data//basel-risk.Rdata')
rt =new_dat
dim(rt)
rt[1:4, 1:4]
table(rt$Risk_level)
scores = as.data.frame(t(scores))

identical(rownames(rt),rownames(scores))
scores$group = rt$Risk_level
#########################################################
library(dplyr)
library(rstatix)
library(reshape2)
colnames(scores)
b1 = ggplot(dat = scores, aes(group,StromaScore))+
  geom_boxplot(aes(fill = group))+
  scale_fill_nejm()+
  stat_compare_means(method='t.test') + labs(x='Score', title = 'StromaScore')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        legend.position = '')
b1
b2 = ggplot(dat = scores,aes(group,ImmuneScore ))+
  geom_boxplot(aes(fill = group))+
  scale_fill_nejm()+
  stat_compare_means(method='kruskal.test')+labs(x='Score', title = 'ImmuneScore')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        legend.position = '')
b2
b3 = ggplot(dat = scores,aes(group,MicroenvironmentScore ))+
  geom_boxplot(aes(fill = group))+
  scale_fill_nejm()+
  stat_compare_means(method='t.test')+labs(x='Score', title = 'MicroenvironmentScore')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5),
        axis.title = element_text(face = 'bold'),
        legend.position = '')
b3
b1+b2
ggsave('06-xcell//estimate_boxplot.pdf', height = 5,width = 8)

