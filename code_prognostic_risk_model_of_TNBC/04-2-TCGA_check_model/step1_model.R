
# rm(list = ls())
# options(stringsAsFactors = F)
# 
# library(data.table)
# library(dplyr)
# library(survminer) 
# library(survival)
# library(loose.rock)
# library(futile.logger) 
# library(glmSparseNet)
# library(ggrisk)  
# library(pheatmap) 
# library(ggplot2)
# library(ggsci)
# library(ggstatsplot)
# 
# 
# load('04-3-TCGA_check_model//out_data/input.Rdata')
# gene_overlap <- read.table('03-doOverlap-DEG-COX/out-data//03_overlap_logFC_t_1.5_DEG_with_Cox_overlap_genes.txt')[, 1]
# gene_overlap
# 
# # 1. tidy data for LASSO----
# ## prepare xdata----
# kp <- gene_overlap
# cdat <- exp.tcga.basal[kp, ]
# cdat <- na.omit(cdat)
# xdata <- t(cdat)
# xdata[1:4,1:4]
# 
# xdata= as.data.frame(xdata)
# colnames(xdata) = gsub('-','_',colnames(xdata))  
# colnames(xdata)
# 
# ## prepare ydata----
# ydata <- phe.basal
# head(ydata)
# boxplot(ydata$time)
# kp <- ydata$time >0
# table(kp)
# xdata <- xdata[kp,]
# ydata <- ydata[kp,]
# identical(rownames(xdata), rownames(ydata))
# 
# data = cbind(ydata,xdata)
# 
# 
# save(xdata,ydata,data,file = '04-3-TCGA_check_model/out_data//input_stepwise.Rdata')
# 
