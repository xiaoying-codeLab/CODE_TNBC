# Settings
# Purpose: Prepare xdata and ydata - then create lasso model
# Data Source:
# Save data:

# Create a folder in the current directory
# Create a folder in the current directory
getwd()
dir.create('04-1-doModel-lasso/out-plot')
dir.create('04-1-doModel-lasso//out-data')

rm(list=ls())
options(stringsAsFactors = F)


library(data.table)
library(dplyr)
library(survminer) 
library(survival)
library(loose.rock)
library(futile.logger) 
library(glmSparseNet)
library(ggrisk)  
library(pheatmap) 
library(ggplot2)
library(ggsci)
library(ggstatsplot)


# 0. Survival data
load(  file = '02-doCOX/out-data/2.input_survival_data.Rdata')
expr = expr.BLIS
phe = phe_BLIS
expr[1:4,1:4] 
head(phe)
identical(colnames(expr), rownames(phe))
fs = list.files(path = '03-doOverlap-DEG-COX/out-data/',pattern = 'Cox_overlap_genes.txt')
fs
gene_overlap <- read.table('03-doOverlap-DEG-COX/out-data//03_overlap_logFC_t_2_DEG_with_Cox_overlap_genes.txt')[, 1]
gene_overlap

# 1. tidy data for LASSO----
## prepare xdata----
kp <- gene_overlap
cdat <- expr[kp, ]
xdata <- t(cdat)
xdata[1:4,1:4]

## prepare ydata----
ydata <- phe
head(ydata)
identical(rownames(xdata), rownames(ydata))
## time transform to year
head(ydata)
boxplot(ydata$time)
kp <- ydata$time >0
table(kp)
xdata <- xdata[kp,]
ydata <- ydata[kp,]


# 2. run lasso model----
source('04-1-doModel-lasso//run_LASSO_cox.R') 

lasso_cox(xdata,
          ydata)


save(xdata,ydata,file = '04-1-doModel-lasso//input_model_data.Rdata')


gene_model <- read.table('04-1-doModel-lasso/out-data//model_lasso_cg_genes.txt')[, 1]
gene_model   






