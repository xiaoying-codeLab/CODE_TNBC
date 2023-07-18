# Visualization, kegg

# Settings
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(ggsci)
library(ggpubr)
library(patchwork)


#读入数据
data.cor = read.csv("11_cor_fuction_lnc_mRNA/cor.csv")
# Do it separately next.
gene_model <- read.table('04-1-doModel-lasso/out-data//model_lasso_cg_genes.txt')[, 1]
gene_model   # 
gene_up = c("C5orf66-AS2","LINC00393","lnc-ERI1-32","lnc-SPARCL1-1")
gene_down = c("DIO3OS","FZD10-DT","HCG23","lnc-FOXO1-2","lnc-MMD-4","lnc-TMEM106C-6")
