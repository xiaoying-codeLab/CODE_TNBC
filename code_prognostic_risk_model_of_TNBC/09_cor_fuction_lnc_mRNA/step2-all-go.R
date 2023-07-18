
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

data.cor = read.csv("11_cor_fuction_lnc_mRNA/cor.csv")
gene_model <- read.table('04-1-doModel-lasso/out-data//model_lasso_cg_genes.txt')[, 1]
gene_model   

########################################################################
#pos。
########################################################################
data.cor.up = data.cor[data.cor$gene %in% gene_model,]
table(data.cor.up$gene)
data.cor.up.pos = data.cor.up[data.cor.up$cor > 0.5,]
data.cor.up.pos = data.cor.up.pos[data.cor.up.pos$p.val <0.01,]
library(clusterProfiler)
library(org.Hs.eg.db)
genenames<-bitr(data.cor.up.pos$Sample,fromType = "SYMBOL",
                toType = "ENTREZID",OrgDb =" org.Hs.eg.db")
go <- enrichGO(genenames$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all") 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave('11_cor_fuction_lnc_mRNA/last_plot/go_up_nes.pdf') 


################################################
#   1、下调基因----pos。
########################################################################
rm(list=ls())
options(stringsAsFactors = F)
#
data.cor = read.csv("11_cor_fuction_lnc_mRNA/cor.csv")
gene_model <- read.table('04-1-doModel-lasso/out-data//model_lasso_cg_genes.txt')[, 1]

data.cor.up = data.cor[data.cor$gene %in% gene_model,]
table(data.cor.up$gene)
data.cor.up.pos = data.cor.up[data.cor.up$cor< -0.4,]
data.cor.up.pos = data.cor.up.pos[data.cor.up.pos$p.val <0.01,]
#
library(clusterProfiler)
library(org.Hs.eg.db)
genenames<-bitr(data.cor.up.pos$Sample,fromType = "SYMBOL",
                toType = "ENTREZID",OrgDb =" org.Hs.eg.db")

# 
go <- enrichGO(genenames$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all",
               pvalueCutoff = 0.01,
               qvalueCutoff = 0.2, 
               minGSSize = 100, maxGSSize = 500) 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave('11_cor_fuction_lnc_mRNA/last_plot/go_nes.pdf') 


