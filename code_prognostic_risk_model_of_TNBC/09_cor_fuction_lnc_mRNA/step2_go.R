
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
#  
gene_model <- read.table('04-1-doModel-lasso/out-data//model_lasso_cg_genes.txt')[, 1]
gene_model   # 
gene_up = c("C5orf66-AS2","LINC00393","lnc-ERI1-32","lnc-SPARCL1-1")
gene_down = c("DIO3OS","FZD10-DT","HCG23","lnc-FOXO1-2","lnc-MMD-4","lnc-TMEM106C-6")

########################################################################
#   1. Up-regulated gene ----pos.
########################################################################
data.cor.up = data.cor[data.cor$gene %in% gene_up,]
table(data.cor.up$gene)
data.cor.up.pos = data.cor.up[data.cor.up$cor > 0.5,]
data.cor.up.pos = data.cor.up.pos[data.cor.up.pos$p.val <0.01,]

library(clusterProfiler)
library(org.Hs.eg.db)
genenames<-bitr(data.cor.up.pos$Sample,fromType = "SYMBOL",
                toType = "ENTREZID",OrgDb =" org.Hs.eg.db")

# go
go <- enrichGO(genenames$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all") 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave('11_cor_fuction_lnc_mRNA/go_up_pos.pdf') 



########################################################################
#  1. 上调基因----pos。
########################################################################
data.cor.up = data.cor[data.cor$gene %in% gene_up,]
table(data.cor.up$gene)
data.cor.up.pos = data.cor.up[data.cor.up$cor < -0.4,]
data.cor.up.pos = data.cor.up.pos[data.cor.up.pos$p.val <0.01,]
#id
library(clusterProfiler)
library(org.Hs.eg.db)
genenames<-bitr(data.cor.up.pos$Sample,fromType = "SYMBOL",
                toType = "ENTREZID",OrgDb =" org.Hs.eg.db")

# g0
go <- enrichGO(genenames$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all") 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave('11_cor_fuction_lnc_mRNA/go_up_nes.pdf') 









################################################
#   1、下调基因----pos。
########################################################################
rm(list=ls())
options(stringsAsFactors = F)
#读入数据
data.cor = read.csv("11_cor_fuction_lnc_mRNA/cor.csv")
gene_model <- read.table('04-1-doModel-lasso/out-data//model_lasso_cg_genes.txt')[, 1]
gene_up = c("C5orf66-AS2","LINC00393","lnc-ERI1-32","lnc-SPARCL1-1")
gene_down = c("DIO3OS","FZD10-DT","HCG23","lnc-FOXO1-2","lnc-MMD-4","lnc-TMEM106C-6")

data.cor.up = data.cor[data.cor$gene %in% gene_down,]
table(data.cor.up$gene)
data.cor.up.pos = data.cor.up[data.cor.up$cor > 0.5,]
data.cor.up.pos = data.cor.up.pos[data.cor.up.pos$p.val <0.01,]
#id转换
library(clusterProfiler)
library(org.Hs.eg.db)
genenames<-bitr(data.cor.up.pos$Sample,fromType = "SYMBOL",
                toType = "ENTREZID",OrgDb =" org.Hs.eg.db"
                )

# go分析
go <- enrichGO(genenames$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all",
               pAdjustMethod = "none"
               ) 
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave('11_cor_fuction_lnc_mRNA/go_down_pos.pdf') 


########################################################################
#   1、下调基因----nes。
########################################################################
rm(list=ls())
options(stringsAsFactors = F)
#读入数据
data.cor = read.csv("11_cor_fuction_lnc_mRNA/cor.csv")
gene_model <- read.table('04-1-doModel-lasso/out-data//model_lasso_cg_genes.txt')[, 1]
gene_up = c("C5orf66-AS2","LINC00393","lnc-ERI1-32","lnc-SPARCL1-1")
gene_down = c("DIO3OS","FZD10-DT","HCG23","lnc-FOXO1-2","lnc-MMD-4","lnc-TMEM106C-6")


data.cor.up = data.cor[data.cor$gene %in% gene_down,]
table(data.cor.up$gene)
data.cor.up.pos = data.cor.up[data.cor.up$cor < -0.4,]
data.cor.up.pos = data.cor.up.pos[data.cor.up.pos$p.val <0.01,]
#id转换
library(clusterProfiler)
library(org.Hs.eg.db)
genenames<-bitr(data.cor.up.pos$Sample,fromType = "SYMBOL",
                toType = "ENTREZID",OrgDb =" org.Hs.eg.db")

# go分析
go <- enrichGO(genenames$ENTREZID, OrgDb = "org.Hs.eg.db", ont="all") 

pdf(file='11_cor_fuction_lnc_mRNA/go_down_nes.pdf',width = 10,height = 8)
barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
dev.off()

