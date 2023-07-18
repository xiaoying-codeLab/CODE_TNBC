
# Purpose: To correlate lnc-mRNA and find the relevant mRNA- then go and kegg
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


load('00-folder/out-data/0.expr.Rdata')
dim(expr_lnc)
expr_lnc[1:4,1:4]
table(gp)  
expr_lnc = expr_lnc[,gp == "BLIS"]
dim(expr_lnc)
colnames(expr_lnc)
gene_model <- read.table('04-1-doModel-lasso/out-data//model_lasso_cg_genes.txt')[, 1]
gene_model   
expr_lnc = expr_lnc[gene_model,]
library(data.table)
exp_mRNA = fread( 'data_input//tnbc448.counts.GRCh38.99.txt',header = T, data.table = F)
row.names(exp_mRNA) = exp_mRNA$Geneid
exp_mRNA[1:4,1:4]
exp_mRNA<-exp_mRNA[,-c(1:6)]  
colnames(exp_mRNA)<-gsub("^...Lib_","",gsub(".bam","",colnames(exp_mRNA))) 
exp_mRNA=log2(edgeR::cpm(exp_mRNA)+1)
keep_feature <- rowSums(exp_mRNA>1 ) > 10 
exp_mRNA<-exp_mRNA[keep_feature,]
exp_mRNA[1:4,1:4]
dim(exp_mRNA)
library(clusterProfiler)
library(org.Hs.eg.db)
genenames<-bitr(rownames(exp_mRNA),fromType = "ENSEMBL",
                toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
head(genenames)
length(unique(genenames$ENSEMBL))
length(unique(genenames$SYMBOL))
genenames=genenames[!duplicated(genenames$SYMBOL),]
exp_mRNA<-exp_mRNA[match(genenames$ENSEMBL,rownames(exp_mRNA)),]
rownames(exp_mRNA)=genenames$SYMBOL
exp_mRNA[1:4,1:4]
# Find correlation
## Merge the two data. Two sets of data. mRNA+lncRNA
dim(exp_mRNA)
dim(expr_lnc)
table(colnames(expr_lnc) %in% colnames(exp_mRNA))
colnames(exp_mRNA)
colnames(exp_mRNA)<-gsub("*.rep","",colnames(exp_mRNA)) 
exp_mRNA = exp_mRNA[,colnames(expr_lnc)]
identical(colnames(exp_mRNA),colnames(expr_lnc))
expr_all=rbind(exp_mRNA,expr_lnc)
dim(expr_all)
gene_model <- read.table('04-1-doModel-lasso/out-data//model_lasso_cg_genes.txt')[, 1]
gene_model  
data.cor = read.csv("11_cor_fuction_lnc_mRNA/cor.csv")
data.cor.up = data.cor[data.cor$gene %in% gene_model,]
table(data.cor.up$gene)


data.cor.pos = data.cor.up[data.cor.up$cor > 0.5,]
data.cor.pos = data.cor.pos[data.cor.pos$p.val <0.01,]
data.cor.pos = na.omit(data.cor.pos)

data.cor.nes = data.cor.up[data.cor.up$cor < -0.4,]
data.cor.nes = data.cor.nes[data.cor.nes$p.val <0.01,]


gene_all = c(data.cor.pos$Sample,data.cor.pos$gene,data.cor.nes$Sample,data.cor.nes$gene)
gene_all = gene_all[!duplicated(gene_all)]
gene_all

expr_all2 = expr_all[gene_all,]
expr_all2 = as.data.frame(expr_all2)
save(expr_all2,file = "11_cor_fuction_lnc_mRNA/expr_all2.Rdata")

library(Hmisc)
load('11_cor_fuction_lnc_mRNA/cor_result.Rdata')
data = apply(cor_df, 2,function(x)
  x>0.5)
table(data)
data = as.data.frame(data)
data = ifelse(data == "FALSE",0,1)
data = as.data.frame(data)
data$sum = apply(data, 1,sum)
data = data[data$sum > 0,]
dim(data)
gene_pos = rownames(data)


data = apply(cor_df, 2,function(x)
  x < -0.4)
table(data)
data = as.data.frame(data)
data = ifelse(data == "FALSE",0,1)
data = as.data.frame(data)
data$sum = apply(data, 1,sum)
data = data[data$sum > 0,]
dim(data)
gene_nes = rownames(data)
gene_cor = unique(c(gene_pos,gene_nes,gene_model)) 

#   然后用这个包来计算。
expr_all[1:4,1:4]
expr_all = expr_all[gene_cor,]
data.cor <- rcorr(as.matrix(t(expr_all)))
data.cor.P = as.data.frame(data.cor$P)
data.cor.R = as.data.frame(data.cor$r)
data.cor.R <- data.cor.R %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = gene,value = cor,-Sample)

data.cor.P <- data.cor.P %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = gene,value = P.val,-Sample)
identical(data.cor.R$Sample,data.cor.P$Sample)
data.cor.R$p.val = data.cor.P$P.val

write.csv(data.cor.R,file = "11_cor_fuction_lnc_mRNA/cor.csv")

