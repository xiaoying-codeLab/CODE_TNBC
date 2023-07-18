# ceRNA Network building

# 1. Lncrnas with differences between high-risk and low-risk groups were screened
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

# Read data
# lnc represents the matrix
# Settings
# Purpose: Preprocess data. Find the matrix of bssel and normal
# Data source: input
# Save data:

# Create a folder in the current directory
getwd()
dir.create('./00-folder/out-plot')
dir.create('./00-folder/out-data')

rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(stringi)
library(data.table)

library(data.table)
a = fread( '../data_input/FUSCC//tnbc448.counts.lncipedia_5_2_hc_hg38.txt',
           header = T, data.table = F)
dim(a)  #
a[1:4,1:4]
expr= a[,7:ncol(a)]   #
rownames(expr)=a$Geneid
expr[1:4,1:4]
colnames(expr)
colnames(expr)=gsub('^...Lib_','',
                    gsub('.bam','',  colnames(expr)))
colnames(expr)
colnames(expr) = gsub('.rep','',colnames(expr))  
dim(expr)
expr[1:4,1:4]


lnc.anno = read.table(file = '../data_input/LNCipedia/lncipedia_5_2_ensembl_92_genes.txt',header = T)
dim(lnc.anno)
gene
write.table(gene,file = "F4_ceRNA//lncRNA.txt",
            row.names = FALSE,col.names = FALSE,quote=FALSE)




table(rownames(expr) %in% lnc.anno$lncipediaGeneID)

expr = expr[rownames(expr) %in% lnc.anno$lncipediaGeneID,]
dim(expr)

lnc.anno=lnc.anno[!duplicated(lnc.anno$lncipediaGeneID),]  
rownames(lnc.anno) =lnc.anno$lncipediaGeneID
lnc.anno = lnc.anno[rownames(expr),]
identical(rownames(expr),rownames(lnc.anno))
rownames(expr) = stri_sub(lnc.anno$ensemblGeneID,1,15)#
expr[1:4,1:4]
expr=log2(edgeR::cpm(expr)+1)#
boxplot(expr[,1:4],las=2) 
keep_feature <- rowSums (expr > 1) > 10
table(keep_feature)  #
expr <- expr[keep_feature, ]
expr[1:4,1:4] 
dim(expr)
expr_lnc = expr


# 
load('../04-1-doModel-lasso/out-data//basel-risk.Rdata')
group = new_dat
expr_lnc = expr_lnc[,rownames(new_dat)]
group = group[colnames(expr_lnc),]
identical(rownames(group),colnames(expr_lnc))
group = group$Risk_level
table(group)
group <- factor(group, levels = c('High','Low'))
table(group)

# 
dim(expr_lnc)
expr_lnc[1:4,1:4]
table(group)  #
levels(group)
expr_lnc = expr_lnc[]


library(limma)
design=model.matrix(~factor( group ))
fit=lmFit(expr_lnc,design)
fit=eBayes(fit)
## 
options(digits = 4) 
deg = topTable(fit,coef=2,adjust='BH', n=Inf) 
deg[rownames(deg)%in% c("ENSG00000249001"),]




colnames()
save(deg, file='01-doGEG-normal-basel//out-data//1.limma_deg.Rdata') 




####################################################################################
# 2、
library(data.table)
exp_mRNA = fread( '../data_input/FUSCC//tnbc448.counts.GRCh38.99.txt',header = T, data.table = F)
row.names(exp_mRNA) = exp_mRNA$Geneid
exp_mRNA[1:4,1:4]

exp_mRNA<-exp_mRNA[,-c(1:6)]  
colnames(exp_mRNA)<-gsub("^...Lib_","",gsub(".bam","",colnames(exp_mRNA)))  #
exp_mRNA=log2(edgeR::cpm(exp_mRNA)+1)#
keep_feature <- rowSums(exp_mRNA>1 ) > 10 #
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


load('../04-1-doModel-lasso/out-data//basel-risk.Rdata')
group = new_dat
exp_mRNA = as.data.frame(exp_mRNA)
colnames(exp_mRNA)
colnames(exp_mRNA)<-gsub(".rep","",colnames(exp_mRNA))  #
exp_mRNA = exp_mRNA[,rownames(group)]
identical(rownames(group),colnames(exp_mRNA))
group = group$Risk_level
table(group)
group <- factor(group, levels = c('High','Low'))
table(group)


dim(exp_mRNA)
exp_mRNA[1:4,1:4]
table(group)  #
levels(group)


library(limma)
design=model.matrix(~factor( group ))
fit=lmFit(exp_mRNA,design)
fit=eBayes(fit)

deg_mrna = topTable(fit,coef=2,adjust='BH', n=Inf) 

deg_mrna[rownames(deg_mrna) %in% c("KRT5","NPR3","MFAP5","DSG1","CCND1"),]


deg_mrna$g=ifelse(deg_mrna$P.Value>0.05,'stable', 
            ifelse( deg_mrna$logFC >1,'up', 
                    ifelse( deg_mrna$logFC < -1,'down','stable') )
)
table(deg_mrna$g)
save(deg_mrna, file='F4_ceRNA//limma_deg_mrna.Rdata') 

gene_up = deg_mrna[deg_mrna$g == "up",]
gene_down = deg_mrna[deg_mrna$g == "down",]
gene_diff = c(rownames(gene_up),rownames(gene_down))
write.table(gene_diff,file = "../F4_ceRNA//diff_mrna3.txt",
            row.names = FALSE,col.names = FALSE,quote=FALSE)





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
gene_model   #

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



library(Hmisc)#加载包
# 
load('11_cor_fuction_lnc_mRNA/cor_result.Rdata')
data = apply(cor_deg_mrna, 2,function(x)
  x>0.5)
table(data)
data = as.data.frame(data)
data = ifelse(data == "FALSE",0,1)
data = as.data.frame(data)
data$sum = apply(data, 1,sum)
data = data[data$sum > 0,]
dim(data)
gene_pos = rownames(data)


data = apply(cor_deg_mrna, 2,function(x)
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

#   
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

# 
identical(data.cor.R$Sample,data.cor.P$Sample)
data.cor.R$p.val = data.cor.P$P.val

write.csv(data.cor.R,file = "11_cor_fuction_lnc_mRNA/cor.csv")

