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
# # # # # # # # # # # # # # # # # # # # #
#1. Process data
#bssel and normal matrix
# # # # # # # # # # # # # # # # # # # # #
library(data.table)
a = fread( 'data_input/FUSCC//tnbc448.counts.lncipedia_5_2_hc_hg38.txt',
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
colnames(expr) = gsub('.rep','',colnames(expr))  #去掉rep才好匹配.
dim(expr)
expr[1:4,1:4]

# Match lnc
#matches the transformation table provided by the LNCipedia database to convert to Ensembl
lnc.anno = read.table(file = 'data_input/LNCipedia/lncipedia_5_2_ensembl_92_genes.txt',header = T)
dim(lnc.anno)
table(rownames(expr) %in% lnc.anno$lncipediaGeneID)
# Filter the lnc that matches
expr = expr[rownames(expr) %in% lnc.anno$lncipediaGeneID,]
dim(expr)

expr=log2(edgeR::cpm(expr)+1)
boxplot(expr[,1:4],las=2) 
keep_feature <- rowSums (expr > 1) > 10
table(keep_feature)  
expr <- expr[keep_feature, ]
expr[1:4,1:4] 
dim(expr)


#################
#  Read the TCGA data. And we need to match each other.
################
##############################################################################################
#Read the expression matrix representing TCGA
#tidy
exp.tcga=fread(file = 'data_input/TCGA/TCGA-BRCA.htseq_counts.tsv.gz',data.table = F)
dim(exp.tcga)  #
exp.tcga[1:4,1:4]
exp.tcga$Ensembl_ID=stri_sub(exp.tcga$Ensembl_ID,1,15)
rownames(exp.tcga)=exp.tcga[,1]
exp.tcga=exp.tcga[,-1]
exp.tcga[1:4,1:4]


# Grow a Basal tumor
# To phenotype data:
# Filter directly from clinicalMatrix
subtype =fread('data_input/TCGA/TCGA.BRCA.sampleMap_BRCA_clinicalMatrix',data.table = F)
Basal = subtype[subtype$PAM50Call_RNAseq == 'Basal',]
dim(Basal)  
colnames(exp.tcga) =stri_sub(colnames(exp.tcga),1,15)
exp.tcga.basal = exp.tcga[,colnames(exp.tcga) %in% Basal$sampleID]
dim(exp.tcga.basal)

# Then convert the id to lnc
# Match lnc
#1, matches the transformation table provided by the LNCipedia database to convert to Ensembl
lnc.anno = read.table(file = 'data_input/LNCipedia//lncipedia_5_2_ensembl_92_genes.txt',header = T)
dim(lnc.anno)
lnc.anno=lnc.anno[!duplicated(lnc.anno$lncipediaGeneID),]   
lnc.anno$ensemblGeneID = stri_sub(lnc.anno$ensemblGeneID,1,15)
table(rownames(exp.tcga.basal) %in% lnc.anno$ensemblGeneID)


exp.tcga.basal = exp.tcga.basal[rownames(exp.tcga.basal) %in% lnc.anno$ensemblGeneID,]
lnc.anno = lnc.anno[match(rownames(exp.tcga.basal),lnc.anno$ensemblGeneID),]
identical(rownames(exp.tcga.basal),lnc.anno$ensemblGeneID)
rownames(exp.tcga.basal) = lnc.anno$lncipediaGeneID
exp.tcga.basal[1:4,1:4]


boxplot(exp.tcga.basal[,1:4],las=2) 
keep_feature <- rowSums (exp.tcga.basal > 1) > 10
table(keep_feature) 
exp.tcga.basal <- exp.tcga.basal[keep_feature, ]
exp.tcga.basal[1:4,1:4] 
dim(exp.tcga.basal)


############
#  Then take the expression matrix, which needs to be in TCGA as well
###########
table(rownames(exp.tcga.basal) %in%  row.names(expr))
expr = expr[row.names(expr) %in% rownames(exp.tcga.basal),]
dim(expr)


#####################
#2、subtyp-
#####################
library(xlsx)
pd = read.xlsx("data_input/FUSCC/mmc2.xlsx",sheetIndex = 1,startRow = 2)
colnames(pd) 
subtype = pd[,c("Project_ID","mRNA_Subtype")]
head(subtype) 
dim(subtype)

expr_normal = expr[,grep('.PT',colnames(expr))]
colnames(expr_normal)
dim(expr_normal)

BLIS_pd = subtype[subtype$mRNA_Subtype == 'BLIS',]
BLIS_pd =na.omit(BLIS_pd)
dim(BLIS_pd)
head(BLIS_pd)

#3. Select the matrix with only BLIS and normal
table(colnames(expr) %in% BLIS_pd$Project_ID)
expr_BLIS = expr[,BLIS_pd$Project_ID]
dim(expr_BLIS)
dim(expr_normal)
expr_lnc = cbind(expr_normal,expr_BLIS) 
dim(expr_lnc)

#4. Group information
library(stringr)
gp = ifelse(grepl('PT', colnames(expr_lnc)),'Normal','BLIS') 
table(gp)
gp <- factor(gp, levels = c('Normal','BLIS'))
table(gp)
gp 

# Save data: Matrix + grouping information
expr_lnc[1:4,1:4]
table(gp)
save(expr_lnc,gp,file='00-folder/out-data/0.expr.Rdata')


###########
#Survival data
#########
library(data.table)
phe = pd[,c("Project_ID","RFS_Status","RFS_time_Days","RFS_time_Months")]
# phe <- fread('data_input//clinical.txt',header = T, data.table = F) 
head(phe)
rownames(phe) = phe$Project_ID
phe <- phe[, c(4,2)]
colnames(phe) <- c('time', 'event')
head(phe)
boxplot(phe$time)  #月的时间
phe$time <- phe$time/12
head(phe)
boxplot(phe$time) 
kp <- phe$time > 0
table(kp)
phe <- phe[kp,]



# Match the sample of survival data with the sample of representation matrix
## Just sample BLIS
table(colnames(expr_BLIS) %in% row.names(phe))
phe_BLIS = phe[match(colnames(expr_BLIS),row.names(phe)),]
dim(phe_BLIS)
head(phe_BLIS)

save(phe_BLIS,file = '00-folder/out-data/0.survival.Rdata')







