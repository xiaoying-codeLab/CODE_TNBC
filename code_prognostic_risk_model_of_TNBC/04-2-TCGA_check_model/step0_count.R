# Settings
# Purpose: To process TCGA data first. Convert id to lncRNA- all Ensembl does the conversion centrally
# Data Source:
# Save data:
rm(list = ls())
options(stringsAsFactors = F)

library(stringi)
library(data.table)


##############################################################################################
#tidy
exp.tcga=fread(file = 'data_input/TCGA//TCGA-BRCA.htseq_counts.tsv.gz',data.table = F)
dim(exp.tcga)  #
exp.tcga[1:4,1:4]
exp.tcga$Ensembl_ID=stri_sub(exp.tcga$Ensembl_ID,1,15)
exp.tcga=exp.tcga[,-1]
exp.tcga[1:4,1:4]

subtype =fread('data_input/TCGA//TCGA.BRCA.sampleMap_BRCA_clinicalMatrix',data.table = F)
colnames(subtype)
subtype = subtype[subtype$ER_Status_nature2012 == "Negative",]
subtype = subtype[subtype$HER2_Final_Status_nature2012 == "Negative",]
subtype = subtype[subtype$PR_Status_nature2012 == "Negative",]
dim(subtype)

colnames(exp.tcga) = stri_sub(colnames(exp.tcga),1,15)   
tumor = subtype$sampleID[subtype$sampleID %in% colnames(exp.tcga)]
tumor
exp.tcga.basal = exp.tcga[,c(tumor)]
exp.tcga.basal[1:4,1:4]
dim(exp.tcga.basal)

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
dim(exp.tcga.basal)


# Survival data
# Survival analysis data
library(data.table)
phe=fread('data_input/TCGA//TCGA.BRCA.sampleMap_BRCA_clinicalMatrix',data.table = F)
phe[1:4,1:4]
colnames(phe)
phe =phe[,c('sampleID',"OS_Time_nature2012","OS_event_nature2012")]
phe=na.omit(phe) 
dim(phe)    
head(phe)  
colnames(phe) <- c('sampleID','time','status')
row.names(phe) <-  phe$sampleID
phe=phe[phe$time>0,]
head(phe)
dim(phe)  

phe.basal = phe[subtype$sampleID,]
phe.basal = na.omit(phe.basal) 
phe.basal$time <- phe.basal$time/365 
boxplot(phe.basal$time)
head(phe.basal)

# Match the expression matrix with the survival information
exp.tcga.basal = exp.tcga.basal[,colnames(exp.tcga.basal) %in% rownames(phe.basal)]
phe.basal = phe.basal[match(colnames(exp.tcga.basal),rownames(phe.basal)),]
identical(colnames(exp.tcga.basal),rownames(phe.basal))


# Save data
dim(exp.tcga.basal)
dim(phe.basal)
save(exp.tcga.basal,phe.basal,file = '04-3-TCGA_check_model//out_data/input.Rdata')




