# Settings
# Objective: To prepare data for immunoassay: RNA expression matrix of basel tumor sample is required
# Data Source:
# Save data:

# Create a folder in the current directory
getwd()
dir.create('./05-doEstimate//out-plot')
dir.create('./05-doEstimate//out-data')

rm(list=ls())
options(stringsAsFactors = F)
library(stringr) 

# 0. load data----
library(data.table)
rna_exp = fread( 'data_input/FUSCC//tnbc448.counts.GRCh38.99.txt',header = T, data.table = F)
row.names(rna_exp) = rna_exp$Geneid
rna_exp[1:4,1:4]

rna_exp<-rna_exp[,-c(1:6)] 
colnames(rna_exp)<-gsub("^...Lib_","",gsub(".bam","",colnames(rna_exp))) 
keep_feature <- rowSums(rna_exp>1 ) > 10 
rna_exp<-rna_exp[keep_feature,]
rna_exp[1:4,1:4]
dim(rna_exp)

#PT','tumor'
colnames(rna_exp)
gp=ifelse(grepl('PT',colnames(rna_exp)),'PT','tumor')
table(gp)

#2、subtyp
a=read.table('data_input//subtype.txt',header = T)
head(a) 

a$Project_ID   
pid=gsub('.rep','',gsub('.PT','',colnames(rna_exp)))
pos=match(pid,a$Project_ID  ) 
subtype=a[pos,2]
table(subtype)
d=data.frame(gp=gp,
             pid=pid,
             subtype=subtype)
d.tumor = d[d$gp == 'tumor',]
d.tumor.BLIS = d.tumor[d.tumor$subtype == 'BLIS',]
dim(d.tumor.BLIS)

colnames(rna_exp) = gsub('.rep','',colnames(rna_exp))  
exp.basel = rna_exp[,d.tumor.BLIS$pid]
dim(exp.basel)
exp.basel[1:4,1:4]
library(clusterProfiler)
library(org.Hs.eg.db)
genenames<-bitr(rownames(rna_exp),fromType = "ENSEMBL",
                toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
head(genenames)
length(unique(genenames$ENSEMBL))
length(unique(genenames$SYMBOL))

genenames=genenames[!duplicated(genenames$SYMBOL),]
exp.basel<-exp.basel[match(genenames$ENSEMBL,rownames(exp.basel)),]
rownames(exp.basel)=genenames$SYMBOL
exp.basel[1:4,1:4]
rownames(exp.basel)[1:50]
dim(exp.basel)
save(exp.basel,file = '05-doEstimate/out-data/nolog_exp.basel.Rdata')



