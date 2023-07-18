# Settings
Objective: Whether there are differences in the expression levels of these genes in TCGA. # Let's make a bar chart first:
# Data Source:
# Save data:
# Let's make a bar chart first:

rm(list = ls())
options(stringsAsFactors = F)

library(stringi)
library(data.table)

#loda data
cg_bedtool <- read.csv('07_TCGA_check_model//input_data/cg_bedtool_change.csv',header = T,fill = T)
colnames(cg_bedtool) = c('lnc_id','Ensembl_id')
cg_bedtool$Ensembl_id = stri_sub(cg_bedtool$Ensembl_id,1,15) 
head(cg_bedtool)
##############################################################################################
#tidy
exp.tcga=fread(file = '07_TCGA_check_model/input_data/TCGA-BRCA.htseq_counts.tsv.gz',data.table = F)
dim(exp.tcga)  #
exp.tcga[1:4,1:4]
exp.tcga$Ensembl_ID=stri_sub(exp.tcga$Ensembl_ID,1,15)
rownames(exp.tcga)=exp.tcga[,1]
exp.tcga=exp.tcga[,-1]
exp.tcga[1:4,1:4]
table(cg_bedtool$Ensembl_id%in% row.names(exp.tcga))   
cg_bedtool[cg_bedtool$Ensembl_id%in% row.names(exp.tcga),]
cg_bedtool[!(cg_bedtool$Ensembl_id%in% row.names(exp.tcga)),]

tcga.cg = exp.tcga[cg_bedtool$Ensembl_id,]
tcga.cg = na.omit(tcga.cg )
dim(tcga.cg)
tcga.cg[1:4,1:4]

colnames(tcga.cg)
colnames(tcga.cg) = stri_sub(colnames(tcga.cg),1,15)
subtype =fread('07_TCGA_check_model/input_data/TCGA.BRCA.sampleMap_BRCA_clinicalMatrix',data.table = F)
Basal = subtype[subtype$PAM50Call_RNAseq == 'Basal',]
dim(Basal)  
normal = subtype[stri_sub(subtype$sampleID,14) == '11',]
dim(normal)
sampel = c(Basal$sampleID,normal$sampleID)
sampel[1:4]
expr.tcga = tcga.cg[,colnames(tcga.cg) %in% sampel]
dim(expr.tcga)
group_list=ifelse(colnames(expr.tcga) %in% Basal$sampleID,'Basal','normal')
table(group_list)

#Draw a bar chart of total expressions
library(ggplot2)
library(ggstatsplot) 
dat =t(expr.tcga)
dat = as.data.frame(dat)
dat[1:4,1:4]
dat$group = group_list
library(reshape2) 
dat2 <- melt(dat)
colnames(dat2)
##########################################
library(ggpubr)
p <- ggboxplot(dat2, x="variable", y="value", color = "group", 
               palette = "jco", add = "jitter")#添加p-valuep+stat_compare_means()
p


ggboxplot(dat2, x="variable", y="value", color = "group", palette = "jco")+
  stat_compare_means()


p <- ggboxplot(dat2, x="variable", y="value", color = "group", 
               palette = "jco", add = "jitter")
p
p+stat_compare_means(aes(group=group))
p+stat_compare_means(aes(group=group), label = "p.format")
p+stat_compare_means(aes(group=group), label = "p.signif")


setwd('07_TCGA_check_model/out_plot/')
ggsave('expr.DEG.pdf')
setwd('../../')
getwd()
