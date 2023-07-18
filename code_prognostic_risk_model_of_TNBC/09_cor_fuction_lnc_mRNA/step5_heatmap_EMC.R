# gese
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


load("04-1-doModel-lasso/out-data/basel-risk.Rdata")
library(data.table)
exp_mRNA = fread( 'data_input/FUSCC//tnbc448.counts.GRCh38.99.txt',header = T, data.table = F)
row.names(exp_mRNA) = exp_mRNA$Geneid
exp_mRNA[1:4,1:4]

exp_mRNA<-exp_mRNA[,-c(1:6)]  
colnames(exp_mRNA)<-gsub("^...Lib_","",gsub(".bam","",colnames(exp_mRNA)))
exp_mRNA=log2(edgeR::cpm(exp_mRNA)+1)
keep_feature <- rowSums(exp_mRNA>1 ) > 10 
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


# Read in strom-related genes
gene_emc = read.table("data_input/GSEA/heatmap//KEGG_ECM_RECEPTOR_INTERACTION.txt",sep = ',',header = T)
gene_emc =colnames(gene_emc)
gene_ADHESION = read.table("data_input/GSEA/heatmap/KEGG_FOCAL_ADHESION.txt",sep = ',',header = T)
gene_ADHESION =colnames(gene_ADHESION)
gene_all = c(gene_emc,gene_ADHESION)


# Matching data
colnames(exp_mRNA)
colnames(exp_mRNA) = gsub('.rep','',gsub('.PT','',colnames(exp_mRNA)))
rownames(new_dat)
table(rownames(new_dat) %in% colnames(exp_mRNA))
exp_mRNA = exp_mRNA[,rownames(new_dat)]
exp_mRNA = as.data.frame(exp_mRNA)
exp_EMC = exp_mRNA[gene_all,]
exp_EMC = na.omit(exp_EMC)


# # # # # # # # # # #
# Immune checkpoint heat map
# # # # # # # # # # # # #
library(pheatmap)
annotation_row = data.frame(
  gene = c(gene_emc,gene_ADHESION),
  group_gene = c(rep("ECM_RECEPTOR_INTERACTION",time = 84),rep("FOCAL_ADHESION",time = 201))
)
annotation_row = annotation_row[annotation_row$gene %in% rownames(exp_EMC),]
table(annotation_row$group_gene) 

# 
sampel_H = new_dat[new_dat$Risk_level == "High",]
data_H = exp_EMC[,rownames(sampel_H)]
data_H = data.frame(
  gene = row.names(data_H),
  mean_H = apply(data_H,1,mean)
)

sampel_L = new_dat[new_dat$Risk_level == "Low",]
data_L = exp_EMC[,rownames(sampel_L)]
data_L = data.frame(
  gene = row.names(data_L),
  mean_L = apply(data_L,1,mean)
)
data_mean = cbind(data_H,data_L)
data_mean = data_mean[,c(2,4)]

# 
pheatmap(data_mean,show_colnames =F,show_rownames = F) 
n=scale(data_mean)
n[n>2]=2 
n[n< -2]= -2
pheatmap(n,show_colnames =F,show_rownames = F,cluster_row = F,cluster_col = F)

#
ac=data.frame(group = c("C1","C2"))
rownames(ac) = colnames(data_mean)
ac

pheatmap(data_mean,show_colnames =F,show_rownames = T,
         annotation_row = annotation_row,
         annotation_col = ac,
         cluster_row = F,cluster_col = F,
         cellheight = 8.7)
pheatmap(data_mean,show_colnames =F,show_rownames = T,
         annotation_col=ac,
         annotation_row = annotation_row,
         cluster_row = F,cluster_col = F,
         cellheight = 8.7,
         filename = '10-F4//1-heatmap_immune.pdf',
         width = 5,
         height = 15)


pheatmap(n,show_colnames =F,show_rownames = T,
         annotation_row = annotation_row,
         annotation_col = ac,
         cluster_row = F,cluster_col = F,
         cellheight = 8.7)
pheatmap(n,show_colnames =F,show_rownames = T,
         annotation_col=ac,
         annotation_row = annotation_row,
         cluster_row = F,cluster_col = F,
         cellheight = 8.7,
         filename = '10-F4//1-heatmap_immune—2.pdf',
         width = 5,
         height = 15)








###################
# 
data_expr  = expr.cu
expr.cu = data_expr[1:30,]


#
require(tidyr)
require(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tibble)
mypalette <- brewer.pal(3,'Set1')

# 
data = t(expr.cu)
dat <- data %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = type,value = Proportion,-Sample)
#
data_C1 =pd.CC[pd.CC$cluster =="C1",] 
dat$Group = ifelse(dat$Sample %in% data_C1$sample,"C1","C2")
table(dat$Group)

#
mypalette <- brewer.pal(3,'Set1')
ggplot(dat,aes(type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Type", y = "Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = c("#377EB8","#E41A1C"))+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")
ggsave('10-F4/boxplot-heatmap_immune-1-30.pdf',height = 10,width = 30)






