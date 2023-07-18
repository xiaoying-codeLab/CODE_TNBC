## Supplementary diagram: Box diagram of tcga's 10lncRNA.
# Read data
# Objective: To compare the difference of 10lncRNA gene between tumor and normal
# Settings
rm(list = ls())  ## 
options(stringsAsFactors = F)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(ggsci)
library(ggpubr)
library(patchwork)
library(reshape2)
library(plyr)
library(stringi)
suppressMessages(library(ggpubr))
suppressMessages(library(dplyr))

gene_overlap <- read.table('03-doOverlap-DEG-COX/out-data//03_overlap_logFC_t_2_DEG_with_Cox_overlap_genes.txt')[, 1]
gene_overlap

#tidy
exp.tcga=fread(file = 'data_input/TCGA//TCGA-BRCA.htseq_counts.tsv.gz',data.table = F)
dim(exp.tcga)  #
exp.tcga[1:4,1:4]
exp.tcga$Ensembl_ID=stri_sub(exp.tcga$Ensembl_ID,1,15)
rownames(exp.tcga)=exp.tcga[,1]
exp.tcga=exp.tcga[,-1]
exp.tcga[1:4,1:4]
colnames(exp.tcga) = stri_sub(colnames(exp.tcga),1,15)  
exp.tcga[1:4,1:4]

subtype =fread('data_input/TCGA//TCGA.BRCA.sampleMap_BRCA_clinicalMatrix',data.table = F)
subtype$sampleID
table(subtype$sampleID %in% colnames(exp.tcga))

colnames(subtype)
tnbc = subtype
tnbc = tnbc[tnbc$ER_Status_nature2012 == "Negative",]
tnbc = tnbc[tnbc$HER2_Final_Status_nature2012 == "Negative",]
tnbc = tnbc[tnbc$PR_Status_nature2012 == "Negative",]
dim(tnbc)

normal = subtype[stri_sub(subtype$sampleID,14) == '11',]
dim(normal) 
sampel = c(tnbc$sampleID,normal$sampleID)
sampel[1:4]
expr.tcga = exp.tcga[,colnames(exp.tcga) %in% sampel]
dim(expr.tcga)
sampel[!sampel %in% colnames(exp.tcga)]
group_list=ifelse(colnames(expr.tcga) %in% tnbc$sampleID,'tnbc','normal')
table(group_list)
lnc.anno = read.table(file = 'data_input/LNCipedia//lncipedia_5_2_ensembl_92_genes.txt',header = T)
dim(lnc.anno)
lnc.anno=lnc.anno[!duplicated(lnc.anno$lncipediaGeneID),]  
lnc.anno$ensemblGeneID = stri_sub(lnc.anno$ensemblGeneID,1,15)
table(rownames(expr.tcga) %in% lnc.anno$ensemblGeneID)


expr.tcga = expr.tcga[rownames(expr.tcga) %in% lnc.anno$ensemblGeneID,]
lnc.anno = lnc.anno[match(rownames(expr.tcga),lnc.anno$ensemblGeneID),]
identical(rownames(expr.tcga),lnc.anno$ensemblGeneID)
rownames(expr.tcga) = lnc.anno$lncipediaGeneID
expr.tcga[1:4,1:4]
dim(expr.tcga)


gene_overlap 
gene_overlap %in% rownames(expr.tcga)
data <- expr.tcga[gene_overlap,]

#3.
table(group_list)
group=as.data.frame(group_list)
rownames(group) = colnames(data)
group$sample = rownames(group)

###################################################
data = t(data)
rownames(data)

library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(ggsci)
library(ggpubr)
library(patchwork)
options(stringsAsFactors = FALSE)



dat <- data %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

#
group_tnbc = group[group$group_list == 'tnbc',]
dat$Group = ifelse(dat$Sample %in% group_tnbc$sample,"tnbc","normal")


library(ggpubr)
library(RColorBrewer)
mypalette <- brewer.pal(3,'Set1')
dat = dat[order(dat$Group,decreasing = F),]

ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Expression") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette)+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..),method = "wilcox.test")
ggsave('04-3-TCGA_check_model/10lnc_tcga_boxpolt.pdf')
##########################################



############################
#
group_tnbc = group[group$group_list == 'tnbc',]
dat$Group = ifelse(dat$Sample %in% group_tnbc$sample,"tnbc","normal")

dat$Group = ifelse(dat$Sample %in% group_tnbc$sample,"normal","tnbc")

library(ggpubr)
library(RColorBrewer)
mypalette <- brewer.pal(3,'Set1')
dat = dat[order(dat$Group,decreasing = F),]

ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Expression") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette)+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..),method = "wilcox.test")
ggsave('04-3-TCGA_check_model/10lnc_tcga_boxpolt2.pdf')
































#4.绘图数据
#4.1加入分组信息
data_new <- data.frame(t(data))
data_new$sample = row.names(data_new)
data_new <- merge(data_new,group,by.x = "sample",by.y = 0)
#4.2融合数据
data_new = melt(data_new)
colnames(data_new) = c("sample","group","gene","expression")
data_new$subject=data_new$sample

#5.加载绘图函数
source("data_input/fustion_viopolt//Function_for_violin_plot.R")

#6.绘制小提琴图
## 6. 绘图
# 6.1 这里注意到原图用的是误差线，这里用步骤三加载的函数，计算一下误差信息
#install.packages("Rmisc")
library(Rmisc)
Data_summary <- summarySE(
  data_new, measurevar="expression", groupvars=c("group","gene"))
head(Data_summary)

# 4.2. 出图
# 这个是我自己写的一个ggplot2的主题，可以自定义修改其中的参数
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) 
}

# 自行调整下面的参数
gene_split_violin <- ggplot(data_new,aes(x= gene,y= expression,fill= group))+
  geom_split_violin(trim= F,color="white",scale = "area") + #绘制分半的小提琴图
  geom_point(data = Data_summary,aes(x= gene, y= expression),pch=19,
             position=position_dodge(0.5),size= 1)+ #绘制均值为点图
  geom_errorbar(data = Data_summary,aes(ymin = expression-ci, ymax= expression+ci), 
                width= 0.05, 
                position= position_dodge(0.5), 
                color="black",
                alpha = 0.8,
                size= 0.5) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ 
  labs(y=("Log2 expression"),x=NULL,title = "Split violin") + 
  theme_bw()+ mytheme +
  stat_compare_means(aes(group = group),
                     label = "p.signif",
                     method = "anova",
                     label.y = max(data_new$expression),
                     hide.ns = T)
gene_split_violin;
ggsave(gene_split_violin,
       filename = "04-3-TCGA_check_model/out_plot//10lncgene_split_violin_tcga.pdf",height = 4.5,width = 6)









