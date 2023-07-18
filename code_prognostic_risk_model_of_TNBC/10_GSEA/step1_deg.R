# Objective: To analyze the difference between high risk group and low risk group.

# 设置
rm(list = ls()) 
options(stringsAsFactors = F)
library(ggsci)
library(limma)


# 读入mRNA矩阵
#   mrna表达矩阵
library(data.table)
exp_mRNA = fread( 'data_input/FUSCC//tnbc448.counts.GRCh38.99.txt',header = T, data.table = F)
row.names(exp_mRNA) = exp_mRNA$Geneid
exp_mRNA[1:4,1:4]

exp_mRNA<-exp_mRNA[,-c(1:6)]  #删去前面几行
colnames(exp_mRNA)<-gsub("^...Lib_","",gsub(".bam","",colnames(exp_mRNA)))  #删去后缀，前缀这些
exp_mRNA=log2(edgeR::cpm(exp_mRNA)+1)#标准化矩阵
keep_feature <- rowSums(exp_mRNA>1 ) > 10 #筛选feature
exp_mRNA<-exp_mRNA[keep_feature,]
exp_mRNA[1:4,1:4]
dim(exp_mRNA)
#   基因注释，ENSEMBL换为SYMBOL
library(clusterProfiler)
library(org.Hs.eg.db)
genenames<-bitr(rownames(exp_mRNA),fromType = "ENSEMBL",
                toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
#再处理和转换一下
head(genenames)
length(unique(genenames$ENSEMBL))
length(unique(genenames$SYMBOL))
genenames=genenames[!duplicated(genenames$SYMBOL),]
exp_mRNA<-exp_mRNA[match(genenames$ENSEMBL,rownames(exp_mRNA)),]
rownames(exp_mRNA)=genenames$SYMBOL
exp_mRNA[1:4,1:4]


# 读入高风险和低风险分组
load("04-1-doModel-lasso/out-data/basel-risk.Rdata")
new_dat[1:4,1:4]
colnames(exp_mRNA)
colnames(exp_mRNA)<-gsub("*.rep","",colnames(exp_mRNA)) #删去后缀，前缀这些
table(row.names(new_dat) %in% colnames(exp_mRNA))
exp_mRNA = exp_mRNA[,row.names(new_dat)]

#4、分组信息
library(stringr)
identical(row.names(new_dat),colnames(exp_mRNA))
group = new_dat$Risk_level
table(group)
group <- factor(group, levels = c('Low','High'))
table(group)
group  

# 前面标准化了，直接用limma包吧
library(limma)
design=model.matrix(~factor( group ))
fit=lmFit(exp_mRNA,design)
fit=eBayes(fit)
## 上面是limma包用法的一种方式 
options(digits = 4) #设置全局的数字有效位数为4
#topTable(fit,coef=2,adjust='BH') 
deg = topTable(fit,coef=2,adjust='BH', n=Inf) 
save(deg, file='12_GSEA/out_data/limma_deg.Rdata') 

library(EnhancedVolcano)
EnhancedVolcano(deg,
                lab =  rownames(deg),
                x = 'logFC',
                y = 'P.Value') 

## for heatmap 
##检查上调和下调--是否是正确的。有没有颠换
MFAP5
data_MFAP5 = exp_mRNA[rownames(exp_mRNA) %in% "MFAP5",]
data_MFAP5 = as.data.frame(data_MFAP5)
identical(rownames(data_MFAP5),rownames(new_dat))
data_MFAP5$group = new_dat$Risk_level
colnames(data_MFAP5)
ggplot(dat = data_MFAP5, aes(group,data_MFAP5))+
  geom_boxplot(aes(fill = group))
# 这个基因确实是高风险组组表达高




