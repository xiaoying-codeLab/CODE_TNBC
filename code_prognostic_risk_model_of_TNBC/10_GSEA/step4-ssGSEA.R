rm(list = ls())  ## 
options(stringsAsFactors = F)

library(Seurat)
library(gplots)
library(ggplot2) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA) # BiocManager::install('GSVA')
library(GSEABase)

library(Seurat)
library(gplots)
library(ggplot2)


# 
load('04-1-doModel-lasso/out-data//basel-risk.Rdata')
table(new_dat$Risk_level)
group = 
#
# 
#####################
#   
library(data.table)
exp_mRNA = fread( 'data_input/FUSCC//tnbc448.counts.GRCh38.99.txt',header = T, data.table = F)
row.names(exp_mRNA) = exp_mRNA$Geneid
exp_mRNA[1:4,1:4]
exp_mRNA<-exp_mRNA[,-c(1:6)]  #
colnames(exp_mRNA)<-gsub("^...Lib_","",gsub(".bam","",colnames(exp_mRNA)))  #
exp_mRNA=log2(edgeR::cpm(exp_mRNA)+1)#
keep_feature <- rowSums(exp_mRNA>1 ) > 10 #
exp_mRNA<-exp_mRNA[keep_feature,]
exp_mRNA[1:4,1:4]
dim(exp_mRNA)
#   
library(clusterProfiler)
library(org.Hs.eg.db)
genenames<-bitr(rownames(exp_mRNA),fromType = "ENSEMBL",
                toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
#
head(genenames)
length(unique(genenames$ENSEMBL))
length(unique(genenames$SYMBOL))
genenames=genenames[!duplicated(genenames$SYMBOL),]
exp_mRNA<-exp_mRNA[match(genenames$ENSEMBL,rownames(exp_mRNA)),]
rownames(exp_mRNA)=genenames$SYMBOL
exp_mRNA[1:4,1:4]
colnames(exp_mRNA)<-gsub("*.rep","",colnames(exp_mRNA)) #
exp_mRNA = exp_mRNA[,rownames(new_dat)]
dim(exp_mRNA)
exp_mRNA[1:4,1:4]

################
exp_mRNA[1:4,1:4]
dat = exp_mRNA
group_list = new_dat
table(group_list)

# https://www.wikipathways.org/index.php/Pathway:WP113
load('12_GSEA/MSigDB-v5.2/human/human_c2_v5p2.rdata')
Hs.c2
Hs.c2 = lapply(Hs.c2, function(x) bitr(geneID = x,fromType = 'ENTREZID',
                                     toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')[,2])
#geneset <- getGmt(Mm.H)
gsc <- GeneSetCollection(mapply(function(geneIds, keggId) {
  GeneSet(geneIds, geneIdType=EntrezIdentifier(),
          collectionType=KEGGCollection(keggId),
          setName=keggId)
}, Hs.c2, names(Hs.c2)))
gsc
geneset <- gsc

# method=c("gsva", "ssgsea", "zscore", "plage"),
# kcdf=c("Gaussian", "Poisson", "none"),
# mx.diff=FALSE: ES is calculated as the maximum distance of the random walk from 0. 
# mx.diff=TRUE (default): ES is calculated as the magnitude difference between the largest positive and negative random walk deviations.

es.max <- gsva(dat , geneset, 
               mx.diff=FALSE, verbose=FALSE, 
               parallel.sz=8)
save(es.max,file = '12_GSEA/gsva_hallmark.Rdata')
load(file = '12_GSEA/gsva_hallmark.Rdata') 
es.max[1:4,1:4]
ac = as.data.frame(group_list$Risk_level)
rownames(ac) = rownames(group_list)

ac=data.frame(group =substring(  colnames(ensembl_matrix) ,1,1),
               replicate=substring(  colnames(ensembl_matrix) ,2,2))
head(ac) 

identical(rownames(ac),colnames( es.max ))
mat=es.max
mat[1:4,1:4]


gseaNb(object = es.max,
       geneSetID = 'KEGG_ECM_RECEPTOR_INTERACTION')

####################################################################
gseaRes <- GSEA(geneList = geneList,
                TERM2GENE = gmt,
                minGSSize    = 10,
                maxGSSize    = 500,
                pvalueCutoff = 1,
                pAdjustMethod = "BH",
                verbose      = FALSE)





pheatmap::pheatmap(mat,
                   show_rownames = T,show_colnames = F,
                   annotation_col = ac)
pheatmap::pheatmap(mat,
                   show_rownames = T,show_colnames = F,
                   annotation_col = ac,
                   filename = 'gsva_hap_hallmark.pdf',width = 10,height = 10
)
dev.off()

library(ggpubr)
es.max[1:4,1:4]
head(ac)
boxplot(es.max['KEGG_ECM_RECEPTOR_INTERACTION',] ~ ac$group)
df=data.frame(value=es.max['KEGG_ECM_RECEPTOR_INTERACTION',],
              group= ac$group)
ggboxplot(df, "group", "value",
          color = "group", #palette =c("#00AFBB", "#E7B800", "#FC4E07"),
          add = "jitter", shape = "group")
ggsave('gsea_value_for_HALLMARK_TGF_BETA_SIGNALING.pdf',
       width = 3,height = 3)


######################################################################
# devtools::install_github("junjunlab/GseaVis")

library(GseaVis)

# load data
test_data <- system.file("extdata", "gseaRes.RDS", package = "GseaVis")
gseaRes <- readRDS(test_data)
gseaRes

gseaNb(object = gseaRes,
       geneSetID = 'GOBP_NUCLEOSIDE_DIPHOSPHATE_METABOLIC_PROCESS')
