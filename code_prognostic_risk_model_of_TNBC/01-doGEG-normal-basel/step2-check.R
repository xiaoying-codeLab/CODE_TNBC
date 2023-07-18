# Settings
# Purpose: PCA and heat map were made to check the differences between normal and basle. Meanwhile, pca of all-sample was also made
# Data source: Rdata
# Save data:

rm(list = ls()) 
options(stringsAsFactors = F)
library(ggsci)

#0- Enter data
load(  file = '00-folder/out-data/0.expr.Rdata')
expr = expr_lnc
group = gp
dim(expr)
expr[1:4,1:4]
table(group)  

##########
#1-PCA
########
# 1. PCA----balse和normal
library("FactoMineR")
library("factoextra")  
dat <- t(expr) 
dat.pca <- PCA(dat , graph = FALSE)
p1=fviz_pca_ind(dat.pca,
                geom.ind = "point", # show points only (nbut not "text")
                col.ind =  group, # color by groups
                palette = c("#00AFBB", "#E7B800"),
                addEllipses = TRUE, # Concentration ellipses
                legend.title = "Tissue Type",
                title = 'Paired Samples',
                ggtheme = theme_minimal()
)
p1 
ggsave(p1,filename = '01-doGEG-normal-basel//out-plot//step1.pca-basel-normal.png')

####################
#2-heatmap  
##################
## tidy data  
dat <- expr
dat[1:4,1:4] 
table(group)

cg=names(tail(sort(apply(dat,1,sd)),1000))
library(pheatmap)
pheatmap(dat[cg,],show_colnames =F,show_rownames = F) 
n=t(scale(t(dat[cg,]))) 
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac= as.data.frame(group)
head(ac)
rownames(ac)=colnames(n) 

## Specify colors
ann_colors = list(
  Tissue_type = c(Normal="#00AFBB", BLIS="#E7B800"))
pheatmap(n,show_colnames =F,show_rownames = F,cluster_cols = T,
         annotation_col=ac,
         annotation_colors = ann_colors)
p2 = pheatmap(n,show_colnames =F,show_rownames = F,cluster_cols = T,
              annotation_col=ac,
              annotation_colors = ann_colors
)
p2
ggsave(p2,filename = '01-doGEG-normal-basel//out-plot/step1.top1000_pheatmap.png') 
dev.off()

#2、Correlation heat map
M=cor(expr)
head(M)
colD=ac
p3 = pheatmap::pheatmap(M,
                        show_colnames =F,show_rownames = F,
                        annotation_colors = ann_colors,
                        annotation_col = colD
)
p3
ggsave(p3,filename = '01-doGEG-normal-basel//out-plot/step1.all_cor.png') 
dev.off()

exprSet=expr[names(sort(apply(expr, 1,mad),decreasing = T)[1:500]),]
dim(exprSet)
M=cor(log2(exprSet+1)) 
head(M)
colD=ac

## Specify colors
ann_colors = list(
  Tissue_type = c(Normal="#00AFBB", Tumor="#E7B800"))
pheatmap::pheatmap(M,
                   show_colnames =F,show_rownames = F,
                   annotation_colors = ann_colors,
                   annotation_col = colD)
p3 = pheatmap::pheatmap(M,
                        show_colnames =F,show_rownames = F,
                        annotation_colors = ann_colors,
                        annotation_col = colD
)
p3
ggsave(p3,filename = '01-doGEG-normal-basel//out-plot/step1.top500_cor.png') 
dev.off()

