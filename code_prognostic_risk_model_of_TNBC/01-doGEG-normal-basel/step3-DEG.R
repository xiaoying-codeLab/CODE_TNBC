# Settings
# Objective: To analyze the difference between normal and basle.
# Data Source:
# Save data:

rm(list = ls())
options(stringsAsFactors = F)

#################
#1、DEG
#################
load(file = '00-folder/out-data/0.expr.Rdata')
expr =  expr_lnc
group_list = gp
dim(expr)
expr[1:4,1:4]
table(group_list) 
levels(group_list)
library(limma)
design=model.matrix(~factor( group_list ))
fit=lmFit(expr,design)
fit=eBayes(fit)
options(digits = 4) 
#topTable(fit,coef=2,adjust='BH') 
deg = topTable(fit,coef=2,adjust='BH', n=Inf) 
save(deg, file='01-doGEG-normal-basel//out-data//1.limma_deg.Rdata') 


library(EnhancedVolcano)
EnhancedVolcano(deg,
                lab =  rownames(deg),
                x = 'logFC',
                y = 'P.Value') 
ggsave( filename = '01-doGEG-normal-basel/out-plot/step1.volcano.pdf')

## for heatmap 
if(T){ 
  dat= expr
  # 每次都要检测数据
  dat[1:4,1:4]
  table(group_list) 
  cg=c( head(rownames(deg[order(deg$logFC),]),10) ,
        tail(rownames(deg[order(deg$logFC),]),10) )
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F) 
  n=t(scale(t(dat[cg,])))
  
  n[n>2]=2
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F )
  ac=data.frame(g=group_list)
  rownames(ac)=colnames(n) 
  pheatmap(n,show_colnames =F,
           show_rownames = T,
           cluster_cols = T, 
           annotation_col=ac)
  pheatmap(n,show_colnames =F,
           show_rownames = T,
           cluster_cols = T, 
           annotation_col=ac,
           filename = '01-doGEG-normal-basel/out-plot//step1.choose_top_heatmap_top20_DEG.pdf') #列名注释信息为ac即分组信息
  pheatmap(n,show_colnames =F,show_rownames = F,annotation_col=ac,
           filename = '01-doGEG-normal-basel/step1.choose_top_heatmap_top20_DEG_noSYMBOL.pdf')
  
}
dev.off()

