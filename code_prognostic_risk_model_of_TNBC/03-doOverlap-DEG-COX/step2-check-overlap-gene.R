# Settings
# Objective: To examine the intersection of step1's gene set - to make a heat map
# Data Source:
# Save data:

rm(list = ls()) 
options(stringsAsFactors = F)

# Import data
load('00-folder/out-data//0.expr.Rdata')
expr = expr_lnc
group = gp
dim(expr)
expr[1:4,1:4]
table(group)  

## gene_overlap
gene_overlap <- read.table('03-doOverlap-DEG-COX/out-data//03_overlap_logFC_t_2_DEG_with_Cox_overlap_genes.txt')[, 1]
gene_overlap
library(VennDiagram)
library(RColorBrewer)
library(ggvenn)
load(file='01-doGEG-normal-basel/out-data//1.limma_deg.Rdata') 
head(deg)
load(file = '02-doCOX/out-data//2.batch_cox_results.Rdata') 
table(cox_results[,4] < 0.05)  


thred <- c(0.5, 1, 1.5, 2)
i = 4
logFC_t <- thred[i]
head(deg)
deg_paired <- deg
deg_paired$Change <- as.factor(ifelse(deg_paired$P.Value < 0.05 & abs(deg_paired$logFC) > logFC_t,
                                      ifelse(deg_paired$logFC > logFC_t ,'Up','Down'),'Stable'))
head(deg_paired)
table(deg_paired$Change)
dat_paired <- rownames(deg_paired[deg_paired$Change != 'Stable', ]); length(dat_paired)
dat_paired_up <- rownames(deg_paired[deg_paired$Change == 'Up', ]); length(dat_paired_up)
dat_paired_down <- rownames(deg_paired[deg_paired$Change == 'Down', ]); length(dat_paired_down)


## batch cox----
table(cox_results[,4] < 0.05)  
dat_cox <- rownames(cox_results[cox_results[, 4] < 0.05, ]);  
length(dat_cox)

venn_list <- list(`basel up` = dat_paired_up,
                  `basel down` = dat_paired_down,
                  `Cox-gene` = dat_cox,
                  `Diff-gene` = dat_paired
)
vp_all <- ggvenn(venn_list,
                 c("Diff-gene", 'Cox-gene'),
                 fill_color = brewer.pal(3,'Set1'),
                 set_name_size = 5,
                 stroke_size = 0.5)
vp_all
ggsave(vp_all,filename = "03-doOverlap-DEG-COX/Overlap.pdf")

# heatmap
library(pheatmap)
n=t(scale(t(expr[gene_overlap,]))) 
n[n>2]=2 
n[n< -2]= -2
n[1:4,1:4]
pheatmap(n,show_colnames =F,show_rownames = F)
ac=as.data.frame(group)
head(ac)
rownames(ac)=colnames(n) 

p = pheatmap(n,show_colnames =F,show_rownames = T,cluster_cols = F,
             annotation_col=ac
)
p
ggsave(p,filename = '03-doOverlap-DEG-COX//overlap_genes_heatmap.pdf')


