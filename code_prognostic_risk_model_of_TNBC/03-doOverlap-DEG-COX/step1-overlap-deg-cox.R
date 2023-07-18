# Settings
# Objective: To obtain the gene intersection of deg and cox
# Data Source:
# Save data:
# Create a folder in the current directory
getwd()
dir.create('./03-doOverlap-DEG-COX///out-plot')
dir.create('./03-doOverlap-DEG-COX///out-data')


rm(list=ls())
options(stringsAsFactors = F)
library(ggvenn)
library(patchwork)

# 0. load data----
## basle和normal的 DEG
load(file='01-doGEG-normal-basel/out-data//1.limma_deg.Rdata') 
head(deg)

load(file = '02-doCOX/out-data//2.batch_cox_results.Rdata') 
table(cox_results[,4] < 0.05)  

# Set the threshold for logFC
thred <- c(0.5, 1, 1.5, 2)
lapply(1:4, function(i){
  try({
    # i = 2
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
    # 1. overlap----
    dat_overlap <- Reduce(intersect, 
                          list(dat_paired, dat_cox));
    length(dat_overlap)
    # 2. save overlap txt----
    write.table(sort(dat_overlap),
                file = file.path('03-doOverlap-DEG-COX/out-data//',
                                 paste0('03_overlap_logFC_t_',logFC_t, '_DEG_with_Cox_overlap_genes.txt')),
                row.names = F,col.names = F,quote = F)
    
  }, silent = T)
  # 3. visualization
  ## 2.1 venn----
  venn_list <- list(`basel up` = dat_paired_up,
                    `basel down` = dat_paired_down,
                    `Batch Cox` = dat_cox,
                    `basel` = dat_paired
  )
  names(venn_list)
  ## up
  vp_up <- ggvenn(venn_list,
                  c('basel up','Batch Cox'),
                  fill_color = c("#B4679F", "#EC9C29"),
                  set_name_size = 5,
                  stroke_size = 0.5) 
  vp_up 
  ## down
  vp_down <- ggvenn(venn_list,
                    c('basel down','Batch Cox'),
                    fill_color = c("#97915D", "#679880"),
                    set_name_size = 5,
                    stroke_size = 0.5) 
  vp_down
  ## all
  vp_all <- ggvenn(venn_list,
                   c("basel", 'Batch Cox'),
                   fill_color = c("#97915D", "#679880", '#AA1F24'),
                   set_name_size = 5,
                   stroke_size = 0.5)
  vp_all
  vp_up+vp_down+vp_all
  ggsave(file.path('03-doOverlap-DEG-COX/out-plot/',paste0('step3.logFC_t_',logFC_t,'.png')))
  
})


