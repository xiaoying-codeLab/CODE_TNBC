# Settings
# Objective: To cross 37 genes and do a batch of km survival
# Data Source:
# Save data:

rm(list = ls())  
options(stringsAsFactors = F) 

library(survminer) 
library(survival)
library(ggstatsplot)

load(file = '02-doCOX/out-data//2.input_survival_data.Rdata' )
expr = expr.BLIS
phe = phe_BLIS
expr[1:4,1:4] 
head(phe)
identical(colnames(expr), rownames(phe))
gene_overlap <- read.table('03-doOverlap-DEG-COX/out-data//03_overlap_logFC_t_2_DEG_with_Cox_overlap_genes.txt')[, 1]
gene_overlap
# 1. batch km by overlap genes
## tidy data---- 
exprSet <- expr
exprSet[1:4, 1:4]
boxplot(phe$time)
kp <- phe$time >0
table(kp)
exprSet <- exprSet[, kp]
exprSet[1:4, 1:4]
phe <- phe[kp,]
head(phe)
identical(colnames(exprSet), rownames(phe))
for (i in 1:10) {
  # i = 1
  gene_name = gene_overlap[i]
  survival_dat = phe
  gene = as.numeric(exprSet[gene_name,])
  survival_dat$gene = ifelse(gene > median(gene),'high','low')
  table(survival_dat$gene)
  library(survival)
  fit <- survfit(Surv(time, event) ~ gene,
                 data = survival_dat)
  
  survp = ggsurvplot(fit,data = survival_dat, 
                     legend.title = gene_name, 
                     # legend = "right",
                     legend.labs = c('Risk_High', 'Risk_Low'),
                     legend = c(0.7,0.2),
                     font.legend  = 15,
                     pval = T, 
                     risk.table = F, 
                     risk.table.y.text = F,
                     xlab = "Time in years", 
                     # xlim = c(0, 10), 
                     break.time.by = 1, 
                     size = 1,
                     ggtheme = theme_ggstatsplot(),
                     palette="nejm",
                     font.x = 9,  
                     font.y = 9)
  print(survp)
  pdf( paste0("03-doOverlap-DEG-COX//out-plot/",gene_name,"_multicox_KM.pdf"), onefile = F)
  print(survp)
  dev.off()
  
}



# 2. lapply km plot----
overlap_list <- lapply(gene_overlap, function(i){
  #i = gene_overlap[1]
  survival_dat = phe
  gene = as.numeric(exprSet[i,])
  survival_dat$gene = ifelse(gene > median(gene),'high','low')
  table(survival_dat$gene)
  library(survival)
  fit <- survfit(Surv(time, event) ~ gene,
                 data = survival_dat)
  
  survp = ggsurvplot(fit,data = survival_dat, 
                     legend.title = i, 
                     # legend = "top",
                     # legend.labs = c('High', 'Low'),
                     pval = T, 
                     # pval.method = TRUE,
                     risk.table = TRUE, 
                     risk.table.y.text = F,
                     xlab = "Time in years", 
                     # xlim = c(0, 10), 
                     break.time.by = 1, 
                     size = 1.5, 
                     ggtheme = theme_ggstatsplot(),
                     palette="nejm", 
  )
  return(survp)               
}) 

x=2;y=5
all_plot <- arrange_ggsurvplots(overlap_list,print = F,ncol =x, nrow = y,
                                risk.table.height = 0.3,
                                surv.plot.height = 0.7)
#all_plot  
x=35;y=15   
ggsave(all_plot,filename = '03-doOverlap-DEG-COX/out-plot//step3.choose_overlap_genes_KM.pdf' ,
       height = x, width = y ) 

