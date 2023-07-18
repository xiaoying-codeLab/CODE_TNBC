# Settings
# Objective: To do survival analysis diagram KM
# Data Source:
# Save data:

rm(list = ls())  
options(stringsAsFactors = F) 

library(data.table)
library(dplyr)
library(survminer) 
library(survival)
library(loose.rock)
library(futile.logger) 
library(glmSparseNet)
library(ggrisk)  
library(pheatmap) 
library(ggplot2)
library(ggsci)
library(ggstatsplot)

# 0. Survival data
load(file = '02-doCOX/out-data/2.input_survival_data.Rdata' )
expr = expr.BLIS
survdata = phe_BLIS


expr[1:4,1:4] 
head(survdata)
identical(colnames(expr), rownames(survdata)) 

gene_model <- read.table('04-1-doModel-lasso/out-data/model_lasso_cg_genes.txt')[, 1]
gene_model   # 
# 1. batch km by model genes
## tidy data----

expr[1:4, 1:4]
phe <- survdata
identical(colnames(expr), rownames(phe))
boxplot(phe$time)
head(phe)
kp <- phe$time >0
table(kp)

exprSet <- expr[, kp]
exprSet[1:4, 1:4]
phe <- phe[kp,]
head(phe)
identical(colnames(exprSet), rownames(phe))

# 2. lapply km plot----
model_list <- lapply(gene_model, function(i){
  #i = gene_model[1]
  survival_dat = phe
  gene = as.numeric(exprSet[i,])
  survival_dat$gene = ifelse(gene > median(gene),'high','low')
  table(survival_dat$gene)
  library(survival)
  fit <- survfit(Surv(time, event) ~ gene,
                 data = survival_dat)
  survp = ggsurvplot(fit,data = survival_dat, 
                     legend.title = i, 
                     # legend = 
                     # legend.labs = c('High', 'Low'),
                     pval = T,
                     # pval.method = TRUE,
                     risk.table = TRUE, 
                     risk.table.y.text = F,
                     xlab = "Time in years", 
                     # xlim = c(0, 10), 
                     break.time.by = 1, #
                     size = 1.5, #
                     ggtheme = theme_ggstatsplot(),
                     palette="nejm", #
  )
  return(survp)               
}) 


x=3;y=4
all_plot <- arrange_ggsurvplots(model_list,print = F,ncol =x, nrow = y,
                                risk.table.height = 0.3,
                                surv.plot.height = 0.7)
# all_plot  
x=15;y=20
ggsave(all_plot,filename = '04-1-doModel-lasso/out-plot//step4.choose_lasso_genes_KM.pdf' ,
       height = x, width = y ) 
