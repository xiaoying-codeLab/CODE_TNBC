# Settings
# Objective: To prepare data representation matrix and survival data for survival analysis
# Data source: Rdata
# Save data:


rm(list=ls())
options(stringsAsFactors = F)
library(stringr)  
library(data.table) 

# Import data
load('00-folder/out-data//0.expr.Rdata')
expr =  expr_lnc
group = gp
dim(expr)
expr[1:4,1:4]
table(group)  

## survival information
load('00-folder/out-data//0.survival.Rdata')
head(phe_BLIS)
dim(phe_BLIS)
expr.BLIS = expr[,gp == 'BLIS']
dim(expr.BLIS)
expr.BLIS[1:4,1:4]
table(colnames(expr.BLIS) %in% row.names(phe_BLIS) )
expr.BLIS[1:4,1:4]
head(phe_BLIS)
save(expr.BLIS, phe_BLIS, 
     file = '02-doCOX/out-data//2.input_survival_data.Rdata')








