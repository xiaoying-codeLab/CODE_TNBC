# Settings
# Purpose: Check the data of normal and basle
# Data source: Rdata
# Save data:

# Create a folder in the current directory
getwd() 
rm(list=ls())
options(stringsAsFactors = F)
library(stringr) 
dir.create('./01-doGEG-normal-basel/out-plot')
dir.create('./01-doGEG-normal-basel/out-data')

#倒入数据
load('00-folder/out-data/0.expr.Rdata')

#检查数据：
dim(expr_lnc)
table(gp)  #139个basel



