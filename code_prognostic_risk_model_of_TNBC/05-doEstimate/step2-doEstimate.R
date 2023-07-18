# Settings
# Objective: To analyze the immunoinfiltration. -- From the immune perspective
# Data Source:
# Save data:

rm(list=ls())
options(stringsAsFactors = F)
library(stringr) 
library(limma)
library(estimate)


#1、load 
load('05-doEstimate/out-data/exp.basel.Rdata')
dim(exp.basel)
exp.basel[1:4,1:4]
rownames(exp.basel)[1:50]


# 2. do estimate----
estimate_RNAseq <- function(RNAseq_logCPM,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(RNAseq_logCPM,file = input.f,sep = '\t',quote = F)
  
  library(estimate)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="illumina")
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  scores=data.frame(  scores)
  return(scores)
}

pro <- 'CheckCode'
dat <- exp.basel
scores <- estimate_RNAseq(dat, pro)
head(scores)
save(scores, file = '05-doEstimate/out-data/estimate_results.Rdata')




