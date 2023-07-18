# Settings
# Objective: To do immune pathway enrichment, ssGSVA
# Data Source:
# Save data:
getwd()
dir.create('./06-Immune-ssGSEA//out-plot')
dir.create('./06-Immune-ssGSEA///out-data')

rm(list = ls())
options(stringsAsFactors = F)
library(GSVA)
library(limma)
library(GSEABase)
#1. Pour the data. Representation matrix
load('05-doEstimate/out-data/exp.basel.Rdata')
dim(exp.basel)
exp.basel[1:4,1:4]
exp = exp.basel


# 2. gmt
gmtFile <- "06-Immune-ssGSEA/immune.gmt"   
dimnames <- list(rownames(exp),colnames(exp))
mat <- matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
mat <- avereps(mat)
mat <- mat[rowMeans(mat)>0,]
geneSet=getGmt(gmtFile, 
               geneIdType=SymbolIdentifier())

# 2. do gsva----
ssgseaScore <- gsva(mat, geneSet, 
                    method='ssgsea', 
                    kcdf='Gaussian', abs.ranking=TRUE)

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}

ssgseaOut <- normalize(ssgseaScore)
ssgseaOut <- rbind(id=colnames(ssgseaOut),ssgseaOut)

# 3. save data----
write.table(ssgseaOut,file="06-Immune-ssGSEA/out-data//step1_ssgseaOut.txt",
            sep="\t",quote=F,col.names=F)






