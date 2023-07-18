# Objective: Drug analysis

# Installation package:
library('oncoPredict')

rm(list = ls())  
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)

dir= "10_drug/DataFiles/Training Data/"
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 


# Read your own expression matrix
load('05-doEstimate/out-data/exp.basel.Rdata')
dim(exp.basel)
exp.basel[1:4,1:4]
rownames(exp.basel)
# testExpr = exp.basel[1:40,1:40]
testExpr = exp.basel
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

library(data.table)
testPtype <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F,header = T)
testPtype[1:4, 1:4]
rownames(testPtype) = testPtype$V1
testPtype =  testPtype[,-1]
testPtype[1:4, 1:4]

sum = apply(testPtype,2,mean)
table(sum > 30)
testPtype = testPtype[,sum > 30]

####################
library(foreign)
library(car) 
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
load('04-1-doModel-lasso/out-data/basel-risk.Rdata')
testPtype = testPtype[rownames(new_dat),]
identical(rownames(new_dat),rownames(testPtype))
testPtype$group = new_dat$Risk_level
dim(testPtype)
pval_pearson <- apply(testPtype[1:88], 2, function(x) 
  (wilcox.test(formula = x ~ group, data = testPtype))[3])
library (plyr)
df <- ldply (pval_pearson, data.frame)

options(digits = 4)
df$p.value
df$p= round(df$p.value, 3)
table(df$p < 0.05)


write.csv(df,file = '10_drug/dug_p_1.value.csv')
kp_dug = c("Cisplatin",
           "Paclitaxel",
           "Epirubicin",
           "Nilotinib",
           "Dasatinib",
           "Cytarabine",
           "Sinularin",
           "Leflunomide",
           "XAV939",
           "Fulvestrant",
           "Entinostat",
           "PRIMA-1MET"
           
)
table(kp_dug %in% colnames(testPtype))

colnames(testPtype) <- gsub('_****','',colnames(testPtype)) 
exp.dug = testPtype[,colnames(testPtype)  %in% kp_dug]
write.csv(exp.dug,file = "04-dug/dug.csv")
write.csv(pd.CC,file = "04-dug/分组.csv")
data = exp.dug
dat <- data %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
pd.c1 = pd.CC[pd.CC$cluster == 'C1',]
dat$Group = ifelse(dat$Sample %in% pd.c1$sample,"C1","C2")
write.csv(testPtype,file = '10_drug/dug_value.csv')
