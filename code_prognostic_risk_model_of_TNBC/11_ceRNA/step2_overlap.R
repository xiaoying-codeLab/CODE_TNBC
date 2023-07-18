
data_lncRNA = read.csv("../F4_ceRNA/lnc_mirna.csv")
data_mrna = read.csv("../F4_ceRNA/miRWalk_miRNA_Targets3.csv")
data_mrna$validated

data_mrna$validated = ifelse(data_mrna$validated %in% data_mrna$validated[grep('MIR',data_mrna$validated)],"1","0")
table(data_mrna$validated)
colnames(data_mrna)
data_mrna$sum = rowSums(data_mrna[,c(19,20,21)])
table(data_mrna$sum)
data_mrna = data_mrna[data_mrna$sum >0,]





class(data_mrna[,19])
data_mrna[,19] = as.numeric(data_mrna[,c(19)])
data_mrna[,c(19,20,21)] = as.numeric(data_mrna[,c(19,20,21)])


data_mrna$validated[grep('MIR',data_mrna$validated)]




table(data_lncRNA$miRNA %in% data_mrna$锘縨irnaid)


data_lncRNA = data_lncRNA[data_lncRNA$miRNA %in% data_mrna$锘縨irnaid,]
data_mrna = data_mrna[data_mrna$锘縨irnaid %in% data_lncRNA$miRNA,]

write.csv(data_lncRNA,file = "data_lncRNA.csv")

write.csv(data_mrna,file = "data_mrna.csv")







