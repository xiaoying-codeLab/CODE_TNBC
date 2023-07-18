# Objective: Visualization of drug sensitivity

rm(list = ls())  
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)


library(data.table)
testPtype <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F,header = T)
testPtype[1:4, 1:4]
rownames(testPtype) = testPtype$V1
testPtype =  testPtype[,-1]
testPtype[1:4, 1:4]



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
df  = read.csv('04-dug/C1-C2-dug2.csv')
table(df$p < 0.05)
kp.dug = df[df$p < 0.05,]
data.dug = testPtype[,kp.dug$.id]
data.dug <- data.dug[,apply(data.dug, 2, sum) > 500] 

rownames(data.dug.p)  = data.dug.p$.id
data.dug.kp = data.dug.p[colnames(data.dug),]
write.csv(data.dug.kp,file = "04-dug/dug3.csv")
kp_dug = c("BMS-536924_1091",
           "Dihydrorotenone_1827",
           "WIKI4_1940",
           "GSK343_1627",
           "IRAK4_4710_1716",
           "GSK1904529A_1093",
           "Lapatinib_1558",
           "Leflunomide_1578",
           "Nilotinib_1013",
           "YK-4-279_1239",
           "GSK591_2110",
           "OSI-027_1594",
           "Picolinici-acid_1635",
           "GSK2578215A_1927",
           "AZD5991_1720",
           "Sorafenib_1085",
           "Temozolomide_1375",
           "SB505124_1194",
           "Cisplatin_1005",
           "VSP34_8731_1734"
           
)
table(kp_dug %in% colnames(testPtype))
exp.dug = testPtype[,colnames(testPtype)  %in% kp_dug]

write.csv(exp.dug,file = "04-dug/dug.csv")
write.csv(pd.CC,file = "04-dug/分组.csv")

data = exp.dug


dat <- data %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
pd.c1 = pd.CC[pd.CC$cluster == 'C1',]
dat$Group = ifelse(dat$Sample %in% pd.c1$sample,"C1","C2")
library(ggpubr)
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")
ggsave('04-dug//dug1.pdf')
kp_dug = c("IRAK4_4710_1716",
           "GSK1904529A_1093",
           "Leflunomide_1578",
           "OSI-027_1594",
           "Picolinici-acid_1635",
           "GSK2578215A_1927",
           "AZD5991_1720",
           "Temozolomide_1375"
)
table(kp_dug %in% colnames(testPtype))
exp.dug = testPtype[,colnames(testPtype)  %in% kp_dug]
data = exp.dug


dat <- data %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)
pd.c1 = pd.CC[pd.CC$cluster == 'C1',]
dat$Group = ifelse(dat$Sample %in% pd.c1$sample,"C1","C2")
library(ggpubr)
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_violin()(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "DUG", y = "score") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ 
  stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")
ggsave('04-dug//dug2.pdf')
















