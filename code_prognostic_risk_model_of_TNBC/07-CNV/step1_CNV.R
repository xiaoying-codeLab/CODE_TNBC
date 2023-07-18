# Drawing of CNV

# Settings
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(dplyr)

# Read CNV data
rt=read.table("07-CNV/TCGA-BRCA.gistic.tsv.gz", header=T, sep="\t", check.names=F, row.names=1)
rownames(rt) = substr(rownames(rt),1,15)

# Then convert the id to lnc
# Match lnc
#1, matches the transformation table provided by the LNCipedia database to convert to Ensembl
lnc.anno = read.table(file = 'data_input/LNCipedia//lncipedia_5_2_ensembl_92_genes.txt',header = T)
dim(lnc.anno)
lnc.anno=lnc.anno[!duplicated(lnc.anno$lncipediaGeneID),]   
lnc.anno$ensemblGeneID = stri_sub(lnc.anno$ensemblGeneID,1,15)#
table(rownames(rt) %in% lnc.anno$ensemblGeneID)


expr.tcga = expr.tcga[rownames(expr.tcga) %in% lnc.anno$ensemblGeneID,]
lnc.anno = lnc.anno[match(rownames(expr.tcga),lnc.anno$ensemblGeneID),]
identical(rownames(expr.tcga),lnc.anno$ensemblGeneID)
rownames(expr.tcga) = lnc.anno$lncipediaGeneID
expr.tcga[1:4,1:4]
dim(expr.tcga)

# Conversion id
# # # # # # # # # # # # # # # # # # #
# Gene annotation, ENSEMBL changed to SYMBOL
# # # # # # # # # # # # # # # # # # #
library(clusterProfiler)
library(org.Hs.eg.db)
genenames<-bitr(rownames(rt),fromType = "ENSEMBL",
                toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
head(genenames)
length(unique(genenames$ENSEMBL))
length(unique(genenames$SYMBOL))
genenames=genenames[!duplicated(genenames$SYMBOL),]
rt<-rt[match(genenames$ENSEMBL,rownames(rt)),]
rownames(rt)=genenames$SYMBOL
rt[1:4,1:4]
rownames(rt)[1:50]
dim(rt)
rt = as.data.frame(rt)

# And then screen for triple-negative breast cancer.
# Filter directly from clinicalMatrix
subtype =fread('07-CNV//TCGA.BRCA.sampleMap_BRCA_clinicalMatrix',data.table = F)
colnames(subtype)
subtype = subtype[subtype$ER_Status_nature2012 == "Negative",]
subtype = subtype[subtype$HER2_Final_Status_nature2012 == "Negative",]
subtype = subtype[subtype$PR_Status_nature2012 == "Negative",]
dim(subtype)
subtype$sampleID

colnames(rt)
colnames(rt) = substr(colnames(rt),1,15) 
table(subtype$sampleID %in% colnames(rt))
rt = rt[colnames(rt) %in% subtype$sampleID]
dim(rt)




# And then there are the nine iron-death genes.
# But choose one of them first
gene_overlap <- read.table('03-doOverlap-DEG-COX/out-data//03_overlap_logFC_t_2_DEG_with_Cox_overlap_genes.txt')[, 1]
gene_overlap
gene_overlap %in% rownames(rt)
rt = rt[gene_fe,]
write.csv(rt,file = "07-CNV/CNV.csv")


##### 
GAIN=rowSums(rt> 0)       
LOSS=rowSums(rt< 0)       
GAIN=GAIN/ncol(rt)*100      
LOSS=LOSS/ncol(rt)*100      
data=cbind(GAIN, LOSS)
data=data[order(data[,"GAIN"],decreasing = T),]
write.csv(data,file = "07-CNV/CNVfreq.csv")

#
data.max = apply(data, 1, max)
pdf(file="07-CNV/CNVfreq.pdf", width=9, height=6)
cex=1.3
par(cex.lab=cex, cex.axis=cex, font.axis=2, las=1, xpd=T)
bar=barplot(data.max, col="grey80", border=NA,
            xlab="", ylab="CNV.frequency(%)", space=1.5,
            xaxt="n", ylim=c(0,1.2*max(data.max)))
points(bar,data[,"GAIN"], pch=20, col=2, cex=3)
points(bar,data[,"LOSS"], pch=20, col=3, cex=3)
legend("top", legend=c('GAIN','LOSS'), col=2:3, pch=20, bty="n", cex=2, ncol=2)
par(srt=45)
text(bar, par('usr')[3]-0.2, rownames(data), adj=1)
dev.off()




