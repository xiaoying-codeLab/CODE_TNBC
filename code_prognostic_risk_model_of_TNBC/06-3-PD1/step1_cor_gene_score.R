# Correlation between score and gene expression.

# Calculate correlation
rm(list=ls())
options(stringsAsFactors = F)
library(reshape2)
library(cowplot)
library(plyr)
library(ggplot2)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

# 0. load data----
library(data.table)
rna_exp = fread( 'data_input//tnbc448.counts.GRCh38.99.txt',header = T, data.table = F)
row.names(rna_exp) = rna_exp$Geneid
rna_exp[1:4,1:4]
rna_exp<-rna_exp[,-c(1:6)]  
colnames(rna_exp)<-gsub("^...Lib_","",gsub(".bam","",colnames(rna_exp))) 
colnames(rna_exp)<-gsub("*.rep","",colnames(rna_exp)) 

rna_exp = log2(rna_exp+1)
keep_feature <- rowSums(rna_exp>1 ) > 10 
rna_exp<-rna_exp[keep_feature,]
rna_exp[1:4,1:4]
dim(rna_exp)
rna_exp = as.data.frame(rna_exp)
#id也转换一下吧
library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(unique(rownames(rna_exp)), fromType = "ENSEMBL",
           toType = c( "SYMBOL"),
           OrgDb = org.Hs.eg.db)
head(df)
df=df[!duplicated(df$ENSEMBL),]
df=df[!duplicated(df$SYMBOL),]
rna_exp = rna_exp[df$ENSEMBL,]
identical(rownames(rna_exp),df$ENSEMBL)
rownames(rna_exp) = df$SYMBOL
rna_exp[1:4,1:4]
dim(rna_exp)

# Read the pd-1 gene.
# Read the immunoassay locus gene
gene.ICB = c("VEGFA",
             "CX3CL1",
             "TNFSF4",
             "PDCD1",
             "CD274", 
             'CTLA4')
table(gene.ICB %in% rownames(rna_exp))
gene.ICB[!gene.ICB %in% rownames(rna_exp)]
# Expression matrix of ICB gene
exp.ICB = rna_exp[rownames(rna_exp)  %in% gene.ICB,]
# Read score.
#1. Pour the data
load('05-doEstimate/out-data/estimate_results.Rdata')
dim(scores)
# Merge the expression matrix and the score
exp.ICB = as.data.frame(t(exp.ICB))
exp.ICB = exp.ICB[rownames(scores),]
identical(rownames(scores),rownames(exp.ICB))
data.all = cbind(exp.ICB,scores)


# Start drawing
colnames(data.all)
symbol <- "StromalScore"
topnumber = 6
pointnumber = 20 
tcga_expr <- t(data.all)
tcga_expr[1:3,1:3]
tcga_expr = as.data.frame(tcga_expr)
target.exps <- tcga_expr[symbol,]
other.expr <- tcga_expr[-which(rownames(tcga_expr)==symbol),]

# pearson
sgcor <- as.data.frame(cor(t(other.expr), t(target.exps))) 
sgcor$pval_pearson <- apply(other.expr, 1, function(x) (cor.test(x, t(target.exps))$p.value))

# spearman
sgcor_spearman <- as.data.frame(cor(t(other.expr), t(target.exps), method = "spearman"))
colnames(sgcor_spearman) <- "r_spearman"
sgcor_spearman$pval_spearman <- apply(other.expr, 1, function(x)(cor.test(x, t(target.exps), method = "spearman")$p.value))
#
cors <- cbind(sgcor, sgcor_spearman)
cors$gene <- rownames(other.expr)
head(cors)

newcor <- cors[!(is.na(cors$r_pearson)),]
dim(newcor)

sortcor <- newcor[order(newcor$r_pearson, newcor$r_spearman, decreasing = T),]
topcor <- sortcor[c(1:topnumber, #
                    (nrow(sortcor) - topnumber + 1):nrow(sortcor)),] #
rownames(topcor)

genes <- c(symbol,rownames(topcor))
genesexps <- as.data.frame(t(tcga_expr[genes,]))
sortgenesexps <- genesexps[order(genesexps[,1]),]
samplenum <- nrow(sortgenesexps)

if(is.null(pointnumber)){
  pointnumber=samplenum
}

# plot with special pointnumber by groupping gene expression levels into equally sized bins
group <- as.integer(cut(1:samplenum, breaks=c(seq(from=0.5, to=samplenum, by=samplenum/pointnumber), samplenum+0.5)))
ddf <- data.frame(row.names = 1:pointnumber)

for( i in 1:(1 + topnumber*2)){
  ddf <- cbind(ddf,tapply(sortgenesexps[,i],group,median))
}

colnames(ddf) <- c(symbol,topcor$gene)
mddf <- melt(ddf,id.vars=symbol)

mddf$r <- topcor[mddf$variable,]$r_pearson
mddf$P <- topcor[mddf$variable,]$pval_pearson 
mddf$P = round(mddf$P, 5)

friend <- "EIF3I"
df <- mddf[mddf$variable == friend,]

# Put r and P in italics
rvalue <- as.character(as.expression(substitute(~~italic(r)~"="~rvalue, list(rvalue = format(round(unique(df$r),2), nsmall= 2)))))
pvalue <- as.character(as.expression(substitute(~~italic(P)~"="~pvalue, list(pvalue = format(sprintf("%1.1e", unique(df$P)), nsmall= 2)))))

ggplot(df, aes_string(x=symbol, y="value")) +
  geom_point() +
  ylab(friend) +
  annotate("text", x = 334, y = 350, label = rvalue, parse = TRUE) +
  annotate("text", x = 340, y = 300, label = pvalue, parse = TRUE) +
  geom_smooth(method = "lm", se = F, colour = "#206BB5")

plist <- dlply(mddf, .(variable), function(trig){ggplot(trig, aes_string(x=symbol, y="value")) +
    geom_point() +
    ylab(unique(trig$variable)) +
    ggtitle(paste0("r = ", round(unique(trig$r),2),
                   "\nP = ", sprintf("%1.1e", unique(trig$P)),
                   trig$variable
                   )) +
    geom_smooth(method = "lm", se=F, colour = "#206BB5")})
plist[3]



plist <- dlply(mddf, .(variable), function(trig){ggplot(trig, aes_string(x=symbol, y="value")) +
    geom_point() +
    ylab(unique(trig$variable)) +
    ggtitle(paste0("r = ", round(unique(trig$r),2),
                   "\nP = ", sprintf("%1.1e", unique(trig$P)),
                   paste0("        ",trig$variable)
    )) +
    geom_smooth(method = "lm", se=F, colour = "#206BB5")})
plist[1]

pg <- plot_grid(plotlist = plist, ncol=1, align = "hv")
pg
ggsave("06-PD1/cor_stromalscore.pdf", width = 8, height = 10)

#####
# One picture, one picture save
length(plist)
for (i in 1:length(plist)) {
  # i = 8
  pdf(paste0("06-PD1/cor/",names(plist[i]),"_cor.pdf"))
  plist[i]
  dev.off()
  
}










