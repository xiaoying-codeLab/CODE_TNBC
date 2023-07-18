# go and kegg analysis
# Purpose: go and kegg
# Data: deg difference analysis


# Read data
rm(list = ls()) 
options(stringsAsFactors = F)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)


gene_overlap <- read.table('03-doOverlap-DEG-COX/out-data//03_overlap_logFC_t_1.5_DEG_with_Cox_overlap_genes.txt')[, 1]
gene_overlap

gene_overlap <- bitr(gene_overlap, fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
gene_overlap = gene_overlap$ENTREZID
go <- enrichGO(gene_overlap, OrgDb = "org.Hs.eg.db", ont="all") 
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free") 

kk.up <- enrichKEGG(gene         = gene_up,
                    organism     = 'hsa',
                    #universe     = gene_all,
                    pvalueCutoff = 0.9,
                    qvalueCutoff =0.9)
head(kk.up)[,1:6]
kk=kk.up
g_kegg = dotplot(kk)
g_kegg
ggsave(g_kegg,filename = "up_kegg.pdf")







## KEGG pathway analysis
### 做KEGG数据集超几何分布检验分析，重点在结果的可视化及生物学意义的理解。
run_kegg <- function(gene_up,gene_down,geneList=F,pro='test'){
  gene_up=unique(gene_up)
  gene_down=unique(gene_down)
  gene_diff=unique(c(gene_up,gene_down))
  ###   over-representation test
  # 下面把3个基因集分开做超几何分布检验
  # 首先是上调基因集。
  kk.up <- enrichKEGG(gene         = gene_up,
                      organism     = 'hsa',
                      #universe     = gene_all,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
  head(kk.up)[,1:6]
  kk=kk.up
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'_kk.up.csv'))
  
  # 首先是下调基因集。
  kk.down <- enrichKEGG(gene         =  gene_down,
                        organism     = 'hsa',
                        #universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9)
  head(kk.down)[,1:6]
  kk=kk.down
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'_kk.down.csv'))
  
  # 最后是上下调合并后的基因集。
  kk.diff <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
  head(kk.diff)[,1:6]
  kk=kk.diff
  dotplot(kk)
  kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
  write.csv(kk@result,paste0(pro,'_kk.diff.csv'))
  
  
  kegg_diff_dt <- as.data.frame(kk.diff)
  kegg_down_dt <- as.data.frame(kk.down)
  kegg_up_dt <- as.data.frame(kk.up)
  down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.01,];down_kegg$group=-1
  up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.01,];up_kegg$group=1
  #画图设置, 这个图很丑，大家可以自行修改。
  g_kegg=kegg_plot(up_kegg,down_kegg)
  print(g_kegg)
  ggsave(g_kegg,filename = 'kegg_up_down.pdf',width = 7,height = 10)
  
  ggsave(g_kegg,filename = paste0(pro,'_kegg_up_down.pdf') )
  
  if(geneList){
    ###  GSEA 
    ## GSEA算法跟上面的使用差异基因集做超几何分布检验不一样。
    kk_gse <- gseKEGG(geneList     = geneList,
                      organism     = 'hsa',
                      nPerm        = 1000,
                      minGSSize    = 20,
                      pvalueCutoff = 0.9,
                      verbose      = FALSE)
    head(kk_gse)[,1:6]
    gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
    gseaplot(kk_gse, 'hsa04110',title = 'Cell cycle') 
    kk=DOSE::setReadable(kk_gse, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
    tmp=kk@result
    write.csv(kk@result,paste0(pro,'_kegg.gsea.csv'))
    
    
    # 这里找不到显著下调的通路，可以选择调整阈值，或者其它。
    down_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore < 0,];down_kegg$group=-1
    up_kegg<-kk_gse[kk_gse$pvalue<0.05 & kk_gse$enrichmentScore > 0,];up_kegg$group=1
    
    g_kegg=kegg_plot(up_kegg,down_kegg)
    print(g_kegg)
    ggsave(g_kegg,filename = paste0(pro,'_kegg_gsea.pdf'))
    # 
  }
  
}

kegg_plot <- function(up_kegg,down_kegg){
  dat=rbind(up_kegg,down_kegg)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  
  dat=dat[order(dat$pvalue,decreasing = F),]
  
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing = F)), y=pvalue, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="log10P-value") +
    coord_flip() + theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
    ggtitle("Pathway Enrichment") 
}

load(file = '../../01-doGEG-normal-basel/out-data/1.limma_deg.Rdata')
head(deg)
## 不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
logFC_t=0.58
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)
head(deg)
deg$symbol=rownames(deg)




library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
head(df)
DEG=deg
head(DEG)

DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
head(DEG)
save(DEG,file = 'anno_DEG.Rdata')


gene_up= DEG[DEG$g == 'UP','ENTREZID'] 
gene_down=DEG[DEG$g == 'DOWN','ENTREZID'] 
gene_diff=c(gene_up,gene_down)
gene_all=as.character(DEG[ ,'ENTREZID'] )
data(geneList, package="DOSE")
head(geneList)
boxplot(geneList)
boxplot(DEG$logFC)

geneList=DEG$logFC
names(geneList)=DEG$ENTREZID
geneList=sort(geneList,decreasing = T)

run_kegg(gene_up,gene_down,pro='comp1')
# 需要多go数据库的3个条目进行3次富集分析，非常耗时。
# run_go(gene_up,gene_down,pro='suiyue')




gene_overlap <- read.table('03-doOverlap-DEG-COX/out-data//03_overlap_logFC_t_1.5_DEG_with_Cox_overlap_genes.txt')[, 1]
gene_overlap



go <- enrichGO(gene_up, OrgDb = "org.Hs.eg.db", ont="all") 
library(ggplot2)
library(stringr)
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free") 



barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave('gene_up_GO_all_barplot.pdf') 
go=DOSE::setReadable(go, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
write.csv(go@result,file = 'gene_up_GO_all_.csv')

go <- enrichGO(gene_down, OrgDb = "org.Hs.eg.db", ont="all") 

barplot(go, split="ONTOLOGY",font.size =10)+ 
  facet_grid(ONTOLOGY~., scale="free") + 
  scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave('gene_down_GO_all_barplot.pdf')
go=DOSE::setReadable(go, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
write.csv(go@result,file = 'gene_down_GO_all_.csv')