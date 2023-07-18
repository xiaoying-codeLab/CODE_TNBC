# Map the network correlation of 9 genes
# Different from PPI network diagram

#seting
rm(list=ls())
options(stringsAsFactors = F)
library(data.table)
library(dplyr)
load('00-folder/out-data//0.expr.Rdata')
gene_overlap <- read.table('03-doOverlap-DEG-COX/out-data//03_overlap_logFC_t_2_DEG_with_Cox_overlap_genes.txt')[, 1]
gene_overlap
expr_tumor = expr_lnc[gene_overlap,]

# Calculate correlation
library(corrplot)
M = cor(t(expr_tumor))
testRes = cor.mtest(t(expr_tumor), conf.level = 0.95)$p
library(tidyverse)
g = pivot_longer(rownames_to_column(as.data.frame(M),var = "from"),
                 cols = 2:(ncol(M)+1),
                 names_to = "to",
                 values_to = "cor")
gp = pivot_longer(rownames_to_column(as.data.frame(testRes)),
                  cols = 2:(ncol(M)+1),
                  names_to = "gene",
                  values_to = "p")
g$p = gp$p
g = g[g$from!=g$to,]
g$group = case_when(g$cor>0.3 & g$p<0.05 ~ "positive",
                    g$cor< -0.3 & g$p<0.05 ~ "negative",
                    T~"not" )
head(g)

# Drawing
library(igraph)
network =  graph_from_data_frame(d=g[g$group!="not",c(1,2,3,5)], directed=F) 
my_color = c("#2874C5","#f87669")[as.numeric(as.factor(E(network)$group))]
par(bg="white", mar=c(0,0,0,0))
pdf("06-cor/cor-2.pdf")
plot(network,
     vertex.size=40,
     layout=layout.circle,
     vertex.label.cex=0.7,
     vertex.frame.color="transparent",
     edge.width=abs(E(network)$cor)*10,
     edge.color=my_color,
     edge.curved = 0.2)
dev.off()



# Draw something else
t(expr_tumor) %>% corrr::correlate() %>%
  corrr::network_plot(colours = c("#2874C5", "white", "#f87669"),repel = F,min_cor = .5,
                      legend = TRUE,
                      curved = TRUE)

ggsave("06-cor/10lnccor4.pdf")

data = corrr::correlate(t(expr_tumor))
write.csv(data,file = "06-cor/cor.data.csv")
data2 = read.csv("06-cor/cor.data2.csv")
data2 = data2[,-1]
data2 = data2[,-2]
data2 = data2[-1,]

# And finally, 8 genes:
corrr::network_plot(data2,colours = c("#2874C5", "white", "#f87669"),repel = F,min_cor = .6,
                    legend = TRUE,
                    curved = TRUE)
ggsave("06-cor/cor-8gene.pdf")

library(ggplot2)
