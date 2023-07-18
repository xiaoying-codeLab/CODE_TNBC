
# 设置
rm(list=ls())
options(stringsAsFactors = F)
library(stringr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(ggsci)
library(ggpubr)
library(patchwork)


load("11_cor_fuction_lnc_mRNA/expr_all2.Rdata")
library(corrplot)
M = cor(t(expr_all2))
testRes = cor.mtest(t(expr_all2), conf.level = 0.95)$p
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
library(igraph)
network =  graph_from_data_frame(d=g[g$group!="not",c(1,2,3,5)], directed=F) 
my_color = c("#2874C5","#f87669")[as.numeric(as.factor(E(network)$group))]
par(bg="white", mar=c(0,0,0,0))
pdf("11_cor_fuction_lnc_mRNA//cor-2.pdf")
plot(network,
     vertex.size=10,
     layout=layout.circle,
     vertex.label.cex=0.2,
     vertex.frame.color="transparent",
     edge.width=abs(E(network)$cor)*10,
     edge.color=my_color,
     edge.curved = 0.2)
dev.off()
t(expr_all2) %>% corrr::correlate() %>%
  corrr::network_plot(colours = c("#2874C5", "white", "#f87669"),repel = F,min_cor = .8,
                      legend = TRUE,
                      curved = TRUE)

ggsave("06-cor/cor.pdf")

data = corrr::correlate(t(expr_all2))


data2 = data[1:30,1:31]

corrr::network_plot(data2,colours = c("#2874C5", "white", "#f87669"),repel = F,min_cor = .4,
                    legend = TRUE,
                    curved = TRUE)
ggsave("11_cor_fuction_lnc_mRNA//cor-1.pdf")
