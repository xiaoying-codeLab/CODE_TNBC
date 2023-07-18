# gase
# 
rm(list = ls())
options(stringsAsFactors = F)
library(clusterProfiler)
library(msigdbr)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(plyr)
library(dplyr)


## 
merge_result2 <- function(enrichResultList, output = "compareClusterResult") {
  if ( !is(enrichResultList, "list")) {
    stop("input should be a name list...")
  }
  if ( is.null(names(enrichResultList))) {
    stop("input should be a name list...")
  }
  x <- lapply(enrichResultList, as.data.frame)
  names(x) <- names(enrichResultList)
  y <- ldply(x, "rbind")   
  
  
  if (output == "compareClusterResult") {
    y <- plyr::rename(y, c(.id="Cluster"))
    y$Cluster = factor(y$Cluster, levels=names(enrichResultList))
    return(new("compareClusterResult",
               compareClusterResult = y))
  } 
  
  y <- plyr::rename(y, c(.id="Category"))
  if (output == "enrichResult") {
    return(new("enrichResult",
               result = y))        
  }
  
  if (output == "gseaResult") {
    return(new("gseaResult",
               result = y))        
  }   
}

keep_category <- function(em_ORA, n) {
  table_em <- as.numeric(table(em_ORA$Category))
  start <- rep(0, length(table_em) - 1)
  for(i in seq_len(length(table_em) - 1)) {
    start[i] <- sum(table_em[seq_len(i)])   
  }
  showCategorys <- sapply(table_em, function(x) min(n, x))
  start <- c(0, start) + 1
  end <- start + showCategorys - 1
  keep <- NULL
  for(i in seq_len(length(start))) {
    keep <- c(keep, c(start[i] : end[i]))
  } 
  return(keep)
}


enrich_filter <- function(em_result, showCategory) {
  keep <- keep_category(em_result, showCategory)
  em_result <- em_result[keep, ]
  if ("NES" %in% colnames(em_result))
    em_result$Count <- em_result$core_enrichment %>% 
    strsplit(split = "/")  %>%  
    vapply(length, FUN.VALUE = 1)
  return(em_result)
}

## 
em_plot <- function(em_1 = NULL, em_2 = NULL, showCategory = 2, fill = "p.adjust", hjust = 1) {
  
  fill <- match.arg(fill, c("Category", "p.adjust", "log10_p.adjust"))
  result1 <- enrich_filter(em_1, showCategory)     
  if (is.null(em_2)) { 
    result <- result1 
  } else {
    result2 <- enrich_filter(em_2, showCategory) 
    result2$Count <- -result2$Count          
    result <- rbind(result1, result2)
  }
  result$Category <- gsub("\n.*", "", result$Category)     
  result$log10_p.adjust <- log10(result$p.adjust)
  
  data_plot <- result[, c("ID", "Category", "p.adjust", "log10_p.adjust", "Count")]
  data_plot2 <- data_plot
  data_plot2$ID <- factor(data_plot2$ID, levels = unique(data_plot2$ID))
  data_plot2 <- plyr::rename(data_plot2, c("Count" = "gene_number"))
  h_just <- ifelse(data_plot2$gene_number < 0, -hjust, hjust)
  ggplot(data_plot2, aes_string(x = "gene_number", y = "ID", fill = fill)) + 
    geom_col() +       
    geom_text(aes_(x =~ gene_number + h_just, label =~ abs(gene_number)), 
              color="black") + 
    scale_x_continuous(label = abs,
                       expand = expansion(mult = c(.01, .01))) + #两侧留空
    theme_classic() + 
    ylab("") + 
    theme(axis.title.x = element_text(size = 15)) +     
    facet_grid(Category ~ ., scales="free", space="free") 
}

# 
msigdbr_species()
gmt <- msigdbr(species = "Homo sapiens") 
gmt2 <- gmt%>%
  dplyr::select(gs_name, entrez_gene)
gmts <- split(gmt2, gmt$gs_cat)

# 
load("12_GSEA/out_data/limma_deg.Rdata")
gsym.fc <- deg
gsym.fc$SYMBOL = rownames(gsym.fc)
gsym.fc[1:3,]
gsym.id <- bitr(gsym.fc$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
#
gsym.fc.id <- merge(gsym.fc, gsym.id, by="SYMBOL", all=F)
#
gsym.fc.id.sorted <- gsym.fc.id[order(gsym.fc.id$logFC, decreasing = T),]
#
id.fc <- gsym.fc.id.sorted$logFC
names(id.fc) <- gsym.fc.id.sorted$ENTREZID


# 
gsea_func <- function(x, genelist, readable = FALSE) {
  GSEA_result <- GSEA(genelist, TERM2GENE = x, eps = 0)
  if (readable & nrow(GSEA_result) > 0) 
    GSEA_result <- setReadable(GSEA_result, 'org.Hs.eg.db', #物种
                               'ENTREZID') #转回gene symbol
  return(GSEA_result)
}
em <- setNames(lapply(gmts, gsea_func, id.fc, readable = TRUE), names(gmts)) %>% 
  merge_result2(output = "gseaResult")

# 
write.csv(em@result, "12_GSEA/out_data/output_GSEA.csv", quote = F)

# 
em_plot(em, showCategory = 5, fill = "log10_p.adjust", hjust = 5)
em_plot(em, showCategory = 5, fill = "Category", hjust = 5) + scale_fill_brewer(palette="Set1")
em_GSEA1 <- em_GSEA2 <- em
res <- em@result
em_GSEA1@result <- res[which(res$NES > 0), ] # 
em_GSEA2@result <- res[which(res$NES < 0), ] # 
em_plot(em_GSEA1, em_GSEA2, showCategory = 2, fill = "log10_p.adjust", hjust = 10)
em_plot(em_GSEA1, em_GSEA2, showCategory = 2, 
        fill = "Category", hjust = 12) +
  scale_fill_brewer(palette="Paired")
ggsave("12_GSEA/out_polt/batchGSEA_sep2.pdf", width = 12, height = 12)

