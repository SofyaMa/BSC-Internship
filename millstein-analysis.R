library(fgsea)

millstein <- read.csv("/home/marchenko/For_Eliane/101-millstein.csv")
load(file= "/home/marchenko/For_Eliane/DESeq_TCGA_matched_samples.RData")
load(file = "/home/marchenko/For_Eliane/orderd_metaData.RData")
pathways.millstein <- gmtPathways("/home/marchenko/For_Eliane/genesets.gmt")

run_GSEA_bulk <- function(data, pathways) {
  
  for_gsea <- data %>% rownames_to_column()
  res2 <- for_gsea %>% 
    dplyr::select(rowname, logFC) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(rowname) %>% 
    summarize(stat = mean(logFC))
  
  ranks <- deframe(res2)
  
  fgseaRes <- fgsea(pathways = pathways, stats = ranks)
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>% head(50)
  fgseaResTidy2 <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>% tail(50)
  
  fgseaResTidy_final <- rbind(fgseaResTidy,fgseaResTidy2)
  fgsea_plot <- ggplot(fgseaResTidy_final, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = padj < 0.05), width = 0.7) +  # Adjust the width of the bars
    coord_flip() +
    labs(x = NULL, y = "Normalized Enrichment Score",  # Remove x-axis label
         title = "Hallmark Pathways NES from GSEA") +
    scale_fill_manual(values = c("TRUE" = "#FF5733", "FALSE" = "#3C8DAD"),  # Custom color palette
                      guide = guide_legend(title = "adjusted p-val < 0.05")) +  # Custom legend title
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold"),  # Increase title size and make it bold
          axis.text.y = element_text(size = 9),  # Increase y-axis label size
          axis.title.x = element_blank(),  # Remove x-axis title
          panel.grid = element_blank())  # Remove grid lines
  
  return(list( table = fgseaResTidy_final, plot = fgsea_plot, ranks = ranks))
}

res <- results(dds)
res_sorted <- res[order((res$log2FoldChange), decreasing = TRUE), ]
tmp <- readRDS(file="/home/marchenko/For_Eliane/diffexp-LTS-STS.rds")
top.table_for_gsea <- topTable(tmp, sort.by = "logFC",n = 100000) 
top.table_for_gsea <- top.table_for_gsea[order(top.table_for_gsea$logFC, decreasing = TRUE), ]
millstein_genes <- millstein$Gene

hallmarks_tcga <- run_GSEA_bulk(as.data.frame(top.table_for_gsea),pathways.millstein)



# Find the ranks of millstein genes in the sorted DESeq2 results
ranks_tcga <- match(millstein_genes, rownames(res_sorted))
ranks <- match(genes, rownames(top.table_for_gsea))
genes <- c(
  "TAP1", "ZFHX4", "CXCL9", "PTGER3", "GMNN", "CD38", "SNRPA1", "SVIL", 
  "ADH1B", "DUSP1", "TRIM27", "PLK2", "ESD", "OLFML3", "RBMS3", "RASA1", 
  "ENOX1", "PCK2", "PARP4", "MRPS27", "GJB1", "GALNT6", "FEN1", "SAC3D1", 
  "MCM3", "RB1", "APBB2", "IGFBP4", "UCP2"
)
# Create a data frame to show the gene names and their corresponding ranks
ranked_genes <- data.frame(
  Gene = genes,
  Rank_in_DESeq = ranks_tcga
)
plot(ranked_genes$Rank_in_DESeq)

# Display the result
print(ranked_genes)
top_30_percent_cutoff <- ceiling(0.5 * nrow(top.table_for_gsea))

# Count how many of the millstein genes are within the top 30%
top_30_genes <- sum(ranks <= top_30_percent_cutoff, na.rm = TRUE)

#Analysis of single genes in the Charité Cohort 

genes <- c("RB1", "CCNE1", "BRCA1", "BRCA2", "KRAS")
titles <- c("Boxplot for RB1", "Boxplot for CCNE1", "Boxplot for BRCA1", "Boxplot for BRCA2", "Boxplot for KRAS")
ylabs <- c("RB1 Expression", "CCNE1 Expression", "BRCA1 Expression", "BRCA2 Expression", "KRAS Expression")

# Loop through the genes to generate the boxplots
for (i in seq_along(genes)) {
  gene <- genes[i]
  title <- titles[i]
  ylab <- ylabs[i]
  
  # Generate boxplot
  boxplot(log(as.numeric(cleaned[gene, ])) + 1 ~ metaData$Survival,
          main = title,
          xlab = "Survival-Groups",
          ylab = ylab,
          col = "lightblue",
          border = "black")
}

