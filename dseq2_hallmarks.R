#hgsoc from tcga 
#match clinical data to have 20 vs. 20 patients (balanced design for diff exp)
library(tidyr)
library(DESeq2)
library(openxlsx)
library(GEOquery)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(data.table)
library(biomaRt)
library(stringr)
library(edgeR)
library(fgsea)
library(tibble)
library(dplyr)
library(ggpubr)
library(reshape2)
library(org.Hs.eg.db)
library(pheatmap)

pathways.hallmark <- gmtPathways("/home/marchenko/h.all.v2023.1.Hs.symbols.gmt.txt")
setwd("/home/marchenko/fuer_leon/")
# Functions ----
# Function to convert Ensembl Gene IDs to Gene Symbols
IDtoSymbol <- function(expdata) {
  # Set up the connection to Ensembl
  gene_symbols <- AnnotationDbi::select(org.Hs.eg.db, 
                                        keys = as.character(expdata$GeneID),
                                        columns = "SYMBOL", 
                                        keytype = "ENTREZID")
  my_data_with_symbols <- merge(expdata, gene_symbols, by.x = "GeneID", by.y = "ENTREZID", all.x = TRUE)
  
  return(my_data_with_symbols)
}
# read data from the Mao et all paper, which was downloaded from GEO ---- 
info <- read.delim("/home/marchenko/fuer_leon/sample.csv",sep = ",")
raw_counts <- read.delim("/home/marchenko/fuer_leon/GSE153509_raw_counts_GRCh38.p13_NCBI.tsv",sep = "\t")

raw_counts <- IDtoSymbol(expdata = raw_counts) 
raw_counts$GeneID <- NULL
dt <- as.data.table(raw_counts)
dt <- dt[!is.na(SYMBOL)]  # Remove rows with NA SYMBOL
dt <- dt[, lapply(.SD, mean, na.rm = TRUE), by = SYMBOL]  # Summarize by SYMBOL
raw_counts <- as.data.frame(dt) %>% column_to_rownames("SYMBOL")

#prepare matdata
sample_data <- info %>% mutate(Title = gsub("\\.\\.S\\d+$", "", Title)) %>%
  mutate(
    cell_line = case_when(
      str_detect(Title, "^FGFR2_") ~ str_extract(Title, "^FGFR2_[^_]+"),
      str_detect(Title, "^FGF3") ~ "FGF3",
      TRUE ~ str_extract(Title, "^[^_]+")
    ),
    treatment = case_when(
      str_detect(Title, "^FGFR2_") ~ str_replace(Title, "^FGFR2_[^_]+_", ""),
      str_detect(Title, "^FGF3") ~ str_replace(Title, "^FGF3_", ""),
      TRUE ~ str_replace(Title, "^[^_]+_", "")
    )
  )

rownames(sample_data) <- sample_data$Accession
sample_data$cell_line <- factor(sample_data$cell_line)
sample_data$treatment <- factor(sample_data$treatment)
sample_data <- sample_data %>% dplyr::select(cell_line, treatment)

# If you want to set a reference level for the treatment factor ----
sample_data$treatment <- relevel(sample_data$treatment, ref = "DMSO")
metadata_DMSO <- sample_data[sample_data$treatment == "DMSO", ]
counts_DMSO <- raw_counts[, colnames(raw_counts) %in% rownames(metadata_DMSO)]
metadata_DMSO <- metadata_DMSO[match(colnames(counts_DMSO), rownames(metadata_DMSO)), ]
all(colnames(counts_DMSO) == rownames(metadata_DMSO)) 
counts_DMSO <- round(counts_DMSO)
#perform DSEQ2
dds <- DESeqDataSetFromMatrix(countData = counts_DMSO,
                              colData = metadata_DMSO,
                              design = ~ cell_line)
dds <- DESeq(dds)


#for a cell line 
# Define the function
process_cell_line_data <- function(cell_line_name, sample_data, raw_counts) {
  
  # Subset the metadata for the specific cell line
  metadata_cell_line <- sample_data[sample_data$cell_line == cell_line_name, ]
  
  # Subset the raw counts for the specific cell line
  counts_cell_line <- raw_counts[, colnames(raw_counts) %in% rownames(metadata_cell_line)]
  
  # Match column names of counts and row names of metadata
  metadata_cell_line <- metadata_cell_line[match(colnames(counts_cell_line), rownames(metadata_cell_line)), ]
  
  # Check if the column names match
  if (!all(colnames(counts_cell_line) == rownames(metadata_cell_line))) {
    stop("Column names of counts do not match row names of metadata")
  }
  
  # Round the counts data
  counts_cell_line <- round(counts_cell_line)
  
  # Perform DESeq2 analysis
  dds <- DESeqDataSetFromMatrix(countData = counts_cell_line,
                                colData = metadata_cell_line,
                                design = ~ treatment)
  dds <- DESeq(dds)
  
  return(dds)
}

# Example usage: 
# To process for the get DEGs for each cell line
dds_FGF3 <- process_cell_line_data("FGF3", sample_data, raw_counts)
dds_FGFR1 <- process_cell_line_data("FGFR1", sample_data, raw_counts)
dds_FGFR2_WT <- process_cell_line_data("FGFR2_WT", sample_data, raw_counts)
dds_FGFR2_M538I <- process_cell_line_data("FGFR2_M538I", sample_data, raw_counts)
dds_FGFR2_N550K <- process_cell_line_data("FGFR2_N550K", sample_data, raw_counts)
dds_GFP <- process_cell_line_data("GFP", sample_data, raw_counts)
dds_FGFR2_K660N <- process_cell_line_data("FGFR2_K660N", sample_data, raw_counts)
dds_parental <- process_cell_line_data("parental", sample_data, raw_counts)


# Extract DEGs for each comparison (GFP vs parental, FGFR1 vs parental)
cell_lines <- levels(dds$cell_line)

# Remove "parental" from the list (as it is the reference)
cell_lines <- cell_lines[cell_lines != "parental"]

# Create an empty list to store results
results_list <- list()

# Loop through each cell line and extract results
for (line in cell_lines) {
  res <- results(dds, contrast = c("cell_line", line, "parental"))
  results_list[[line]] <- res
}

for (line in cell_lines) {
  write.csv(as.data.frame(results_list[[line]]),
            file = paste0("DEGs_", line, "_vs_parental.csv"))
}


extract_DEGs_for_treatments <- function(dds, treatments_column = "treatment", cell_line_name = "cell_line") {
  # Extract unique treatments (excluding "DMSO" as the reference)
  treatments <- levels(dds[[treatments_column]])
  treatments <- treatments[treatments != "DMSO"]
  
  # Create an empty list to store results
  results_list <- list()
  
  # Loop through each treatment and extract results
  for (treatment in treatments) {
    res <- results(dds, contrast = c(treatments_column, treatment, "DMSO"))
    results_list[[treatment]] <- res
  }
  
  # Save the results to CSV files with the cell line name included in the file name
  for (treatment in treatments) {
    write.csv(as.data.frame(results_list[[treatment]]),
              file = paste0("DEGs_", cell_line_name, "_", treatment, "_vs_DMSO.csv"))
  }
  
  # Return the results_list
  return(results_list)
}

DEGs_FGF3 <- extract_DEGs_for_treatments(dds_FGF3, "treatment",cell_line_name = "FGF3")
DEGs_FGFR1 <- extract_DEGs_for_treatments(dds_FGFR1, "treatment", cell_line_name = "FGFR1")
DEGs_FGFR2_WT <- extract_DEGs_for_treatments(dds_FGFR2_WT, "treatment", cell_line_name = "FGFR2_WT")
DEGs_FGFR2_M538I <- extract_DEGs_for_treatments(dds_FGFR2_M538I, "treatment", cell_line_name = "FGFR2_M538I")
DEGs_FGFR2_N550K <- extract_DEGs_for_treatments(dds_FGFR2_N550K, "treatment", cell_line_name = "FGFR2_N550K")
DEGs_GFP <- extract_DEGs_for_treatments(dds_GFP, "treatment",cell_line_name = "DEGs_GFP")
DEGs_FGFR2_K660N <- extract_DEGs_for_treatments(dds_FGFR2_K660N, "treatment", cell_line_name = "FGFR2_K660N")
DEGs_parental <- extract_DEGs_for_treatments(dds_parental, "treatment",cell_line_name = "parental")


#perform hallmarks 

run_GSEA <- function(data, pathways) {
  
  for_gsea <- data %>% as.data.frame() %>% rownames_to_column()
  res2 <- for_gsea %>% 
    dplyr::select(rowname, log2FoldChange) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(rowname) %>% 
    summarize(stat = mean(log2FoldChange)) %>% 
    ungroup() %>%
    arrange(desc(stat)) 
  
  ranks <- deframe(res2)
  
  fgseaRes <- fgsea(pathways = pathways, stats = ranks,nPermSimple = 10000)
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>% dplyr::filter(.,padj <= 0.05)
  
  fgsea_table <- fgseaResTidy %>% 
    #dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
    arrange(padj)
  
  fgsea_plot <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = padj < 0.05), width = 0.7) +  # Adjust the width of the bars
    coord_flip() +
    labs(x = NULL, y = "Normalized Enrichment Score",  # Remove x-axis label
         title = "Hallmark Pathways NES from GSEA") +
    scale_fill_manual(values = c("TRUE" = "#FF5733", "FALSE" = "#3C8DAD"),  # Custom color palette
                      guide = guide_legend(title = "adjusted p-val < 0.05")) +  # Custom legend title
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold"),  # Increase title size and make it bold
          axis.text.y = element_text(size = 12),  # Increase y-axis label size
          axis.title.x = element_blank(),  # Remove x-axis title
          panel.grid = element_blank())  # Remove grid lines
  
  return(list( table = fgsea_table, plot = fgsea_plot))
}

#hallmarks
#@dmso 
# Run GSEA for each condition and assign the results to variables with names based on the condition
gsea_FGF3  <- run_GSEA(results_list$FGF3, pathways.hallmark)
gsea_FGFR1 <- run_GSEA(results_list$FGFR1, pathways.hallmark)
gsea_FGFR2_K660N <- run_GSEA(results_list$FGFR2_K660N, pathways.hallmark)
gsea_FGFR2_N550K <- run_GSEA(results_list$FGFR2_N550K, pathways.hallmark)
gsea_FGFR2_M538I <- run_GSEA(results_list$FGFR2_M538I, pathways.hallmark)
gsea_FGFR2_WT <- run_GSEA(results_list$FGFR2_WT, pathways.hallmark)
gsea_GFP <- run_GSEA(results_list$GFP, pathways.hallmark)
gsea_parental <- run_GSEA(DEGs, pathways.hallmark)

# Add the source column to each table
gsea_FGF3$table$source  <- "gsea_FGF3"
gsea_FGFR1$table$source  <- "gsea_FGFR1"
gsea_FGFR2_K660N$table$source  <- "gsea_FGFR2_K660N"
gsea_FGFR2_N550K$table$source <- "gsea_FGFR2_N550K"
gsea_FGFR2_M538I$table$source  <- "gsea_FGFR2_M538I"
gsea_FGFR2_WT$table$source  <- "gsea_FGFR2_WT"
gsea_GFP$table$source   <- "gsea_GFP"

# Combine the data frames into one list
df_list <- bind_rows(
  gsea_FGF3$table,
  gsea_FGFR1$table,
  gsea_FGFR2_K660N$table,
  gsea_FGFR2_N550K$table,
  gsea_FGFR2_M538I$table,
  gsea_FGFR2_WT$table,
  gsea_GFP$table
)


generate_balloon_plot <- function(DEGs_list, pathways, output_csv = "gsea_results.csv") {
  # Input: 
  #   DEGs_list: List of DEGs for different treatments or cell lines
  #   pathways: Pathways object (e.g., pathways.hallmark) for GSEA
  #   output_csv: Path to save the combined GSEA results as CSV
  
  # Create an empty list to store GSEA results
  gsea_results <- list()
  
  # Extract treatments or conditions
  treatments <- names(DEGs_list)
  
  # Loop through each entry in the DEGs list
  for (line in treatments) {
    # Run GSEA
    gsea_result <- run_GSEA(DEGs_list[[line]], pathways)
    
    # Check if the table exists in gsea_result and add a source column
    if (!is.null(gsea_result$table)) {
      gsea_result$table$source <- paste0("gsea_", line)
      gsea_results[[line]] <- gsea_result$table
    } else {
      warning(paste("No GSEA table for", line, "- skipping."))
    }
  }
  
  # Combine all GSEA tables into one data frame
  if (length(gsea_results) > 0) {
    df_list <- bind_rows(gsea_results)
    df_list_save <- df_list %>% dplyr::select(pathway, source, NES, padj)
  } else {
    stop("No valid GSEA results found.")
  }
  
  # Check for padj column and calculate -log10(padj)
  if ("padj" %in% colnames(df_list)) {
    
    df_list$log_padj <- -log10(df_list$padj)
  } else {
    stop("Missing 'padj' column in GSEA results.")
  }
  
  # Write the combined results to a CSV file
#  write.csv2(df_list, output_csv, row.names = FALSE)
  
  # Create the balloon plot
  balloon_plot <- ggplot(df_list, aes(x = source, y = pathway, size = log_padj, color = NES)) +
    geom_point() +
    scale_size_continuous(range = c(3, 12)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(
      x = "Source",
      y = "Pathway",
      size = "-log10(padj)",
      color = "NES"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Return the plot and the combined data frame
  list(plot = balloon_plot, gsea_results = df_list_save)
}


# Define a list of DEGs for each cell line
cell_lines_DEGs <- list(
  FGF3 = DEGs_FGF3,
  FGFR1 = DEGs_FGFR1,
  FGFR2_K660N = DEGs_FGFR2_K660N,
  FGFR2_N550K = DEGs_FGFR2_N550K,
  FGFR2_M538I = DEGs_FGFR2_M538I,
  FGFR2_WT = DEGs_FGFR2_WT,
  GFP = DEGs_GFP,
  parental = DEGs_parental
)

# Initialize an empty list to store results
hallmark_results <- list()

for (line in names(cell_lines_DEGs)) {
  # Define output file names for the CSV and PNG
  output_png <- paste0(line, "_balloon_plot.png")
  output_csv <- paste0("Hallmarks", "_", line, ".csv")
  # Run the generate_balloon_plot function
  result <- generate_balloon_plot(cell_lines_DEGs[[line]], pathways.hallmark, output_csv)
  
  # Store the result (plot and data) in the results list
  hallmark_results[[line]] <- result
  
  # Save the plot as PNG
  ggsave(output_png, result$plot, width = 10, height = 8, dpi = 300)
  write.csv2(result$df_list_save, output_csv, row.names = FALSE)
}

