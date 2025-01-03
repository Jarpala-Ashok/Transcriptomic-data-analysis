# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2"))

# Load the libraries
library(clusterProfiler)
library(org.Mm.eg.db)  # 'Mm' is for Mus musculus
library(enrichplot)
library(ggplot2)

# Assume your data is in a CSV file
data <- read.csv("C:/Users/ashok/Downloads/filtered_logFC_±5.csv")

head(data)
str(data)

# Convert gene symbols to Entrez IDs for Mus musculus
genes <- bitr(data$Gene_ID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Check for unmapped genes
unmapped_genes <- data$Gene_ID[!data$Gene_ID %in% genes$SYMBOL]

# Print unmapped genes to review
print(paste("Unmapped Genes:", paste(unmapped_genes, collapse = ", ")))

# Merge the Entrez IDs with your original data, keeping only successfully mapped genes
data_mapped <- merge(data, genes, by.x = "Gene_ID", by.y = "SYMBOL")

# Proceed with the mapped genes for further analysis
gene_list <- data_mapped$ENTREZID


# Gene list with Entrez IDs for the enrichment analysis
gene_list <- data$ENTREZID

# Perform GO enrichment analysis for Mus musculus
go_bp <- enrichGO(gene          = gene_list,
                  OrgDb         = org.Mm.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)

go_cc <- enrichGO(gene          = gene_list,
                  OrgDb         = org.Mm.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)

go_mf <- enrichGO(gene          = gene_list,
                  OrgDb         = org.Mm.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)

barplot(go_bp, showCategory = 10, title = "Biological Processes") + 
  theme(text = element_text(size = 14), # Adjust overall text size
        axis.title = element_text(size = 16),  # Axis titles
        axis.text = element_text(size = 14),   # Axis text
        plot.title = element_text(size = 18, face = "bold"),  # Plot title
        legend.text = element_text(size = 16),  # Legend text
        legend.title = element_text(size = 14)) # Legend title

dotplot(go_bp, showCategory = 10, title = "Biological Processes") + 
  theme(text = element_text(size = 14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))

barplot(go_cc, showCategory = 10, title = "Cellular Components") + 
  theme(text = element_text(size = 14), # Adjust overall text size
        axis.title = element_text(size = 16),  # Axis titles
        axis.text = element_text(size = 14),   # Axis text
        plot.title = element_text(size = 18, face = "bold"),  # Plot title
        legend.text = element_text(size = 16),  # Legend text
        legend.title = element_text(size = 14)) # Legend title

dotplot(go_cc, showCategory = 10, title = "Cellular Components") + 
  theme(text = element_text(size = 14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))


barplot(go_mf, showCategory = 10, title = "Molecular Functions") + 
  theme(text = element_text(size = 20), # Adjust overall text size
        axis.title = element_text(size = 20),  # Axis titles
        axis.text = element_text(size = 20),   # Axis text
        plot.title = element_text(size = 20, face = "bold"),  # Plot title
        legend.text = element_text(size = 20),  # Legend text
        legend.title = element_text(size = 20)) # Legend title

dotplot(go_mf, showCategory = 10, title = "Molecular Functions") + 
  theme(text = element_text(size = 14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))


#barplot(go_bp, showCategory = 10, title = "Biological Processes")
#barplot(go_cc, showCategory = 10, title = "Cellular Components")
#barplot(go_mf, showCategory = 10, title = "Molecular Functions")

#dotplot(go_bp, showCategory = 10, title = "Biological Processes")
#dotplot(go_cc, showCategory = 10, title = "Cellular Components")

#dotplot(go_bp, showCategory = 10, title = "Biological Processes") + 
 # theme(text = element_text(size = 20), 
  #      axis.title = element_text(size = 20),
   #     axis.text = element_text(size = 10),
    #    plot.title = element_text(size = 20, face = "bold"),
     #   legend.text = element_text(size = 20),
      #  legend.title = element_text(size = 20))
#dotplot(go_mf, showCategory = 10, title = "Molecular Functions")

# Print the first few entries in the gene_list
print(head(gene_list))

}
# Convert gene symbols to Entrez IDs for Mus musculus
genes <- bitr(data$Gene_ID, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Print the mapping result to check
print(head(genes))

# Merge the Entrez IDs with your original data
data_mapped <- merge(data, genes, by.x = "Gene_ID", by.y = "SYMBOL", all.x = TRUE)

# Ensure you have a non-empty list of Entrez IDs
gene_list <- na.omit(data_mapped$ENTREZID)

# Proceed with GO enrichment if gene_list is not empty
#if (length(gene_list) > 0) {
 # go_bp <- enrichGO(gene          = gene_list,
  #                  OrgDb         = org.Mm.eg.db,
   #                 ont           = "BP",
    #                pAdjustMethod = "BH",
     #               pvalueCutoff  = 0.05,
      #              qvalueCutoff  = 0.2,
       #             readable      = TRUE)
  #go_cc <- enrichGO(gene          = gene_list,
   #                 OrgDb         = org.Mm.eg.db,
    #                ont           = "CC",
     #               pAdjustMethod = "BH",
      #              pvalueCutoff  = 0.05,
       #             qvalueCutoff  = 0.2,
       #             readable      = TRUE)
  
  #go_mf <- enrichGO(gene          = gene_list,
                  # OrgDb         = org.Mm.eg.db,
                  # ont           = "MF",
                  # pvalueCutoff  = 0.05,
                  # qvalueCutoff  = 0.2,
                  # readable      = TRUE)
#} else {
 # print("No valid Entrez IDs found for GO enrichment.")
#}


# Compute the pairwise term similarity matrix
go_bp <- pairwise_termsim(go_bp)


#barplot(go_bp, showCategory = 10, title = "Biological Processes")
#barplot(go_cc, showCategory = 10, title = "Cellular Components")
#barplot(go_mf, showCategory = 10, title = "Molecular Functions")

#dotplot(go_bp, showCategory = 10, title = "Biological Processes")
#dotplot(go_cc, showCategory = 10, title = "Cellular Components")
#dotplot(go_mf, showCategory = 10, title = "Molecular Functions")





# Compute the pairwise term similarity matrix
go_bp <- pairwise_termsim(go_bp)


# Create and save an emapplot with larger dimensions
emapplot(go_bp) + 
  theme(text = element_text(size = 28),
        axis.title = element_text(size = 28),
        axis.text = element_text(size = 25),
        plot.title = element_text(size = 28, face = "bold"),
        legend.text = element_text(size = 2),
        legend.title = element_text(size = 20),
        plot.margin = margin(5, 5, 5, 5))  # Add extra space around the plot

# Save the plot with increased size
ggsave("emapplot_large.png", plot = last_plot(), width = 20, height = 20, dpi = 300)


cnetplot(go_bp, categorySize = "pvalue", foldChange = data$LogFC) + 
  theme(text = element_text(size = 14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
plot.margin = margin(10, 10, 10, 10)  # Add extra space around the plot
theme(
  # Customize cluster label colors
  plot.title = element_text(color = "red"),  # Change cluster label color here
  # To customize other elements, you can use the following approach:
  # Assuming you have a way to identify cluster labels (adjust as needed)
  # Use `geom_text()` or similar functions if necessary to modify text elements.
  # Note: Actual adjustment may require manual editing or further plotting functions.
)


# Generate the cnetplot and assign it to a variable
cnet <- cnetplot(go_bp, categorySize = "pvalue", foldChange = data$LogFC) 

# Extract the data layers of the cnetplot
cnet_data <- ggplot_build(cnet)


# Recreate the plot with custom colors
cnet_custom <- cnet +
  geom_text(data = cnet_data$data[[1]], 
            aes(label = name, color = "red"), # Change color of clusters to red
            size = 5, vjust = -1) +
  geom_text(data = cnet_data$data[[3]], 
            aes(label = name, color = "black"), # Keep genes' labels black
            size = 4) +
  theme_minimal() +
  theme(text = element_text(size = 14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.margin = margin(10, 10, 10, 10))





# Save the plot with increased size
ggsave("cnetplot_large.png", plot = last_plot(), width = 20, height = 20, dpi = 300)


# Compute the pairwise term similarity matrix
go_cc <- pairwise_termsim(go_cc)


# Create and save an emapplot with larger dimensions
emapplot(go_cc) + 
  theme(text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 16),
        plot.margin = margin(5, 5, 5, 5))  # Add extra space around the plot

# Save the plot with increased size
ggsave("emapplot_large.png", plot = last_plot(), width = 20, height = 20, dpi = 300)


cnetplot(go_cc, categorySize = "pvalue", foldChange = data$LogFC) + 
  theme(text = element_text(size = 14), 
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16))
plot.margin = margin(10, 10, 10, 10)  # Add extra space around the plot

# Save the plot with increased size
ggsave("cnetplot_large.png", plot = last_plot(), width = 20, height = 20, dpi = 300)



# Compute the pairwise term similarity matrix
go_mf <- pairwise_termsim(go_mf)


# Create and save an emapplot with larger dimensions
emapplot(go_mf) + 
  theme(text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        plot.title = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        plot.margin = margin(5, 5, 5, 5))  # Add extra space around the plot

# Save the plot with increased size
ggsave("emapplot_large.png", plot = last_plot(), width = 20, height = 20, dpi = 300)


cnetplot(go_mf, categorySize = "pvalue", foldChange = data$LogFC) + 
  theme(text = element_text(size = 20), 
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25))
plot.margin = margin(10, 10, 10, 10)  # Add extra space around the plot

# Save the plot with increased size
ggsave("cnetplot_large.png", plot = last_plot(), width = 20, height = 20, dpi = 300)


  
  
  
  
}


