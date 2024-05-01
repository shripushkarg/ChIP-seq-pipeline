library(readr)
library(dplyr)
library(ggplot2)

deseq_results <- read.table("GSE75070_MCF7_shRUNX1_shNS_RNAseq_log2_foldchange.txt", header = TRUE)

annotated_peaks <- read.table("annotated_peaks.txt", header = TRUE, sep = "\t")

filtered_results <- deseq_results %>%
  filter(padj < 0.01, abs(log2FoldChange) > 1)

num_DE_genes <- nrow(filtered_results)

# Print number of DE genes
print(paste("Number of DE genes:", num_DE_genes))

duplicate_genes <- annotated_peaks %>% 
  group_by(`Gene.Name`) %>% 
  filter(n() > 1) %>% 
  pull(`Gene.Name`)

annotated_peaks_unique <- annotated_peaks %>% 
  filter(!`Gene.Name` %in% duplicate_genes)

# Merge DE genes with annotated peak file based on gene ID
merged_filtered_data <- merge(filtered_results, annotated_peaks_unique, by.x = "genename", by.y = "Gene.Name", all.x = TRUE)

calculate_counts_by_distance <- function(data, distance_threshold, total_upregulated, total_downregulated) {
  upregulated_data <- data %>%
    filter(log2FoldChange > 0)
  
  downregulated_data <- data %>%
    filter(log2FoldChange < 0)
  
  bound_upregulated <- upregulated_data %>%
    filter(abs(`Distance.to.TSS`) <= distance_threshold)
  
  bound_downregulated <- downregulated_data %>%
    filter(abs(`Distance.to.TSS`) <= distance_threshold)
  
  count_upregulated_bound <- nrow(bound_upregulated)
  count_downregulated_bound <- nrow(bound_downregulated)
  
  count_upregulated_unbound <- total_upregulated - count_upregulated_bound
  count_downregulated_unbound <- total_downregulated - count_downregulated_bound
  
  list(
    upregulated = list(
      bound = count_upregulated_bound,
      unbound = count_upregulated_unbound
    ),
    downregulated = list(
      bound = count_downregulated_bound,
      unbound = count_downregulated_unbound
    )
  )
}

# Calculate totals for each group of data (upregulated/ downregulated)
total_upregulated <- sum(merged_filtered_data$log2FoldChange > 0)
total_downregulated <- sum(merged_filtered_data$log2FoldChange < 0)

# Counts for 5 kb
counts_5kb <- calculate_counts_by_distance(merged_filtered_data, 5000, total_upregulated, total_downregulated)

# Counts for 20 kb
counts_20kb <- calculate_counts_by_distance(merged_filtered_data, 20000, total_upregulated, total_downregulated)

# Counts for 100 kb
counts_100kb <- calculate_counts_by_distance(merged_filtered_data, 100000, total_upregulated, total_downregulated)

# Function to generate data frame for plotting
generate_kb_data <- function(count_upregulated_bound, count_downregulated_bound, count_upregulated_unbound, count_downregulated_unbound, total_upregulated, total_downregulated) {
  data.frame(
    Group = factor(c("Upregulated", "Downregulated", "Upregulated", "Downregulated")),
    State = factor(c("RUNX1 bound", "RUNX1 bound", "Not bound", "Not bound")),
    Percent = c(count_upregulated_bound / total_upregulated, 
                count_downregulated_bound / total_downregulated,
                count_upregulated_unbound / total_upregulated, 
                count_downregulated_unbound / total_downregulated),
    Counts = c(count_upregulated_bound, count_downregulated_bound, 
               count_upregulated_unbound, count_downregulated_unbound)
  )
}

# Generate data for +/- 5kb of TSS
kb5_data <- generate_kb_data(counts_5kb$upregulated$bound, counts_5kb$downregulated$bound, 
                             counts_5kb$upregulated$unbound, counts_5kb$downregulated$unbound, 
                             total_upregulated, total_downregulated)

# Generate data for +/- 20kb of TSS
kb20_data <- generate_kb_data(counts_20kb$upregulated$bound, counts_20kb$downregulated$bound, 
                              counts_20kb$upregulated$unbound, counts_20kb$downregulated$unbound, 
                              total_upregulated, total_downregulated)

# Generate data for +/- 100kb of TSS
kb100_data <- generate_kb_data(counts_100kb$upregulated$bound, counts_100kb$downregulated$bound, 
                               counts_100kb$upregulated$unbound, counts_100kb$downregulated$unbound, 
                               total_upregulated, total_downregulated)

# Function to plot stacked bar chart
plot_stacked_barchart <- function(data, title, x_label) {
  ggplot(data, aes(fill = State, y = Percent, x = Group)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(title = title, x = x_label, y = "Percentage of genes") +
    theme_minimal() +
    geom_text(aes(label = Counts), position = position_fill(vjust = 0.5), size = 3, color = "white")
}

# Plotting for +/- 5kb of TSS
plot_stacked_barchart(kb5_data, "RUNX1 peak binding +/− 5kb of transcriptional start site (TSS)", "+/- 5kb of TSS")

# Plotting for +/- 20kb of TSS
plot_stacked_barchart(kb20_data, "RUNX1 peak binding +/− 20kb of transcriptional start site (TSS)", "+/- 20kb of TSS")

# Plotting for +/- 100kb of TSS
plot_stacked_barchart(kb100_data, "RUNX1 peak binding +/− 100kb of transcriptional start site (TSS)", "+/- 100kb of TSS")