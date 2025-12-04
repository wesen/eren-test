#!/usr/bin/env Rscript
# Generate test datasets and create ggplot2 visualizations for gene expression comparison
# Based on reference: 01-ggplot2-for-comparing-gene-expression-datasets.md
# Uses library functions from lib/

library(ggplot2)
library(ggrepel)

# Source library functions
source("lib/gene_expr_data.R")
source("lib/gene_expr_plots.R")

# Set seed for reproducibility
set.seed(42)

# Create output directory for plots
output_dir <- "various/plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# Generate Test Datasets
# ============================================================================
n_genes <- 1000

dataset1 <- generate_gene_dataset(n_genes = n_genes, seed = 42)
dataset2 <- generate_correlated_dataset(dataset1, correlation = 0.7, seed = 42)
merged_data <- merge_datasets(dataset1, dataset2)

# Print summary
print_dataset_summary(dataset1, dataset2, merged_data)

# ============================================================================
# Generate All Plots
# ============================================================================

# Plot 1: Basic Scatter Plot
p1 <- plot_scatter_basic(merged_data)
ggsave(file.path(output_dir, "01_scatter_basic.png"), p1, width = 10, height = 8, dpi = 300)
cat("Saved: 01_scatter_basic.png\n")

# Plot 2: Enhanced Scatter Plot with Gene Labels
p2 <- plot_scatter_with_labels(merged_data, fc_threshold = 2, max_labels = 20)
ggsave(file.path(output_dir, "02_scatter_with_labels.png"), p2, width = 10, height = 8, dpi = 300)
cat("Saved: 02_scatter_with_labels.png\n")

# Plot 3: MA Plot (Dataset 1)
p3 <- plot_ma(dataset1, title = "MA Plot - Dataset 1")
ggsave(file.path(output_dir, "03_ma_plot_dataset1.png"), p3, width = 10, height = 8, dpi = 300)
cat("Saved: 03_ma_plot_dataset1.png\n")

# Plot 4: MA Plot (Dataset 2)
p4 <- plot_ma(dataset2, title = "MA Plot - Dataset 2")
ggsave(file.path(output_dir, "04_ma_plot_dataset2.png"), p4, width = 10, height = 8, dpi = 300)
cat("Saved: 04_ma_plot_dataset2.png\n")

# Plot 5: Volcano Plot (Dataset 1)
p5 <- plot_volcano(dataset1, fc_threshold = 1, title = "Volcano Plot - Dataset 1")
ggsave(file.path(output_dir, "05_volcano_plot_dataset1.png"), p5, width = 10, height = 8, dpi = 300)
cat("Saved: 05_volcano_plot_dataset1.png\n")

# Plot 6: Volcano Plot (Dataset 2)
p6 <- plot_volcano(dataset2, fc_threshold = 1, title = "Volcano Plot - Dataset 2")
ggsave(file.path(output_dir, "06_volcano_plot_dataset2.png"), p6, width = 10, height = 8, dpi = 300)
cat("Saved: 06_volcano_plot_dataset2.png\n")

# Plot 7: Multi-Category Significance Scatter Plot
p7 <- plot_scatter_multi_category(merged_data)
ggsave(file.path(output_dir, "07_scatter_multi_category.png"), p7, width = 10, height = 8, dpi = 300)
cat("Saved: 07_scatter_multi_category.png\n")

# Plot 8: ggpubr MA Plot (if available)
if (requireNamespace("ggpubr", quietly = TRUE)) {
  library(ggpubr)
  
  # Prepare data for ggpubr (needs 'name' column)
  dataset1_ggpubr <- dataset1
  dataset1_ggpubr$name <- dataset1_ggpubr$gene_id
  
  p8 <- ggmaplot(
    dataset1_ggpubr,
    main = expression("Condition A" %->% "Condition B"),
    fdr = 0.05,
    fc = 2,
    size = 0.4,
    palette = c("#B31B21", "#1465AC", "darkgray"),
    genenames = as.vector(dataset1_ggpubr$name),
    legend = "top",
    top = 20,
    font.label = c("bold", 11),
    font.legend = "bold",
    font.main = "bold",
    ggtheme = theme_minimal()
  )
  
  ggsave(file.path(output_dir, "08_ggpubr_ma_plot.png"), p8, width = 10, height = 8, dpi = 300)
  cat("Saved: 08_ggpubr_ma_plot.png\n")
} else {
  cat("Skipping ggpubr plot (package not installed)\n")
  cat("Install with: install.packages('ggpubr')\n")
}

# ============================================================================
# Save test datasets as CSV files
# ============================================================================
write.csv(dataset1, file.path(output_dir, "dataset1.csv"), row.names = FALSE)
write.csv(dataset2, file.path(output_dir, "dataset2.csv"), row.names = FALSE)
write.csv(merged_data, file.path(output_dir, "merged_data.csv"), row.names = FALSE)

cat("\nAll plots saved to:", output_dir, "\n")
cat("Test datasets saved as CSV files\n")
cat("Done!\n")
