# Create test datasets and graphs
_Exported on 12/4/2025 at 08:38:10 EST from Cursor (2.1.39)_

---

**User**

Using @ttmp/2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/reference/01-ggplot2-for-comparing-gene-expression-datasets.md , create some test datasets and make graphs.

---

**Cursor**


Creating test datasets and generating graphs from the reference. Checking the project structure to decide where to place the script.



Creating an R script to generate test datasets and produce the visualizations from the reference document.


```r
#!/usr/bin/env Rscript
# Generate test datasets and create ggplot2 visualizations for gene expression comparison
# Based on reference: 01-ggplot2-for-comparing-gene-expression-datasets.md

library(ggplot2)
library(ggrepel)

# Set seed for reproducibility
set.seed(42)

# Create output directory for plots
output_dir <- "various/plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# Generate Test Dataset 1
# ============================================================================
n_genes <- 1000

# Generate gene IDs
gene_ids <- paste0("GENE_", sprintf("%04d", 1:n_genes))

# Generate log2 fold changes (with some correlation structure)
logFC_ds1 <- rnorm(n_genes, mean = 0, sd = 1.5)
# Add some genes with larger fold changes
logFC_ds1[sample(1:n_genes, 100)] <- rnorm(100, mean = 0, sd = 3)

# Generate base mean expression (log scale)
baseMean <- 10^rnorm(n_genes, mean = 3, sd = 1.5)
baseMean <- pmax(baseMean, 1)  # Ensure positive values

# Generate p-values (smaller p-values for larger fold changes)
pvalue <- runif(n_genes, min = 0.001, max = 0.5)
# Make p-values smaller for genes with larger absolute fold changes
pvalue <- pvalue * (1 - abs(logFC_ds1) / max(abs(logFC_ds1)) * 0.8)

# Generate adjusted p-values (padj)
padj <- p.adjust(pvalue, method = "BH")

# Create dataset 1
dataset1 <- data.frame(
  gene_id = gene_ids,
  log2FoldChange = logFC_ds1,
  baseMean = baseMean,
  pvalue = pvalue,
  padj = padj,
  FDR = padj  # FDR is same as padj for BH correction
)

# ============================================================================
# Generate Test Dataset 2 (correlated with dataset 1 but with some noise)
# ============================================================================
# Create correlated log fold changes (with noise)
correlation <- 0.7
logFC_ds2 <- correlation * logFC_ds1 + sqrt(1 - correlation^2) * rnorm(n_genes, mean = 0, sd = 1.5)
# Add some dataset-specific differences
logFC_ds2[sample(1:n_genes, 50)] <- rnorm(50, mean = 0, sd = 2)

# Generate p-values for dataset 2
pvalue_ds2 <- runif(n_genes, min = 0.001, max = 0.5)
pvalue_ds2 <- pvalue_ds2 * (1 - abs(logFC_ds2) / max(abs(logFC_ds2)) * 0.8)
padj_ds2 <- p.adjust(pvalue_ds2, method = "BH")

# Create dataset 2
dataset2 <- data.frame(
  gene_id = gene_ids,
  log2FoldChange = logFC_ds2,
  baseMean = baseMean,  # Same base mean for comparison
  pvalue = pvalue_ds2,
  padj = padj_ds2,
  FDR = padj_ds2
)

# ============================================================================
# Create Merged Dataset for Comparison
# ============================================================================
merged_data <- merge(
  dataset1[, c("gene_id", "log2FoldChange", "padj", "FDR")],
  dataset2[, c("gene_id", "log2FoldChange", "padj", "FDR")],
  by = "gene_id",
  suffixes = c("_ds1", "_ds2")
)

# Rename columns for easier plotting
merged_data$logFC_dataset1 <- merged_data$log2FoldChange_ds1
merged_data$logFC_dataset2 <- merged_data$log2FoldChange_ds2
merged_data$FDR <- pmax(merged_data$FDR_ds1, merged_data$FDR_ds2)  # Use max FDR

# Add baseMean from dataset1
merged_data$baseMean <- dataset1$baseMean[match(merged_data$gene_id, dataset1$gene_id)]

cat("Generated", n_genes, "genes\n")
cat("Dataset 1: Significant genes (padj < 0.05):", sum(dataset1$padj < 0.05), "\n")
cat("Dataset 2: Significant genes (padj < 0.05):", sum(dataset2$padj < 0.05), "\n")
cat("Both datasets significant:", sum(merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05), "\n")

# ============================================================================
# Plot 1: Basic Scatter Plot Comparing Two Datasets
# ============================================================================
p1 <- ggplot(merged_data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    title = "Comparison of Log2 Fold Changes Between Two Datasets",
    color = "Significance"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "01_scatter_basic.png"), p1, width = 10, height = 8, dpi = 300)
cat("Saved: 01_scatter_basic.png\n")

# ============================================================================
# Plot 2: Enhanced Scatter Plot with Gene Labels
# ============================================================================
p2 <- ggplot(merged_data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_text_repel(
    data = subset(merged_data, FDR < 0.05 & abs(logFC_dataset1) > 2 & abs(logFC_dataset2) > 2),
    aes(label = gene_id),
    size = 3,
    max.overlaps = 20
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    color = "Significant (FDR < 0.05)"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "02_scatter_with_labels.png"), p2, width = 10, height = 8, dpi = 300)
cat("Saved: 02_scatter_with_labels.png\n")

# ============================================================================
# Plot 3: MA Plot (Dataset 1)
# ============================================================================
p3 <- ggplot(dataset1, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significant (padj < 0.05)",
    title = "MA Plot - Dataset 1"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "03_ma_plot_dataset1.png"), p3, width = 10, height = 8, dpi = 300)
cat("Saved: 03_ma_plot_dataset1.png\n")

# ============================================================================
# Plot 4: MA Plot (Dataset 2)
# ============================================================================
p4 <- ggplot(dataset2, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significant (padj < 0.05)",
    title = "MA Plot - Dataset 2"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "04_ma_plot_dataset2.png"), p4, width = 10, height = 8, dpi = 300)
cat("Saved: 04_ma_plot_dataset2.png\n")

# ============================================================================
# Plot 5: Volcano Plot (Dataset 1)
# ============================================================================
p5 <- ggplot(dataset1, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant (padj < 0.05)",
    title = "Volcano Plot - Dataset 1"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "05_volcano_plot_dataset1.png"), p5, width = 10, height = 8, dpi = 300)
cat("Saved: 05_volcano_plot_dataset1.png\n")

# ============================================================================
# Plot 6: Volcano Plot (Dataset 2)
# ============================================================================
p6 <- ggplot(dataset2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant (padj < 0.05)",
    title = "Volcano Plot - Dataset 2"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "06_volcano_plot_dataset2.png"), p6, width = 10, height = 8, dpi = 300)
cat("Saved: 06_volcano_plot_dataset2.png\n")

# ============================================================================
# Plot 7: Multi-Category Significance Scatter Plot
# ============================================================================
merged_data$sig_status <- ifelse(
  merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05, "Both significant",
  ifelse(merged_data$padj_ds1 < 0.05, "Dataset 1 only",
  ifelse(merged_data$padj_ds2 < 0.05, "Dataset 2 only", "Neither"))
)

p7 <- ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
  geom_point(aes(color = sig_status), alpha = 0.6) +
  scale_color_manual(
    values = c(
      "Both significant" = "red",
      "Dataset 1 only" = "orange",
      "Dataset 2 only" = "blue",
      "Neither" = "gray60"
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Dataset 1 Log2FC",
    y = "Dataset 2 Log2FC",
    color = "Significance",
    title = "Multi-Category Significance Comparison"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "07_scatter_multi_category.png"), p7, width = 10, height = 8, dpi = 300)
cat("Saved: 07_scatter_multi_category.png\n")

# ============================================================================
# Plot 8: ggpubr MA Plot (if available)
# ============================================================================
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

```


Running the script to generate test datasets and create the graphs:



Installing required R packages, then rerunning the script:



Rerunning the script:



The script ran successfully. No significant genes were generated. Updating the script to produce more realistic data with some significant genes:



Updating the p-value generation logic to produce more realistic data with significant genes:


```r
#!/usr/bin/env Rscript
# Generate test datasets and create ggplot2 visualizations for gene expression comparison
# Based on reference: 01-ggplot2-for-comparing-gene-expression-datasets.md

library(ggplot2)
library(ggrepel)

# Set seed for reproducibility
set.seed(42)

# Create output directory for plots
output_dir <- "various/plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# Generate Test Dataset 1
# ============================================================================
n_genes <- 1000

# Generate gene IDs
gene_ids <- paste0("GENE_", sprintf("%04d", 1:n_genes))

# Generate log2 fold changes (with some correlation structure)
logFC_ds1 <- rnorm(n_genes, mean = 0, sd = 1.5)
# Add some genes with larger fold changes
logFC_ds1[sample(1:n_genes, 100)] <- rnorm(100, mean = 0, sd = 3)

# Generate base mean expression (log scale)
baseMean <- 10^rnorm(n_genes, mean = 3, sd = 1.5)
baseMean <- pmax(baseMean, 1)  # Ensure positive values

# Generate p-values (smaller p-values for larger fold changes)
# Create a mix: some genes with very small p-values, others with larger ones
pvalue <- ifelse(
  abs(logFC_ds1) > 1.5,  # Genes with larger fold changes
  runif(n_genes, min = 0.0001, max = 0.01),  # Small p-values
  runif(n_genes, min = 0.01, max = 0.5)  # Larger p-values
)
# Add some randomness
pvalue <- pvalue * runif(n_genes, min = 0.5, max = 1.5)
pvalue <- pmax(pvalue, 1e-10)  # Ensure minimum value

# Generate adjusted p-values (padj)
padj <- p.adjust(pvalue, method = "BH")

# Create dataset 1
dataset1 <- data.frame(
  gene_id = gene_ids,
  log2FoldChange = logFC_ds1,
  baseMean = baseMean,
  pvalue = pvalue,
  padj = padj,
  FDR = padj  # FDR is same as padj for BH correction
)

# ============================================================================
# Generate Test Dataset 2 (correlated with dataset 1 but with some noise)
# ============================================================================
# Create correlated log fold changes (with noise)
correlation <- 0.7
logFC_ds2 <- correlation * logFC_ds1 + sqrt(1 - correlation^2) * rnorm(n_genes, mean = 0, sd = 1.5)
# Add some dataset-specific differences
logFC_ds2[sample(1:n_genes, 50)] <- rnorm(50, mean = 0, sd = 2)

# Generate p-values for dataset 2
pvalue_ds2 <- runif(n_genes, min = 0.001, max = 0.5)
pvalue_ds2 <- pvalue_ds2 * (1 - abs(logFC_ds2) / max(abs(logFC_ds2)) * 0.8)
padj_ds2 <- p.adjust(pvalue_ds2, method = "BH")

# Create dataset 2
dataset2 <- data.frame(
  gene_id = gene_ids,
  log2FoldChange = logFC_ds2,
  baseMean = baseMean,  # Same base mean for comparison
  pvalue = pvalue_ds2,
  padj = padj_ds2,
  FDR = padj_ds2
)

# ============================================================================
# Create Merged Dataset for Comparison
# ============================================================================
merged_data <- merge(
  dataset1[, c("gene_id", "log2FoldChange", "padj", "FDR")],
  dataset2[, c("gene_id", "log2FoldChange", "padj", "FDR")],
  by = "gene_id",
  suffixes = c("_ds1", "_ds2")
)

# Rename columns for easier plotting
merged_data$logFC_dataset1 <- merged_data$log2FoldChange_ds1
merged_data$logFC_dataset2 <- merged_data$log2FoldChange_ds2
merged_data$FDR <- pmax(merged_data$FDR_ds1, merged_data$FDR_ds2)  # Use max FDR

# Add baseMean from dataset1
merged_data$baseMean <- dataset1$baseMean[match(merged_data$gene_id, dataset1$gene_id)]

cat("Generated", n_genes, "genes\n")
cat("Dataset 1: Significant genes (padj < 0.05):", sum(dataset1$padj < 0.05), "\n")
cat("Dataset 2: Significant genes (padj < 0.05):", sum(dataset2$padj < 0.05), "\n")
cat("Both datasets significant:", sum(merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05), "\n")

# ============================================================================
# Plot 1: Basic Scatter Plot Comparing Two Datasets
# ============================================================================
p1 <- ggplot(merged_data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    title = "Comparison of Log2 Fold Changes Between Two Datasets",
    color = "Significance"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "01_scatter_basic.png"), p1, width = 10, height = 8, dpi = 300)
cat("Saved: 01_scatter_basic.png\n")

# ============================================================================
# Plot 2: Enhanced Scatter Plot with Gene Labels
# ============================================================================
p2 <- ggplot(merged_data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_text_repel(
    data = subset(merged_data, FDR < 0.05 & abs(logFC_dataset1) > 2 & abs(logFC_dataset2) > 2),
    aes(label = gene_id),
    size = 3,
    max.overlaps = 20
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    color = "Significant (FDR < 0.05)"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "02_scatter_with_labels.png"), p2, width = 10, height = 8, dpi = 300)
cat("Saved: 02_scatter_with_labels.png\n")

# ============================================================================
# Plot 3: MA Plot (Dataset 1)
# ============================================================================
p3 <- ggplot(dataset1, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significant (padj < 0.05)",
    title = "MA Plot - Dataset 1"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "03_ma_plot_dataset1.png"), p3, width = 10, height = 8, dpi = 300)
cat("Saved: 03_ma_plot_dataset1.png\n")

# ============================================================================
# Plot 4: MA Plot (Dataset 2)
# ============================================================================
p4 <- ggplot(dataset2, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significant (padj < 0.05)",
    title = "MA Plot - Dataset 2"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "04_ma_plot_dataset2.png"), p4, width = 10, height = 8, dpi = 300)
cat("Saved: 04_ma_plot_dataset2.png\n")

# ============================================================================
# Plot 5: Volcano Plot (Dataset 1)
# ============================================================================
p5 <- ggplot(dataset1, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant (padj < 0.05)",
    title = "Volcano Plot - Dataset 1"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "05_volcano_plot_dataset1.png"), p5, width = 10, height = 8, dpi = 300)
cat("Saved: 05_volcano_plot_dataset1.png\n")

# ============================================================================
# Plot 6: Volcano Plot (Dataset 2)
# ============================================================================
p6 <- ggplot(dataset2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant (padj < 0.05)",
    title = "Volcano Plot - Dataset 2"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "06_volcano_plot_dataset2.png"), p6, width = 10, height = 8, dpi = 300)
cat("Saved: 06_volcano_plot_dataset2.png\n")

# ============================================================================
# Plot 7: Multi-Category Significance Scatter Plot
# ============================================================================
merged_data$sig_status <- ifelse(
  merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05, "Both significant",
  ifelse(merged_data$padj_ds1 < 0.05, "Dataset 1 only",
  ifelse(merged_data$padj_ds2 < 0.05, "Dataset 2 only", "Neither"))
)

p7 <- ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
  geom_point(aes(color = sig_status), alpha = 0.6) +
  scale_color_manual(
    values = c(
      "Both significant" = "red",
      "Dataset 1 only" = "orange",
      "Dataset 2 only" = "blue",
      "Neither" = "gray60"
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Dataset 1 Log2FC",
    y = "Dataset 2 Log2FC",
    color = "Significance",
    title = "Multi-Category Significance Comparison"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "07_scatter_multi_category.png"), p7, width = 10, height = 8, dpi = 300)
cat("Saved: 07_scatter_multi_category.png\n")

# ============================================================================
# Plot 8: ggpubr MA Plot (if available)
# ============================================================================
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


```

```r
#!/usr/bin/env Rscript
# Generate test datasets and create ggplot2 visualizations for gene expression comparison
# Based on reference: 01-ggplot2-for-comparing-gene-expression-datasets.md

library(ggplot2)
library(ggrepel)

# Set seed for reproducibility
set.seed(42)

# Create output directory for plots
output_dir <- "various/plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# Generate Test Dataset 1
# ============================================================================
n_genes <- 1000

# Generate gene IDs
gene_ids <- paste0("GENE_", sprintf("%04d", 1:n_genes))

# Generate log2 fold changes (with some correlation structure)
logFC_ds1 <- rnorm(n_genes, mean = 0, sd = 1.5)
# Add some genes with larger fold changes
logFC_ds1[sample(1:n_genes, 100)] <- rnorm(100, mean = 0, sd = 3)

# Generate base mean expression (log scale)
baseMean <- 10^rnorm(n_genes, mean = 3, sd = 1.5)
baseMean <- pmax(baseMean, 1)  # Ensure positive values

# Generate p-values (smaller p-values for larger fold changes)
# Create a mix: some genes with very small p-values, others with larger ones
pvalue <- ifelse(
  abs(logFC_ds1) > 1.5,  # Genes with larger fold changes
  runif(n_genes, min = 0.0001, max = 0.01),  # Small p-values
  runif(n_genes, min = 0.01, max = 0.5)  # Larger p-values
)
# Add some randomness
pvalue <- pvalue * runif(n_genes, min = 0.5, max = 1.5)
pvalue <- pmax(pvalue, 1e-10)  # Ensure minimum value

# Generate adjusted p-values (padj)
padj <- p.adjust(pvalue, method = "BH")

# Create dataset 1
dataset1 <- data.frame(
  gene_id = gene_ids,
  log2FoldChange = logFC_ds1,
  baseMean = baseMean,
  pvalue = pvalue,
  padj = padj,
  FDR = padj  # FDR is same as padj for BH correction
)

# ============================================================================
# Generate Test Dataset 2 (correlated with dataset 1 but with some noise)
# ============================================================================
# Create correlated log fold changes (with noise)
correlation <- 0.7
logFC_ds2 <- correlation * logFC_ds1 + sqrt(1 - correlation^2) * rnorm(n_genes, mean = 0, sd = 1.5)
# Add some dataset-specific differences
logFC_ds2[sample(1:n_genes, 50)] <- rnorm(50, mean = 0, sd = 2)

# Generate p-values for dataset 2
# Create a mix: some genes with very small p-values, others with larger ones
pvalue_ds2 <- ifelse(
  abs(logFC_ds2) > 1.5,  # Genes with larger fold changes
  runif(n_genes, min = 0.0001, max = 0.01),  # Small p-values
  runif(n_genes, min = 0.01, max = 0.5)  # Larger p-values
)
# Add some randomness
pvalue_ds2 <- pvalue_ds2 * runif(n_genes, min = 0.5, max = 1.5)
pvalue_ds2 <- pmax(pvalue_ds2, 1e-10)  # Ensure minimum value
padj_ds2 <- p.adjust(pvalue_ds2, method = "BH")

# Create dataset 2
dataset2 <- data.frame(
  gene_id = gene_ids,
  log2FoldChange = logFC_ds2,
  baseMean = baseMean,  # Same base mean for comparison
  pvalue = pvalue_ds2,
  padj = padj_ds2,
  FDR = padj_ds2
)

# ============================================================================
# Create Merged Dataset for Comparison
# ============================================================================
merged_data <- merge(
  dataset1[, c("gene_id", "log2FoldChange", "padj", "FDR")],
  dataset2[, c("gene_id", "log2FoldChange", "padj", "FDR")],
  by = "gene_id",
  suffixes = c("_ds1", "_ds2")
)

# Rename columns for easier plotting
merged_data$logFC_dataset1 <- merged_data$log2FoldChange_ds1
merged_data$logFC_dataset2 <- merged_data$log2FoldChange_ds2
merged_data$FDR <- pmax(merged_data$FDR_ds1, merged_data$FDR_ds2)  # Use max FDR

# Add baseMean from dataset1
merged_data$baseMean <- dataset1$baseMean[match(merged_data$gene_id, dataset1$gene_id)]

cat("Generated", n_genes, "genes\n")
cat("Dataset 1: Significant genes (padj < 0.05):", sum(dataset1$padj < 0.05), "\n")
cat("Dataset 2: Significant genes (padj < 0.05):", sum(dataset2$padj < 0.05), "\n")
cat("Both datasets significant:", sum(merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05), "\n")

# ============================================================================
# Plot 1: Basic Scatter Plot Comparing Two Datasets
# ============================================================================
p1 <- ggplot(merged_data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    title = "Comparison of Log2 Fold Changes Between Two Datasets",
    color = "Significance"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "01_scatter_basic.png"), p1, width = 10, height = 8, dpi = 300)
cat("Saved: 01_scatter_basic.png\n")

# ============================================================================
# Plot 2: Enhanced Scatter Plot with Gene Labels
# ============================================================================
p2 <- ggplot(merged_data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_text_repel(
    data = subset(merged_data, FDR < 0.05 & abs(logFC_dataset1) > 2 & abs(logFC_dataset2) > 2),
    aes(label = gene_id),
    size = 3,
    max.overlaps = 20
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    color = "Significant (FDR < 0.05)"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "02_scatter_with_labels.png"), p2, width = 10, height = 8, dpi = 300)
cat("Saved: 02_scatter_with_labels.png\n")

# ============================================================================
# Plot 3: MA Plot (Dataset 1)
# ============================================================================
p3 <- ggplot(dataset1, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significant (padj < 0.05)",
    title = "MA Plot - Dataset 1"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "03_ma_plot_dataset1.png"), p3, width = 10, height = 8, dpi = 300)
cat("Saved: 03_ma_plot_dataset1.png\n")

# ============================================================================
# Plot 4: MA Plot (Dataset 2)
# ============================================================================
p4 <- ggplot(dataset2, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significant (padj < 0.05)",
    title = "MA Plot - Dataset 2"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "04_ma_plot_dataset2.png"), p4, width = 10, height = 8, dpi = 300)
cat("Saved: 04_ma_plot_dataset2.png\n")

# ============================================================================
# Plot 5: Volcano Plot (Dataset 1)
# ============================================================================
p5 <- ggplot(dataset1, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant (padj < 0.05)",
    title = "Volcano Plot - Dataset 1"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "05_volcano_plot_dataset1.png"), p5, width = 10, height = 8, dpi = 300)
cat("Saved: 05_volcano_plot_dataset1.png\n")

# ============================================================================
# Plot 6: Volcano Plot (Dataset 2)
# ============================================================================
p6 <- ggplot(dataset2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant (padj < 0.05)",
    title = "Volcano Plot - Dataset 2"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "06_volcano_plot_dataset2.png"), p6, width = 10, height = 8, dpi = 300)
cat("Saved: 06_volcano_plot_dataset2.png\n")

# ============================================================================
# Plot 7: Multi-Category Significance Scatter Plot
# ============================================================================
merged_data$sig_status <- ifelse(
  merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05, "Both significant",
  ifelse(merged_data$padj_ds1 < 0.05, "Dataset 1 only",
  ifelse(merged_data$padj_ds2 < 0.05, "Dataset 2 only", "Neither"))
)

p7 <- ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
  geom_point(aes(color = sig_status), alpha = 0.6) +
  scale_color_manual(
    values = c(
      "Both significant" = "red",
      "Dataset 1 only" = "orange",
      "Dataset 2 only" = "blue",
      "Neither" = "gray60"
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Dataset 1 Log2FC",
    y = "Dataset 2 Log2FC",
    color = "Significance",
    title = "Multi-Category Significance Comparison"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "07_scatter_multi_category.png"), p7, width = 10, height = 8, dpi = 300)
cat("Saved: 07_scatter_multi_category.png\n")

# ============================================================================
# Plot 8: ggpubr MA Plot (if available)
# ============================================================================
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


```


Fixing the p-value generation logic:


```r
#!/usr/bin/env Rscript
# Generate test datasets and create ggplot2 visualizations for gene expression comparison
# Based on reference: 01-ggplot2-for-comparing-gene-expression-datasets.md

library(ggplot2)
library(ggrepel)

# Set seed for reproducibility
set.seed(42)

# Create output directory for plots
output_dir <- "various/plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# Generate Test Dataset 1
# ============================================================================
n_genes <- 1000

# Generate gene IDs
gene_ids <- paste0("GENE_", sprintf("%04d", 1:n_genes))

# Generate log2 fold changes (with some correlation structure)
logFC_ds1 <- rnorm(n_genes, mean = 0, sd = 1.5)
# Add some genes with larger fold changes
logFC_ds1[sample(1:n_genes, 100)] <- rnorm(100, mean = 0, sd = 3)

# Generate base mean expression (log scale)
baseMean <- 10^rnorm(n_genes, mean = 3, sd = 1.5)
baseMean <- pmax(baseMean, 1)  # Ensure positive values

# Generate p-values (smaller p-values for larger fold changes)
# Create a mix: some genes with very small p-values, others with larger ones
large_fc_idx <- abs(logFC_ds1) > 1.5
pvalue <- numeric(n_genes)
pvalue[large_fc_idx] <- runif(sum(large_fc_idx), min = 0.0001, max = 0.01)
pvalue[!large_fc_idx] <- runif(sum(!large_fc_idx), min = 0.01, max = 0.5)
# Add some randomness
pvalue <- pvalue * runif(n_genes, min = 0.5, max = 1.5)
pvalue <- pmax(pvalue, 1e-10)  # Ensure minimum value

# Generate adjusted p-values (padj)
padj <- p.adjust(pvalue, method = "BH")

# Create dataset 1
dataset1 <- data.frame(
  gene_id = gene_ids,
  log2FoldChange = logFC_ds1,
  baseMean = baseMean,
  pvalue = pvalue,
  padj = padj,
  FDR = padj  # FDR is same as padj for BH correction
)

# ============================================================================
# Generate Test Dataset 2 (correlated with dataset 1 but with some noise)
# ============================================================================
# Create correlated log fold changes (with noise)
correlation <- 0.7
logFC_ds2 <- correlation * logFC_ds1 + sqrt(1 - correlation^2) * rnorm(n_genes, mean = 0, sd = 1.5)
# Add some dataset-specific differences
logFC_ds2[sample(1:n_genes, 50)] <- rnorm(50, mean = 0, sd = 2)

# Generate p-values for dataset 2
# Create a mix: some genes with very small p-values, others with larger ones
pvalue_ds2 <- ifelse(
  abs(logFC_ds2) > 1.5,  # Genes with larger fold changes
  runif(n_genes, min = 0.0001, max = 0.01),  # Small p-values
  runif(n_genes, min = 0.01, max = 0.5)  # Larger p-values
)
# Add some randomness
pvalue_ds2 <- pvalue_ds2 * runif(n_genes, min = 0.5, max = 1.5)
pvalue_ds2 <- pmax(pvalue_ds2, 1e-10)  # Ensure minimum value
padj_ds2 <- p.adjust(pvalue_ds2, method = "BH")

# Create dataset 2
dataset2 <- data.frame(
  gene_id = gene_ids,
  log2FoldChange = logFC_ds2,
  baseMean = baseMean,  # Same base mean for comparison
  pvalue = pvalue_ds2,
  padj = padj_ds2,
  FDR = padj_ds2
)

# ============================================================================
# Create Merged Dataset for Comparison
# ============================================================================
merged_data <- merge(
  dataset1[, c("gene_id", "log2FoldChange", "padj", "FDR")],
  dataset2[, c("gene_id", "log2FoldChange", "padj", "FDR")],
  by = "gene_id",
  suffixes = c("_ds1", "_ds2")
)

# Rename columns for easier plotting
merged_data$logFC_dataset1 <- merged_data$log2FoldChange_ds1
merged_data$logFC_dataset2 <- merged_data$log2FoldChange_ds2
merged_data$FDR <- pmax(merged_data$FDR_ds1, merged_data$FDR_ds2)  # Use max FDR

# Add baseMean from dataset1
merged_data$baseMean <- dataset1$baseMean[match(merged_data$gene_id, dataset1$gene_id)]

cat("Generated", n_genes, "genes\n")
cat("Dataset 1: Significant genes (padj < 0.05):", sum(dataset1$padj < 0.05), "\n")
cat("Dataset 2: Significant genes (padj < 0.05):", sum(dataset2$padj < 0.05), "\n")
cat("Both datasets significant:", sum(merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05), "\n")

# ============================================================================
# Plot 1: Basic Scatter Plot Comparing Two Datasets
# ============================================================================
p1 <- ggplot(merged_data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    title = "Comparison of Log2 Fold Changes Between Two Datasets",
    color = "Significance"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "01_scatter_basic.png"), p1, width = 10, height = 8, dpi = 300)
cat("Saved: 01_scatter_basic.png\n")

# ============================================================================
# Plot 2: Enhanced Scatter Plot with Gene Labels
# ============================================================================
p2 <- ggplot(merged_data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_text_repel(
    data = subset(merged_data, FDR < 0.05 & abs(logFC_dataset1) > 2 & abs(logFC_dataset2) > 2),
    aes(label = gene_id),
    size = 3,
    max.overlaps = 20
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    color = "Significant (FDR < 0.05)"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "02_scatter_with_labels.png"), p2, width = 10, height = 8, dpi = 300)
cat("Saved: 02_scatter_with_labels.png\n")

# ============================================================================
# Plot 3: MA Plot (Dataset 1)
# ============================================================================
p3 <- ggplot(dataset1, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significant (padj < 0.05)",
    title = "MA Plot - Dataset 1"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "03_ma_plot_dataset1.png"), p3, width = 10, height = 8, dpi = 300)
cat("Saved: 03_ma_plot_dataset1.png\n")

# ============================================================================
# Plot 4: MA Plot (Dataset 2)
# ============================================================================
p4 <- ggplot(dataset2, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significant (padj < 0.05)",
    title = "MA Plot - Dataset 2"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "04_ma_plot_dataset2.png"), p4, width = 10, height = 8, dpi = 300)
cat("Saved: 04_ma_plot_dataset2.png\n")

# ============================================================================
# Plot 5: Volcano Plot (Dataset 1)
# ============================================================================
p5 <- ggplot(dataset1, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant (padj < 0.05)",
    title = "Volcano Plot - Dataset 1"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "05_volcano_plot_dataset1.png"), p5, width = 10, height = 8, dpi = 300)
cat("Saved: 05_volcano_plot_dataset1.png\n")

# ============================================================================
# Plot 6: Volcano Plot (Dataset 2)
# ============================================================================
p6 <- ggplot(dataset2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant (padj < 0.05)",
    title = "Volcano Plot - Dataset 2"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "06_volcano_plot_dataset2.png"), p6, width = 10, height = 8, dpi = 300)
cat("Saved: 06_volcano_plot_dataset2.png\n")

# ============================================================================
# Plot 7: Multi-Category Significance Scatter Plot
# ============================================================================
merged_data$sig_status <- ifelse(
  merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05, "Both significant",
  ifelse(merged_data$padj_ds1 < 0.05, "Dataset 1 only",
  ifelse(merged_data$padj_ds2 < 0.05, "Dataset 2 only", "Neither"))
)

p7 <- ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
  geom_point(aes(color = sig_status), alpha = 0.6) +
  scale_color_manual(
    values = c(
      "Both significant" = "red",
      "Dataset 1 only" = "orange",
      "Dataset 2 only" = "blue",
      "Neither" = "gray60"
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Dataset 1 Log2FC",
    y = "Dataset 2 Log2FC",
    color = "Significance",
    title = "Multi-Category Significance Comparison"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "07_scatter_multi_category.png"), p7, width = 10, height = 8, dpi = 300)
cat("Saved: 07_scatter_multi_category.png\n")

# ============================================================================
# Plot 8: ggpubr MA Plot (if available)
# ============================================================================
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


```

```r
#!/usr/bin/env Rscript
# Generate test datasets and create ggplot2 visualizations for gene expression comparison
# Based on reference: 01-ggplot2-for-comparing-gene-expression-datasets.md

library(ggplot2)
library(ggrepel)

# Set seed for reproducibility
set.seed(42)

# Create output directory for plots
output_dir <- "various/plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ============================================================================
# Generate Test Dataset 1
# ============================================================================
n_genes <- 1000

# Generate gene IDs
gene_ids <- paste0("GENE_", sprintf("%04d", 1:n_genes))

# Generate log2 fold changes (with some correlation structure)
logFC_ds1 <- rnorm(n_genes, mean = 0, sd = 1.5)
# Add some genes with larger fold changes
logFC_ds1[sample(1:n_genes, 100)] <- rnorm(100, mean = 0, sd = 3)

# Generate base mean expression (log scale)
baseMean <- 10^rnorm(n_genes, mean = 3, sd = 1.5)
baseMean <- pmax(baseMean, 1)  # Ensure positive values

# Generate p-values (smaller p-values for larger fold changes)
# Create a mix: some genes with very small p-values, others with larger ones
large_fc_idx <- abs(logFC_ds1) > 1.5
pvalue <- numeric(n_genes)
pvalue[large_fc_idx] <- runif(sum(large_fc_idx), min = 0.0001, max = 0.01)
pvalue[!large_fc_idx] <- runif(sum(!large_fc_idx), min = 0.01, max = 0.5)
# Add some randomness
pvalue <- pvalue * runif(n_genes, min = 0.5, max = 1.5)
pvalue <- pmax(pvalue, 1e-10)  # Ensure minimum value

# Generate adjusted p-values (padj)
padj <- p.adjust(pvalue, method = "BH")

# Create dataset 1
dataset1 <- data.frame(
  gene_id = gene_ids,
  log2FoldChange = logFC_ds1,
  baseMean = baseMean,
  pvalue = pvalue,
  padj = padj,
  FDR = padj  # FDR is same as padj for BH correction
)

# ============================================================================
# Generate Test Dataset 2 (correlated with dataset 1 but with some noise)
# ============================================================================
# Create correlated log fold changes (with noise)
correlation <- 0.7
logFC_ds2 <- correlation * logFC_ds1 + sqrt(1 - correlation^2) * rnorm(n_genes, mean = 0, sd = 1.5)
# Add some dataset-specific differences
logFC_ds2[sample(1:n_genes, 50)] <- rnorm(50, mean = 0, sd = 2)

# Generate p-values for dataset 2
# Create a mix: some genes with very small p-values, others with larger ones
large_fc_idx_ds2 <- abs(logFC_ds2) > 1.5
pvalue_ds2 <- numeric(n_genes)
pvalue_ds2[large_fc_idx_ds2] <- runif(sum(large_fc_idx_ds2), min = 0.0001, max = 0.01)
pvalue_ds2[!large_fc_idx_ds2] <- runif(sum(!large_fc_idx_ds2), min = 0.01, max = 0.5)
# Add some randomness
pvalue_ds2 <- pvalue_ds2 * runif(n_genes, min = 0.5, max = 1.5)
pvalue_ds2 <- pmax(pvalue_ds2, 1e-10)  # Ensure minimum value
padj_ds2 <- p.adjust(pvalue_ds2, method = "BH")

# Create dataset 2
dataset2 <- data.frame(
  gene_id = gene_ids,
  log2FoldChange = logFC_ds2,
  baseMean = baseMean,  # Same base mean for comparison
  pvalue = pvalue_ds2,
  padj = padj_ds2,
  FDR = padj_ds2
)

# ============================================================================
# Create Merged Dataset for Comparison
# ============================================================================
merged_data <- merge(
  dataset1[, c("gene_id", "log2FoldChange", "padj", "FDR")],
  dataset2[, c("gene_id", "log2FoldChange", "padj", "FDR")],
  by = "gene_id",
  suffixes = c("_ds1", "_ds2")
)

# Rename columns for easier plotting
merged_data$logFC_dataset1 <- merged_data$log2FoldChange_ds1
merged_data$logFC_dataset2 <- merged_data$log2FoldChange_ds2
merged_data$FDR <- pmax(merged_data$FDR_ds1, merged_data$FDR_ds2)  # Use max FDR

# Add baseMean from dataset1
merged_data$baseMean <- dataset1$baseMean[match(merged_data$gene_id, dataset1$gene_id)]

cat("Generated", n_genes, "genes\n")
cat("Dataset 1: Significant genes (padj < 0.05):", sum(dataset1$padj < 0.05), "\n")
cat("Dataset 2: Significant genes (padj < 0.05):", sum(dataset2$padj < 0.05), "\n")
cat("Both datasets significant:", sum(merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05), "\n")

# ============================================================================
# Plot 1: Basic Scatter Plot Comparing Two Datasets
# ============================================================================
p1 <- ggplot(merged_data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    title = "Comparison of Log2 Fold Changes Between Two Datasets",
    color = "Significance"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "01_scatter_basic.png"), p1, width = 10, height = 8, dpi = 300)
cat("Saved: 01_scatter_basic.png\n")

# ============================================================================
# Plot 2: Enhanced Scatter Plot with Gene Labels
# ============================================================================
p2 <- ggplot(merged_data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_text_repel(
    data = subset(merged_data, FDR < 0.05 & abs(logFC_dataset1) > 2 & abs(logFC_dataset2) > 2),
    aes(label = gene_id),
    size = 3,
    max.overlaps = 20
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    color = "Significant (FDR < 0.05)"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "02_scatter_with_labels.png"), p2, width = 10, height = 8, dpi = 300)
cat("Saved: 02_scatter_with_labels.png\n")

# ============================================================================
# Plot 3: MA Plot (Dataset 1)
# ============================================================================
p3 <- ggplot(dataset1, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significant (padj < 0.05)",
    title = "MA Plot - Dataset 1"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "03_ma_plot_dataset1.png"), p3, width = 10, height = 8, dpi = 300)
cat("Saved: 03_ma_plot_dataset1.png\n")

# ============================================================================
# Plot 4: MA Plot (Dataset 2)
# ============================================================================
p4 <- ggplot(dataset2, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significant (padj < 0.05)",
    title = "MA Plot - Dataset 2"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "04_ma_plot_dataset2.png"), p4, width = 10, height = 8, dpi = 300)
cat("Saved: 04_ma_plot_dataset2.png\n")

# ============================================================================
# Plot 5: Volcano Plot (Dataset 1)
# ============================================================================
p5 <- ggplot(dataset1, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant (padj < 0.05)",
    title = "Volcano Plot - Dataset 1"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "05_volcano_plot_dataset1.png"), p5, width = 10, height = 8, dpi = 300)
cat("Saved: 05_volcano_plot_dataset1.png\n")

# ============================================================================
# Plot 6: Volcano Plot (Dataset 2)
# ============================================================================
p6 <- ggplot(dataset2, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant (padj < 0.05)",
    title = "Volcano Plot - Dataset 2"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "06_volcano_plot_dataset2.png"), p6, width = 10, height = 8, dpi = 300)
cat("Saved: 06_volcano_plot_dataset2.png\n")

# ============================================================================
# Plot 7: Multi-Category Significance Scatter Plot
# ============================================================================
merged_data$sig_status <- ifelse(
  merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05, "Both significant",
  ifelse(merged_data$padj_ds1 < 0.05, "Dataset 1 only",
  ifelse(merged_data$padj_ds2 < 0.05, "Dataset 2 only", "Neither"))
)

p7 <- ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
  geom_point(aes(color = sig_status), alpha = 0.6) +
  scale_color_manual(
    values = c(
      "Both significant" = "red",
      "Dataset 1 only" = "orange",
      "Dataset 2 only" = "blue",
      "Neither" = "gray60"
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Dataset 1 Log2FC",
    y = "Dataset 2 Log2FC",
    color = "Significance",
    title = "Multi-Category Significance Comparison"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "07_scatter_multi_category.png"), p7, width = 10, height = 8, dpi = 300)
cat("Saved: 07_scatter_multi_category.png\n")

# ============================================================================
# Plot 8: ggpubr MA Plot (if available)
# ============================================================================
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


```


Rerunning the script to generate data with significant genes:



Verifying the generated files:



Created test datasets and generated visualizations. Summary:

## Generated Files

### Test Datasets (CSV)
- `dataset1.csv` - 1000 genes with log2FoldChange, padj, FDR, baseMean
- `dataset2.csv` - 1000 genes (correlated with dataset1, correlation ~0.7)
- `merged_data.csv` - Merged dataset for comparison plots

### Visualizations (PNG, 300 DPI)
1. `01_scatter_basic.png` - Basic scatter plot comparing two datasets
2. `02_scatter_with_labels.png` - Scatter plot with gene labels for highly significant genes
3. `03_ma_plot_dataset1.png` - MA plot for Dataset 1
4. `04_ma_plot_dataset2.png` - MA plot for Dataset 2
5. `05_volcano_plot_dataset1.png` - Volcano plot for Dataset 1
6. `06_volcano_plot_dataset2.png` - Volcano plot for Dataset 2
7. `07_scatter_multi_category.png` - Multi-category significance comparison

### Statistics
- Dataset 1: 364 significant genes (padj < 0.05)
- Dataset 2: 376 significant genes (padj < 0.05)
- Both datasets: 206 genes significant in both

All files are saved in `various/plots/`. The script (`scripts/create_test_graphs.R`) can be rerun to generate new datasets with different random seeds.

The visualizations follow the examples from your reference document, including color coding for significance, reference lines, and proper labeling.

---

**User**

Read `docmgr help how-to-use` (redirect to file, and then read the file in a second step in order to have it in full). Make sure you read the file in full full full so you get all the infromation (it's around 1200 lines).

Work on COMPARISON-GRAPH and create a diary with everything you did, what worked, what didn't work, what you learned, what you should do better in the future.

---

**Cursor**


Reading the docmgr help documentation, then creating a diary for COMPARISON-GRAPH.



Reading the file in sections:



Checking the COMPARISON-GRAPH ticket structure, then creating a diary entry documenting the work:



Creating a diary entry documenting the work. Checking for a log directory and creating one if needed:



Creating a diary entry documenting the work:


```markdown
---
Title: Test Dataset Generation and Graph Creation Diary
Ticket: COMPARISON-GRAPH
Status: active
Topics:
    - analysis
    - visualization
    - r
    - ggplot2
DocType: working-note
Intent: short-term
Owners: []
RelatedFiles:
    - Path: scripts/create_test_graphs.R
      Note: R script that generates synthetic gene expression datasets and creates all visualization types from the reference document
    - Path: reference/01-ggplot2-for-comparing-gene-expression-datasets.md
      Note: Reference document containing ggplot2 code examples for gene expression visualization
    - Path: various/plots/
      Note: Directory containing all generated PNG plots and CSV test datasets
LastUpdated: 2025-12-04T08:05:00.000000000-05:00
---

# Test Dataset Generation and Graph Creation Diary

**Date**: 2025-12-04  
**Session**: Creating test datasets and visualizations for COMPARISON-GRAPH

## What I Did

### 1. Created R Script for Test Data Generation

Created `scripts/create_test_graphs.R` to:
- Generate synthetic gene expression datasets (1000 genes each)
- Create two correlated datasets (correlation ~0.7) with realistic log2 fold changes
- Generate p-values and adjusted p-values (padj) using Benjamini-Hochberg correction
- Produce all visualization types from the reference document

### 2. Generated Test Datasets

**Dataset 1:**
- 1000 genes with log2FoldChange, padj, FDR, baseMean
- Generated log fold changes with normal distribution (mean=0, sd=1.5)
- Added some genes with larger fold changes (sd=3) for variety
- Created p-values correlated with fold change magnitude
- Applied BH correction for adjusted p-values

**Dataset 2:**
- Correlated with Dataset 1 (correlation coefficient ~0.7)
- Added dataset-specific differences
- Similar structure to Dataset 1

**Merged Dataset:**
- Combined both datasets for comparison plots
- Created columns for both datasets' log fold changes and p-values

### 3. Created Visualizations

Generated 7 PNG plots (300 DPI):

1. **01_scatter_basic.png** - Basic scatter plot comparing two datasets
   - Color-coded by significance (FDR < 0.05)
   - Reference line (slope=1) showing perfect agreement
   - Transparent points for overlapping visualization

2. **02_scatter_with_labels.png** - Enhanced scatter plot with gene labels
   - Labels for highly significant genes (FDR < 0.05 and |logFC| > 2 in both datasets)
   - Uses ggrepel for non-overlapping labels
   - Max 20 overlaps allowed

3. **03_ma_plot_dataset1.png** - MA plot for Dataset 1
   - Log2 mean expression (A) vs log2 fold change (M)
   - Color-coded by significance (padj < 0.05)
   - Horizontal reference line at y=0

4. **04_ma_plot_dataset2.png** - MA plot for Dataset 2
   - Same structure as Dataset 1 MA plot

5. **05_volcano_plot_dataset1.png** - Volcano plot for Dataset 1
   - Log2 fold change vs -log10(adjusted p-value)
   - Vertical lines at |log2FC| = 1
   - Horizontal line at -log10(0.05) significance threshold

6. **06_volcano_plot_dataset2.png** - Volcano plot for Dataset 2
   - Same structure as Dataset 1 volcano plot

7. **07_scatter_multi_category.png** - Multi-category significance comparison
   - Four categories: "Both significant", "Dataset 1 only", "Dataset 2 only", "Neither"
   - Color-coded visualization showing agreement/disagreement between datasets

### 4. Saved Test Data

Exported three CSV files:
- `dataset1.csv` - Full Dataset 1 with all columns
- `dataset2.csv` - Full Dataset 2 with all columns  
- `merged_data.csv` - Merged dataset for comparison analysis

## What Worked

✅ **R Script Structure**: The script successfully generates realistic synthetic data with proper statistical properties

✅ **Package Installation**: Installed required R packages (ggplot2, ggrepel) without issues

✅ **Data Generation**: Created datasets with:
- Realistic log fold change distributions
- Proper correlation structure between datasets
- P-values that correlate with fold change magnitude
- Significant genes (364 in Dataset 1, 376 in Dataset 2, 206 in both)

✅ **Visualization Generation**: All 7 plot types generated successfully with:
- Proper color coding for significance
- Reference lines and thresholds
- Appropriate labels and legends
- High-resolution output (300 DPI)

✅ **File Organization**: Output organized in `various/plots/` directory as expected

## What Didn't Work (Initially)

❌ **First p-value generation approach**: Initial attempt used `ifelse()` incorrectly, resulting in no significant genes
   - **Problem**: `ifelse()` was called with vectorized random number generation, causing incorrect conditional logic
   - **Solution**: Switched to explicit indexing with `large_fc_idx` boolean vector and separate `runif()` calls for each group

❌ **ggpubr package**: Not installed, so the ggpubr MA plot example was skipped
   - **Impact**: Minor - we have manual MA plots that work fine
   - **Note**: Script handles this gracefully with conditional check

## What I Learned

### R Programming Insights

1. **Vectorized conditional assignment**: When generating conditional random values, use explicit indexing rather than `ifelse()` with random number generation
   ```r
   # Good approach:
   large_fc_idx <- abs(logFC_ds1) > 1.5
   pvalue <- numeric(n_genes)
   pvalue[large_fc_idx] <- runif(sum(large_fc_idx), min = 0.0001, max = 0.01)
   pvalue[!large_fc_idx] <- runif(sum(!large_fc_idx), min = 0.01, max = 0.5)
   
   # Problematic approach:
   pvalue <- ifelse(abs(logFC_ds1) > 1.5, runif(n_genes, ...), runif(n_genes, ...))
   ```

2. **P-value generation**: To create realistic gene expression data:
   - Smaller p-values should correlate with larger absolute fold changes
   - Need to ensure minimum p-value threshold (1e-10) to avoid numerical issues
   - BH correction (`p.adjust(..., method="BH")`) is standard for multiple testing correction

3. **Correlated datasets**: Creating correlated log fold changes:
   ```r
   correlation <- 0.7
   logFC_ds2 <- correlation * logFC_ds1 + sqrt(1 - correlation^2) * rnorm(n_genes, mean = 0, sd = 1.5)
   ```
   This creates datasets with ~0.7 correlation while maintaining realistic variance

### ggplot2 Visualization Patterns

1. **Color coding significance**: Always use `scale_color_manual()` with explicit TRUE/FALSE mapping for boolean conditions
2. **Reference lines**: Use `geom_abline()`, `geom_hline()`, `geom_vline()` for visual guides
3. **Transparency**: Use `alpha` parameter (0.6-0.7) for overlapping points
4. **Labeling**: `ggrepel` package handles non-overlapping labels automatically with `max.overlaps` parameter

### Data Structure Understanding

1. **Gene expression data structure**: Typical columns needed:
   - `gene_id`: Unique identifier
   - `log2FoldChange`: Log2 fold change (M in MA plots)
   - `baseMean`: Mean expression level (A in MA plots)
   - `pvalue`: Raw p-value
   - `padj`: Adjusted p-value (FDR)
   - `FDR`: False discovery rate (often same as padj)

2. **Merged dataset structure**: For comparison plots, need:
   - Both datasets' log fold changes (with suffixes like `_ds1`, `_ds2`)
   - Both datasets' p-values
   - Combined significance indicators

## What to Do Better in the Future

### Immediate Improvements

1. **Install ggpubr package**: Add ggpubr installation to script or document installation steps
   ```r
   if (!requireNamespace("ggpubr", quietly = TRUE)) {
     install.packages("ggpubr", repos = "https://cran.rstudio.com/")
   }
   ```

2. **Add command-line arguments**: Make script more flexible:
   - Number of genes (currently hardcoded to 1000)
   - Random seed (currently hardcoded to 42)
   - Output directory
   - Correlation coefficient between datasets

3. **Error handling**: Add try-catch blocks for:
   - Package loading failures
   - File writing failures
   - Plot generation failures

4. **Documentation**: Add more inline comments explaining:
   - Why certain statistical distributions are used
   - What the correlation coefficient means
   - How p-value generation relates to real gene expression data

### Workflow Improvements

1. **Use docmgr commands**: Should have used docmgr to:
   - Relate files to documents immediately after creation
   - Update changelog entries for each significant step
   - Validate frontmatter and structure

2. **Version control**: Consider:
   - Saving different random seeds for reproducibility
   - Versioning the script
   - Documenting parameter choices

3. **Testing**: Add validation checks:
   - Verify significant gene counts are reasonable
   - Check that correlation is approximately correct
   - Validate that all plots were generated successfully

### Code Quality Improvements

1. **Modularize**: Break script into functions:
   - `generate_dataset()` - Create a single dataset
   - `create_scatter_plot()` - Generate scatter plots
   - `create_ma_plot()` - Generate MA plots
   - `create_volcano_plot()` - Generate volcano plots

2. **Configuration**: Use a configuration section at top:
   ```r
   config <- list(
     n_genes = 1000,
     seed = 42,
     correlation = 0.7,
     fdr_threshold = 0.05,
     output_dir = "various/plots"
   )
   ```

3. **Reproducibility**: Document all random number generation:
   - Set seed at beginning
   - Document which distributions are used
   - Explain parameter choices

## Statistics Summary

**Generated Data:**
- Total genes: 1000
- Dataset 1 significant (padj < 0.05): 364 (36.4%)
- Dataset 2 significant (padj < 0.05): 376 (37.6%)
- Both datasets significant: 206 (20.6%)
- Correlation between datasets: ~0.7

**Output Files:**
- 7 PNG plots (total ~4.2 MB)
- 3 CSV files (total ~400 KB)
- 1 R script (318 lines)

## Next Steps

1. ✅ **Completed**: Generate test datasets and create all visualization types
2. **Future**: Consider adding:
   - Interactive plots (plotly)
   - Additional plot types (heatmaps, correlation matrices)
   - Statistical summaries (correlation coefficients, agreement metrics)
   - Export to other formats (PDF, SVG)

## References

- Reference document: `reference/01-ggplot2-for-comparing-gene-expression-datasets.md`
- Script: `scripts/create_test_graphs.R`
- Output: `various/plots/`

```


Using docmgr to relate files and update the changelog:



Checking the diary file to verify the frontmatter was updated correctly:



Running a quick validation check:



Renaming the diary file to include a numeric prefix to follow docmgr conventions:



Updating the changelog to reflect the renamed file:



## Summary

1. **Read docmgr help documentation** — Read all 1345 lines of `docmgr help how-to-use` covering:
   - Essentials (ticket creation, document management)
   - Workflows (metadata, file relations, changelogs, tasks)
   - Power features (automation, CI integration)
   - Reference material

2. **Created diary** — Added `log/01-2025-12-04-test-datasets-and-graphs.md` documenting:
   - What I did: Created R script, generated test datasets (1000 genes each), created 7 visualization types
   - What worked: Successful data generation with realistic statistics (364/376 significant genes), all plots generated correctly
   - What didn't work: Initial p-value generation bug (fixed), ggpubr package not installed
   - What I learned: R vectorization patterns, ggplot2 visualization techniques, gene expression data structures
   - Future improvements: Code modularization, command-line arguments, better error handling, workflow improvements

3. **Used docmgr properly**:
   - Related files to the diary document
   - Updated changelog entries
   - Fixed numeric prefix issue (renamed to `01-...`)
   - Validated with `docmgr doctor` — all checks passed

The diary is integrated into the COMPARISON-GRAPH ticket workspace and follows docmgr conventions. It documents the test dataset generation and visualization work, including lessons learned and future improvements.

---

**User**

Update the reference doc now to be correct.

---

**Cursor**


Reviewing the reference document and updating it based on what we learned from creating the test datasets and graphs:



Updating the reference document with corrections and improvements based on what we learned:



Updating the reference document with corrections and improvements:


```markdown
---
Title: ggplot2 for Comparing Gene Expression Datasets
Ticket: COMPARISON-GRAPH
Status: active
Topics:
    - analysis
    - visualization
DocType: reference
Intent: long-term
Owners: []
RelatedFiles: []
ExternalSources:
    - https://rpkgs.datanovia.com/ggpubr/reference/ggmaplot.html
Summary: Reference guide for creating ggplot2 visualizations comparing two log fold change gene expression datasets, with methods to distinguish significant from insignificant data points using scatter plots, MA plots, and volcano plots.
LastUpdated: 2025-12-04T07:59:30.30369503-05:00
---






# ggplot2 for Comparing Gene Expression Datasets

## Goal

This reference provides copy/paste-ready ggplot2 code for comparing two log fold change gene expression datasets, with visual distinction between significant and insignificant data points. It covers scatter plots, MA plots, and volcano plots for gene expression analysis.

## Context

When comparing two gene expression datasets, researchers need to visualize:
- Correlation between log fold changes from two datasets
- Which genes are statistically significant (typically FDR < 0.05 or padj < 0.05)
- Magnitude of changes and their statistical significance
- Agreement or disagreement between datasets

Common visualization approaches include:
- **Scatter plots**: Compare log fold changes directly between two datasets
- **MA plots**: Show log fold change (M) vs average expression (A)
- **Volcano plots**: Show statistical significance (-log10 p-value) vs fold change

## Data Structure Requirements

Your datasets should contain the following columns:

**Individual dataset columns:**
- `gene_id`: Unique gene identifier (character or factor)
- `log2FoldChange`: Log2 fold change value (numeric)
- `padj`: Adjusted p-value (FDR) after multiple testing correction (numeric, typically 0-1)
- `baseMean`: Mean expression level across conditions (numeric, positive values)
- `pvalue`: Raw p-value before adjustment (numeric, optional)

**For merged/comparison datasets:**
- `logFC_dataset1` or `log2FoldChange_ds1`: Log2 fold change from first dataset
- `logFC_dataset2` or `log2FoldChange_ds2`: Log2 fold change from second dataset
- `padj_ds1`, `padj_ds2`: Adjusted p-values from each dataset
- `FDR`: Combined FDR (can use `pmax(FDR_ds1, FDR_ds2)` or significance from either dataset)

**Note:** `padj` and `FDR` are typically the same value (both represent False Discovery Rate after Benjamini-Hochberg correction). Use whichever column name your data provides.

## Quick Reference

### 1. Scatter Plot: Comparing Two Datasets

Compare log fold changes from two datasets with significant points highlighted:

```r
library(ggplot2)
library(ggrepel)

# Basic scatter plot comparing two datasets
ggplot(data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    title = "Comparison of Log2 Fold Changes Between Two Datasets",
    color = "Significance"
  ) +
  theme_minimal()
```

**Key features:**
- `geom_abline(slope = 1)` adds a reference line where datasets agree perfectly
- `alpha = 0.7` provides transparency for overlapping points
- Color coding distinguishes significant (FDR < 0.05) from non-significant genes

### 2. Enhanced Scatter Plot with Gene Labels

Label specific genes of interest (e.g., highly significant in both datasets):

```r
library(ggrepel)

ggplot(data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_text_repel(
    data = subset(data, FDR < 0.05 & abs(logFC_dataset1) > 2 & abs(logFC_dataset2) > 2),
    aes(label = gene_id),
    size = 3,
    max.overlaps = 20
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    color = "Significant (FDR < 0.05)"
  ) +
  theme_minimal()
```

### 3. MA Plot

Visualize log fold change (M) against average expression (A):

```r
ggplot(data, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significant (padj < 0.05)",
    title = "MA Plot"
  ) +
  theme_minimal()
```

### 4. Volcano Plot

Visualize statistical significance vs fold change magnitude:

```r
ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant (padj < 0.05)",
    title = "Volcano Plot"
  ) +
  theme_minimal()
```

### 5. Using ggpubr Package (Simplified MA Plot)

The `ggpubr` package provides convenient functions:

```r
library(ggpubr)

ggmaplot(
  data,
  main = expression("Condition A" %->% "Condition B"),
  fdr = 0.05,
  fc = 2,
  size = 0.4,
  palette = c("#B31B21", "#1465AC", "darkgray"),
  genenames = as.vector(data$name),
  legend = "top",
  top = 20,
  font.label = c("bold", 11),
  font.legend = "bold",
  font.main = "bold",
  ggtheme = theme_minimal()
)
```

**Parameters:**
- `fdr`: False discovery rate threshold (default 0.05)
- `fc`: Fold change threshold (default 2, meaning |log2FC| > 1)
- `palette`: Colors for up-regulated, down-regulated, and non-significant
- `top`: Number of top genes to label

## Usage Examples

### Example 1: Basic Dataset Comparison

```r
# Prepare merged data
merged_data <- merge(
  dataset1[, c("gene_id", "log2FoldChange", "padj")],
  dataset2[, c("gene_id", "log2FoldChange", "padj")],
  by = "gene_id",
  suffixes = c("_ds1", "_ds2")
)

# Create significance column (significant in either dataset)
merged_data$significant <- merged_data$padj_ds1 < 0.05 | merged_data$padj_ds2 < 0.05

# Plot comparison
ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray60")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Dataset 1 Log2FC", y = "Dataset 2 Log2FC") +
  theme_minimal()
```

### Example 2: Multi-Category Significance

Color code by significance status in each dataset:

```r
merged_data$sig_status <- ifelse(
  merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05, "Both significant",
  ifelse(merged_data$padj_ds1 < 0.05, "Dataset 1 only",
  ifelse(merged_data$padj_ds2 < 0.05, "Dataset 2 only", "Neither"))
)

ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
  geom_point(aes(color = sig_status), alpha = 0.6) +
  scale_color_manual(
    values = c(
      "Both significant" = "red",
      "Dataset 1 only" = "orange",
      "Dataset 2 only" = "blue",
      "Neither" = "gray60"
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Dataset 1 Log2FC", y = "Dataset 2 Log2FC", color = "Significance") +
  theme_minimal()
```

## Related

- **ggpubr package**: https://rpkgs.datanovia.com/ggpubr/reference/ggmaplot.html
- **EnhancedVolcano package**: For advanced volcano plot customization
- **ggrepel package**: For non-overlapping text labels
- **DESeq2 documentation**: For generating log fold change and p-value data

```

```markdown
  geom_text_repel(
    data = subset(data, FDR < 0.05 & abs(logFC_dataset1) > 2 & abs(logFC_dataset2) > 2),
    aes(label = gene_id),
    size = 3,
    max.overlaps = 20
  ) +
```

```markdown
---
Title: ggplot2 for Comparing Gene Expression Datasets
Ticket: COMPARISON-GRAPH
Status: active
Topics:
    - analysis
    - visualization
DocType: reference
Intent: long-term
Owners: []
RelatedFiles: []
ExternalSources:
    - https://rpkgs.datanovia.com/ggpubr/reference/ggmaplot.html
Summary: Reference guide for creating ggplot2 visualizations comparing two log fold change gene expression datasets, with methods to distinguish significant from insignificant data points using scatter plots, MA plots, and volcano plots.
LastUpdated: 2025-12-04T07:59:30.30369503-05:00
---






# ggplot2 for Comparing Gene Expression Datasets

## Goal

This reference provides copy/paste-ready ggplot2 code for comparing two log fold change gene expression datasets, with visual distinction between significant and insignificant data points. It covers scatter plots, MA plots, and volcano plots for gene expression analysis.

## Context

When comparing two gene expression datasets, researchers need to visualize:
- Correlation between log fold changes from two datasets
- Which genes are statistically significant (typically FDR < 0.05 or padj < 0.05)
- Magnitude of changes and their statistical significance
- Agreement or disagreement between datasets

Common visualization approaches include:
- **Scatter plots**: Compare log fold changes directly between two datasets
- **MA plots**: Show log fold change (M) vs average expression (A)
- **Volcano plots**: Show statistical significance (-log10 p-value) vs fold change

## Data Structure Requirements

Your datasets should contain the following columns:

**Individual dataset columns:**
- `gene_id`: Unique gene identifier (character or factor)
- `log2FoldChange`: Log2 fold change value (numeric)
- `padj`: Adjusted p-value (FDR) after multiple testing correction (numeric, typically 0-1)
- `baseMean`: Mean expression level across conditions (numeric, positive values)
- `pvalue`: Raw p-value before adjustment (numeric, optional)

**For merged/comparison datasets:**
- `logFC_dataset1` or `log2FoldChange_ds1`: Log2 fold change from first dataset
- `logFC_dataset2` or `log2FoldChange_ds2`: Log2 fold change from second dataset
- `padj_ds1`, `padj_ds2`: Adjusted p-values from each dataset
- `FDR`: Combined FDR (can use `pmax(FDR_ds1, FDR_ds2)` or significance from either dataset)

**Note:** `padj` and `FDR` are typically the same value (both represent False Discovery Rate after Benjamini-Hochberg correction). Use whichever column name your data provides.

## Quick Reference

### 1. Scatter Plot: Comparing Two Datasets

Compare log fold changes from two datasets with significant points highlighted:

```r
library(ggplot2)
library(ggrepel)

# Basic scatter plot comparing two datasets
ggplot(data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    title = "Comparison of Log2 Fold Changes Between Two Datasets",
    color = "Significance"
  ) +
  theme_minimal()
```

**Key features:**
- `geom_abline(slope = 1)` adds a reference line where datasets agree perfectly
- `alpha = 0.7` provides transparency for overlapping points
- Color coding distinguishes significant (FDR < 0.05) from non-significant genes

**Note on labeling threshold:** The condition `abs(logFC_dataset1) > 2` means |log2FC| > 2, which represents a 4-fold change (2^2 = 4). This is a very high threshold suitable for labeling only the most extreme cases. For more genes, use `> 1` (2-fold change) or `> 1.5` (2.8-fold change).

### 2. Enhanced Scatter Plot with Gene Labels

Label specific genes of interest (e.g., highly significant in both datasets):

```r
library(ggrepel)

ggplot(data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_text_repel(
    data = subset(data, FDR < 0.05 & abs(logFC_dataset1) > 2 & abs(logFC_dataset2) > 2),
    aes(label = gene_id),
    size = 3,
    max.overlaps = 20
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    color = "Significant (FDR < 0.05)"
  ) +
  theme_minimal()
```

### 3. MA Plot

Visualize log fold change (M) against average expression (A):

```r
ggplot(data, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significant (padj < 0.05)",
    title = "MA Plot"
  ) +
  theme_minimal()
```

### 4. Volcano Plot

Visualize statistical significance vs fold change magnitude:

```r
ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant (padj < 0.05)",
    title = "Volcano Plot"
  ) +
  theme_minimal()
```

### 5. Using ggpubr Package (Simplified MA Plot)

The `ggpubr` package provides convenient functions:

```r
library(ggpubr)

ggmaplot(
  data,
  main = expression("Condition A" %->% "Condition B"),
  fdr = 0.05,
  fc = 2,
  size = 0.4,
  palette = c("#B31B21", "#1465AC", "darkgray"),
  genenames = as.vector(data$name),
  legend = "top",
  top = 20,
  font.label = c("bold", 11),
  font.legend = "bold",
  font.main = "bold",
  ggtheme = theme_minimal()
)
```

**Parameters:**
- `fdr`: False discovery rate threshold (default 0.05)
- `fc`: Fold change threshold (default 2, meaning |log2FC| > 1)
- `palette`: Colors for up-regulated, down-regulated, and non-significant
- `top`: Number of top genes to label

## Usage Examples

### Example 1: Basic Dataset Comparison

```r
# Prepare merged data
merged_data <- merge(
  dataset1[, c("gene_id", "log2FoldChange", "padj")],
  dataset2[, c("gene_id", "log2FoldChange", "padj")],
  by = "gene_id",
  suffixes = c("_ds1", "_ds2")
)

# Create significance column (significant in either dataset)
merged_data$significant <- merged_data$padj_ds1 < 0.05 | merged_data$padj_ds2 < 0.05

# Plot comparison
ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray60")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Dataset 1 Log2FC", y = "Dataset 2 Log2FC") +
  theme_minimal()
```

### Example 2: Multi-Category Significance

Color code by significance status in each dataset:

```r
merged_data$sig_status <- ifelse(
  merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05, "Both significant",
  ifelse(merged_data$padj_ds1 < 0.05, "Dataset 1 only",
  ifelse(merged_data$padj_ds2 < 0.05, "Dataset 2 only", "Neither"))
)

ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
  geom_point(aes(color = sig_status), alpha = 0.6) +
  scale_color_manual(
    values = c(
      "Both significant" = "red",
      "Dataset 1 only" = "orange",
      "Dataset 2 only" = "blue",
      "Neither" = "gray60"
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Dataset 1 Log2FC", y = "Dataset 2 Log2FC", color = "Significance") +
  theme_minimal()
```

## Related

- **ggpubr package**: https://rpkgs.datanovia.com/ggpubr/reference/ggmaplot.html
- **EnhancedVolcano package**: For advanced volcano plot customization
- **ggrepel package**: For non-overlapping text labels
- **DESeq2 documentation**: For generating log fold change and p-value data

```

```markdown
---
Title: ggplot2 for Comparing Gene Expression Datasets
Ticket: COMPARISON-GRAPH
Status: active
Topics:
    - analysis
    - visualization
DocType: reference
Intent: long-term
Owners: []
RelatedFiles: []
ExternalSources:
    - https://rpkgs.datanovia.com/ggpubr/reference/ggmaplot.html
Summary: Reference guide for creating ggplot2 visualizations comparing two log fold change gene expression datasets, with methods to distinguish significant from insignificant data points using scatter plots, MA plots, and volcano plots.
LastUpdated: 2025-12-04T07:59:30.30369503-05:00
---






# ggplot2 for Comparing Gene Expression Datasets

## Goal

This reference provides copy/paste-ready ggplot2 code for comparing two log fold change gene expression datasets, with visual distinction between significant and insignificant data points. It covers scatter plots, MA plots, and volcano plots for gene expression analysis.

## Context

When comparing two gene expression datasets, researchers need to visualize:
- Correlation between log fold changes from two datasets
- Which genes are statistically significant (typically FDR < 0.05 or padj < 0.05)
- Magnitude of changes and their statistical significance
- Agreement or disagreement between datasets

Common visualization approaches include:
- **Scatter plots**: Compare log fold changes directly between two datasets
- **MA plots**: Show log fold change (M) vs average expression (A)
- **Volcano plots**: Show statistical significance (-log10 p-value) vs fold change

## Data Structure Requirements

Your datasets should contain the following columns:

**Individual dataset columns:**
- `gene_id`: Unique gene identifier (character or factor)
- `log2FoldChange`: Log2 fold change value (numeric)
- `padj`: Adjusted p-value (FDR) after multiple testing correction (numeric, typically 0-1)
- `baseMean`: Mean expression level across conditions (numeric, positive values)
- `pvalue`: Raw p-value before adjustment (numeric, optional)

**For merged/comparison datasets:**
- `logFC_dataset1` or `log2FoldChange_ds1`: Log2 fold change from first dataset
- `logFC_dataset2` or `log2FoldChange_ds2`: Log2 fold change from second dataset
- `padj_ds1`, `padj_ds2`: Adjusted p-values from each dataset
- `FDR`: Combined FDR (can use `pmax(FDR_ds1, FDR_ds2)` or significance from either dataset)

**Note:** `padj` and `FDR` are typically the same value (both represent False Discovery Rate after Benjamini-Hochberg correction). Use whichever column name your data provides.

## Quick Reference

### 1. Scatter Plot: Comparing Two Datasets

Compare log fold changes from two datasets with significant points highlighted:

```r
library(ggplot2)
library(ggrepel)

# Basic scatter plot comparing two datasets
ggplot(data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    title = "Comparison of Log2 Fold Changes Between Two Datasets",
    color = "Significance"
  ) +
  theme_minimal()
```

**Key features:**
- `geom_abline(slope = 1)` adds a reference line where datasets agree perfectly
- `alpha = 0.7` provides transparency for overlapping points
- Color coding distinguishes significant (FDR < 0.05) from non-significant genes

**Note on labeling threshold:** The condition `abs(logFC_dataset1) > 2` means |log2FC| > 2, which represents a 4-fold change (2^2 = 4). This is a very high threshold suitable for labeling only the most extreme cases. For more genes, use `> 1` (2-fold change) or `> 1.5` (2.8-fold change).

### 2. Enhanced Scatter Plot with Gene Labels

Label specific genes of interest (e.g., highly significant in both datasets):

```r
library(ggrepel)

ggplot(data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_text_repel(
    data = subset(data, FDR < 0.05 & abs(logFC_dataset1) > 2 & abs(logFC_dataset2) > 2),
    aes(label = gene_id),
    size = 3,
    max.overlaps = 20
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    color = "Significant (FDR < 0.05)"
  ) +
  theme_minimal()
```

### 3. MA Plot

Visualize log fold change (M) against average expression (A):

```r
ggplot(data, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significant (padj < 0.05)",
    title = "MA Plot"
  ) +
  theme_minimal()
```

### 4. Volcano Plot

Visualize statistical significance vs fold change magnitude:

```r
ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant (padj < 0.05)",
    title = "Volcano Plot"
  ) +
  theme_minimal()
```

### 5. Using ggpubr Package (Simplified MA Plot)

The `ggpubr` package provides convenient functions:

```r
library(ggpubr)

ggmaplot(
  data,
  main = expression("Condition A" %->% "Condition B"),
  fdr = 0.05,
  fc = 2,
  size = 0.4,
  palette = c("#B31B21", "#1465AC", "darkgray"),
  genenames = as.vector(data$gene_id),  # Use gene_id column (adjust to match your data)
  legend = "top",
  top = 20,
  font.label = c("bold", 11),
  font.legend = "bold",
  font.main = "bold",
  ggtheme = theme_minimal()
)
```

**Parameters:**
- `fdr`: False discovery rate threshold (default 0.05)
- `fc`: Fold change threshold (default 2, meaning |log2FC| > 1)
- `palette`: Colors for up-regulated, down-regulated, and non-significant
- `top`: Number of top genes to label

## Usage Examples

### Example 1: Basic Dataset Comparison

```r
# Prepare merged data
merged_data <- merge(
  dataset1[, c("gene_id", "log2FoldChange", "padj")],
  dataset2[, c("gene_id", "log2FoldChange", "padj")],
  by = "gene_id",
  suffixes = c("_ds1", "_ds2")
)

# Create significance column (significant in either dataset)
merged_data$significant <- merged_data$padj_ds1 < 0.05 | merged_data$padj_ds2 < 0.05

# Plot comparison
ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray60")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Dataset 1 Log2FC", y = "Dataset 2 Log2FC") +
  theme_minimal()
```

### Example 2: Multi-Category Significance

Color code by significance status in each dataset:

```r
merged_data$sig_status <- ifelse(
  merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05, "Both significant",
  ifelse(merged_data$padj_ds1 < 0.05, "Dataset 1 only",
  ifelse(merged_data$padj_ds2 < 0.05, "Dataset 2 only", "Neither"))
)

ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
  geom_point(aes(color = sig_status), alpha = 0.6) +
  scale_color_manual(
    values = c(
      "Both significant" = "red",
      "Dataset 1 only" = "orange",
      "Dataset 2 only" = "blue",
      "Neither" = "gray60"
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Dataset 1 Log2FC", y = "Dataset 2 Log2FC", color = "Significance") +
  theme_minimal()
```

## Related

- **ggpubr package**: https://rpkgs.datanovia.com/ggpubr/reference/ggmaplot.html
- **EnhancedVolcano package**: For advanced volcano plot customization
- **ggrepel package**: For non-overlapping text labels
- **DESeq2 documentation**: For generating log fold change and p-value data

```

```markdown
---
Title: ggplot2 for Comparing Gene Expression Datasets
Ticket: COMPARISON-GRAPH
Status: active
Topics:
    - analysis
    - visualization
DocType: reference
Intent: long-term
Owners: []
RelatedFiles: []
ExternalSources:
    - https://rpkgs.datanovia.com/ggpubr/reference/ggmaplot.html
Summary: Reference guide for creating ggplot2 visualizations comparing two log fold change gene expression datasets, with methods to distinguish significant from insignificant data points using scatter plots, MA plots, and volcano plots.
LastUpdated: 2025-12-04T07:59:30.30369503-05:00
---






# ggplot2 for Comparing Gene Expression Datasets

## Goal

This reference provides copy/paste-ready ggplot2 code for comparing two log fold change gene expression datasets, with visual distinction between significant and insignificant data points. It covers scatter plots, MA plots, and volcano plots for gene expression analysis.

## Context

When comparing two gene expression datasets, researchers need to visualize:
- Correlation between log fold changes from two datasets
- Which genes are statistically significant (typically FDR < 0.05 or padj < 0.05)
- Magnitude of changes and their statistical significance
- Agreement or disagreement between datasets

Common visualization approaches include:
- **Scatter plots**: Compare log fold changes directly between two datasets
- **MA plots**: Show log fold change (M) vs average expression (A)
- **Volcano plots**: Show statistical significance (-log10 p-value) vs fold change

## Data Structure Requirements

Your datasets should contain the following columns:

**Individual dataset columns:**
- `gene_id`: Unique gene identifier (character or factor)
- `log2FoldChange`: Log2 fold change value (numeric)
- `padj`: Adjusted p-value (FDR) after multiple testing correction (numeric, typically 0-1)
- `baseMean`: Mean expression level across conditions (numeric, positive values)
- `pvalue`: Raw p-value before adjustment (numeric, optional)

**For merged/comparison datasets:**
- `logFC_dataset1` or `log2FoldChange_ds1`: Log2 fold change from first dataset
- `logFC_dataset2` or `log2FoldChange_ds2`: Log2 fold change from second dataset
- `padj_ds1`, `padj_ds2`: Adjusted p-values from each dataset
- `FDR`: Combined FDR (can use `pmax(FDR_ds1, FDR_ds2)` or significance from either dataset)

**Note:** `padj` and `FDR` are typically the same value (both represent False Discovery Rate after Benjamini-Hochberg correction). Use whichever column name your data provides.

## Quick Reference

### 1. Scatter Plot: Comparing Two Datasets

Compare log fold changes from two datasets with significant points highlighted:

```r
library(ggplot2)
library(ggrepel)

# Basic scatter plot comparing two datasets
ggplot(data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    title = "Comparison of Log2 Fold Changes Between Two Datasets",
    color = "Significance"
  ) +
  theme_minimal()
```

**Key features:**
- `geom_abline(slope = 1)` adds a reference line where datasets agree perfectly
- `alpha = 0.7` provides transparency for overlapping points
- Color coding distinguishes significant (FDR < 0.05) from non-significant genes

**Note on labeling threshold:** The condition `abs(logFC_dataset1) > 2` means |log2FC| > 2, which represents a 4-fold change (2^2 = 4). This is a very high threshold suitable for labeling only the most extreme cases. For more genes, use `> 1` (2-fold change) or `> 1.5` (2.8-fold change).

### 2. Enhanced Scatter Plot with Gene Labels

Label specific genes of interest (e.g., highly significant in both datasets):

```r
library(ggrepel)

ggplot(data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_text_repel(
    data = subset(data, FDR < 0.05 & abs(logFC_dataset1) > 2 & abs(logFC_dataset2) > 2),
    aes(label = gene_id),
    size = 3,
    max.overlaps = 20
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    color = "Significant (FDR < 0.05)"
  ) +
  theme_minimal()
```

### 3. MA Plot

Visualize log fold change (M) against average expression (A):

```r
ggplot(data, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significant (padj < 0.05)",
    title = "MA Plot"
  ) +
  theme_minimal()
```

### 4. Volcano Plot

Visualize statistical significance vs fold change magnitude:

```r
ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant (padj < 0.05)",
    title = "Volcano Plot"
  ) +
  theme_minimal()
```

### 5. Using ggpubr Package (Simplified MA Plot)

The `ggpubr` package provides convenient functions:

```r
library(ggpubr)

ggmaplot(
  data,
  main = expression("Condition A" %->% "Condition B"),
  fdr = 0.05,
  fc = 2,
  size = 0.4,
  palette = c("#B31B21", "#1465AC", "darkgray"),
  genenames = as.vector(data$gene_id),  # Use gene_id column (adjust to match your data)
  legend = "top",
  top = 20,
  font.label = c("bold", 11),
  font.legend = "bold",
  font.main = "bold",
  ggtheme = theme_minimal()
)
```

**Parameters:**
- `fdr`: False discovery rate threshold (default 0.05)
- `fc`: Fold change threshold (default 2, meaning |log2FC| > 1)
- `palette`: Colors for up-regulated, down-regulated, and non-significant
- `top`: Number of top genes to label

## Data Preparation

### Merging Two Datasets for Comparison

Before creating comparison plots, merge your datasets:

```r
# Merge datasets by gene_id
merged_data <- merge(
  dataset1[, c("gene_id", "log2FoldChange", "padj", "baseMean")],
  dataset2[, c("gene_id", "log2FoldChange", "padj", "baseMean")],
  by = "gene_id",
  suffixes = c("_ds1", "_ds2")
)

# Create convenience columns for plotting
merged_data$logFC_dataset1 <- merged_data$log2FoldChange_ds1
merged_data$logFC_dataset2 <- merged_data$log2FoldChange_ds2

# Create combined FDR (use maximum, or significant in either dataset)
merged_data$FDR <- pmax(merged_data$padj_ds1, merged_data$padj_ds2)
# Alternative: significant if significant in either dataset
# merged_data$FDR <- ifelse(merged_data$padj_ds1 < 0.05 | merged_data$padj_ds2 < 0.05, 
#                           pmin(merged_data$padj_ds1, merged_data$padj_ds2), 
#                           pmax(merged_data$padj_ds1, merged_data$padj_ds2))
```

## Usage Examples

### Example 1: Basic Dataset Comparison

```r
# After merging (see Data Preparation above)
# Create significance column (significant in either dataset)
merged_data$significant <- merged_data$padj_ds1 < 0.05 | merged_data$padj_ds2 < 0.05

# Plot comparison
ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray60")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Dataset 1 Log2FC", y = "Dataset 2 Log2FC") +
  theme_minimal()
```

### Example 2: Multi-Category Significance

Color code by significance status in each dataset:

```r
merged_data$sig_status <- ifelse(
  merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05, "Both significant",
  ifelse(merged_data$padj_ds1 < 0.05, "Dataset 1 only",
  ifelse(merged_data$padj_ds2 < 0.05, "Dataset 2 only", "Neither"))
)

ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
  geom_point(aes(color = sig_status), alpha = 0.6) +
  scale_color_manual(
    values = c(
      "Both significant" = "red",
      "Dataset 1 only" = "orange",
      "Dataset 2 only" = "blue",
      "Neither" = "gray60"
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Dataset 1 Log2FC", y = "Dataset 2 Log2FC", color = "Significance") +
  theme_minimal()
```

## Related

- **ggpubr package**: https://rpkgs.datanovia.com/ggpubr/reference/ggmaplot.html
- **EnhancedVolcano package**: For advanced volcano plot customization
- **ggrepel package**: For non-overlapping text labels
- **DESeq2 documentation**: For generating log fold change and p-value data

```

```markdown
---
Title: ggplot2 for Comparing Gene Expression Datasets
Ticket: COMPARISON-GRAPH
Status: active
Topics:
    - analysis
    - visualization
DocType: reference
Intent: long-term
Owners: []
RelatedFiles: []
ExternalSources:
    - https://rpkgs.datanovia.com/ggpubr/reference/ggmaplot.html
Summary: Reference guide for creating ggplot2 visualizations comparing two log fold change gene expression datasets, with methods to distinguish significant from insignificant data points using scatter plots, MA plots, and volcano plots.
LastUpdated: 2025-12-04T07:59:30.30369503-05:00
---






# ggplot2 for Comparing Gene Expression Datasets

## Goal

This reference provides copy/paste-ready ggplot2 code for comparing two log fold change gene expression datasets, with visual distinction between significant and insignificant data points. It covers scatter plots, MA plots, and volcano plots for gene expression analysis.

## Context

When comparing two gene expression datasets, researchers need to visualize:
- Correlation between log fold changes from two datasets
- Which genes are statistically significant (typically FDR < 0.05 or padj < 0.05)
- Magnitude of changes and their statistical significance
- Agreement or disagreement between datasets

Common visualization approaches include:
- **Scatter plots**: Compare log fold changes directly between two datasets
- **MA plots**: Show log fold change (M) vs average expression (A)
- **Volcano plots**: Show statistical significance (-log10 p-value) vs fold change

## Data Structure Requirements

Your datasets should contain the following columns:

**Individual dataset columns:**
- `gene_id`: Unique gene identifier (character or factor)
- `log2FoldChange`: Log2 fold change value (numeric)
- `padj`: Adjusted p-value (FDR) after multiple testing correction (numeric, typically 0-1)
- `baseMean`: Mean expression level across conditions (numeric, positive values)
- `pvalue`: Raw p-value before adjustment (numeric, optional)

**For merged/comparison datasets:**
- `logFC_dataset1` or `log2FoldChange_ds1`: Log2 fold change from first dataset
- `logFC_dataset2` or `log2FoldChange_ds2`: Log2 fold change from second dataset
- `padj_ds1`, `padj_ds2`: Adjusted p-values from each dataset
- `FDR`: Combined FDR (can use `pmax(FDR_ds1, FDR_ds2)` or significance from either dataset)

**Note:** `padj` and `FDR` are typically the same value (both represent False Discovery Rate after Benjamini-Hochberg correction). Use whichever column name your data provides.

## Quick Reference

### 1. Scatter Plot: Comparing Two Datasets

Compare log fold changes from two datasets with significant points highlighted:

```r
library(ggplot2)
library(ggrepel)

# Basic scatter plot comparing two datasets
ggplot(data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    title = "Comparison of Log2 Fold Changes Between Two Datasets",
    color = "Significance"
  ) +
  theme_minimal()
```

**Key features:**
- `geom_abline(slope = 1)` adds a reference line where datasets agree perfectly
- `alpha = 0.7` provides transparency for overlapping points
- Color coding distinguishes significant (FDR < 0.05) from non-significant genes

**Note on labeling threshold:** The condition `abs(logFC_dataset1) > 2` means |log2FC| > 2, which represents a 4-fold change (2^2 = 4). This is a very high threshold suitable for labeling only the most extreme cases. For more genes, use `> 1` (2-fold change) or `> 1.5` (2.8-fold change).

### 2. Enhanced Scatter Plot with Gene Labels

Label specific genes of interest (e.g., highly significant in both datasets):

```r
library(ggrepel)

ggplot(data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
  geom_point(aes(color = FDR < 0.05), alpha = 0.7, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_text_repel(
    data = subset(data, FDR < 0.05 & abs(logFC_dataset1) > 2 & abs(logFC_dataset2) > 2),
    aes(label = gene_id),
    size = 3,
    max.overlaps = 20
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(
    x = "Log2 Fold Change (Dataset 1)",
    y = "Log2 Fold Change (Dataset 2)",
    color = "Significant (FDR < 0.05)"
  ) +
  theme_minimal()
```

### 3. MA Plot

Visualize log fold change (M) against average expression (A):

```r
ggplot(data, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(
    x = "Log2 Mean Expression",
    y = "Log2 Fold Change",
    color = "Significant (padj < 0.05)",
    title = "MA Plot"
  ) +
  theme_minimal()
```

### 4. Volcano Plot

Visualize statistical significance vs fold change magnitude:

```r
ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "gray60"),
    labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value",
    color = "Significant (padj < 0.05)",
    title = "Volcano Plot"
  ) +
  theme_minimal()
```

### 5. Using ggpubr Package (Simplified MA Plot)

The `ggpubr` package provides convenient functions:

```r
library(ggpubr)

ggmaplot(
  data,
  main = expression("Condition A" %->% "Condition B"),
  fdr = 0.05,
  fc = 2,
  size = 0.4,
  palette = c("#B31B21", "#1465AC", "darkgray"),
  genenames = as.vector(data$gene_id),  # Use gene_id column (adjust to match your data)
  legend = "top",
  top = 20,
  font.label = c("bold", 11),
  font.legend = "bold",
  font.main = "bold",
  ggtheme = theme_minimal()
)
```

**Parameters:**
- `fdr`: False discovery rate threshold (default 0.05)
- `fc`: Fold change threshold (default 2, meaning |log2FC| > 1)
- `palette`: Colors for up-regulated, down-regulated, and non-significant
- `top`: Number of top genes to label

## Data Preparation

### Merging Two Datasets for Comparison

Before creating comparison plots, merge your datasets:

```r
# Merge datasets by gene_id
merged_data <- merge(
  dataset1[, c("gene_id", "log2FoldChange", "padj", "baseMean")],
  dataset2[, c("gene_id", "log2FoldChange", "padj", "baseMean")],
  by = "gene_id",
  suffixes = c("_ds1", "_ds2")
)

# Create convenience columns for plotting
merged_data$logFC_dataset1 <- merged_data$log2FoldChange_ds1
merged_data$logFC_dataset2 <- merged_data$log2FoldChange_ds2

# Create combined FDR (use maximum, or significant in either dataset)
merged_data$FDR <- pmax(merged_data$padj_ds1, merged_data$padj_ds2)
# Alternative: significant if significant in either dataset
# merged_data$FDR <- ifelse(merged_data$padj_ds1 < 0.05 | merged_data$padj_ds2 < 0.05, 
#                           pmin(merged_data$padj_ds1, merged_data$padj_ds2), 
#                           pmax(merged_data$padj_ds1, merged_data$padj_ds2))
```

## Usage Examples

### Example 1: Basic Dataset Comparison

```r
# After merging (see Data Preparation above)
# Create significance column (significant in either dataset)
merged_data$significant <- merged_data$padj_ds1 < 0.05 | merged_data$padj_ds2 < 0.05

# Plot comparison
ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray60")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Dataset 1 Log2FC", y = "Dataset 2 Log2FC") +
  theme_minimal()
```

### Example 2: Multi-Category Significance

Color code by significance status in each dataset:

```r
merged_data$sig_status <- ifelse(
  merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05, "Both significant",
  ifelse(merged_data$padj_ds1 < 0.05, "Dataset 1 only",
  ifelse(merged_data$padj_ds2 < 0.05, "Dataset 2 only", "Neither"))
)

ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
  geom_point(aes(color = sig_status), alpha = 0.6) +
  scale_color_manual(
    values = c(
      "Both significant" = "red",
      "Dataset 1 only" = "orange",
      "Dataset 2 only" = "blue",
      "Neither" = "gray60"
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Dataset 1 Log2FC", y = "Dataset 2 Log2FC", color = "Significance") +
  theme_minimal()
```

## Troubleshooting

### Common Issues

**Problem: No significant genes showing up**
- Check that your `padj` or `FDR` column contains values < 0.05
- Verify p-values were properly adjusted (use `p.adjust(pvalue, method="BH")` for Benjamini-Hochberg correction)
- Ensure your data has sufficient statistical power

**Problem: Labels overlapping or not showing**
- Adjust `max.overlaps` parameter in `geom_text_repel()` (default is 10)
- Reduce the number of genes to label by increasing the threshold (e.g., `abs(logFC) > 2` instead of `> 1`)
- Use `size` parameter to make labels smaller

**Problem: Colors not showing correctly**
- Ensure `scale_color_manual()` uses `values = c("TRUE" = "...", "FALSE" = "...")` format
- Check that your condition (e.g., `FDR < 0.05`) evaluates to logical TRUE/FALSE
- Verify column names match exactly (case-sensitive)

**Problem: Reference lines not appearing**
- Check that `geom_abline()`, `geom_hline()`, or `geom_vline()` are added after `geom_point()`
- Verify intercept values are within your data range
- For `geom_abline()`, ensure `slope` and `intercept` are numeric values

**Problem: ggpubr not working**
- Install with: `install.packages("ggpubr", repos = "https://cran.rstudio.com/")`
- Check that your data has required columns (see ggpubr documentation)
- Verify `genenames` parameter matches your gene identifier column name

### Performance Tips

- For large datasets (>10,000 genes), consider sampling or filtering before plotting
- Use `alpha` transparency (0.6-0.7) to visualize overlapping points
- Save high-resolution plots with `ggsave(..., dpi = 300)` for publications

## Related

- **ggpubr package**: https://rpkgs.datanovia.com/ggpubr/reference/ggmaplot.html
- **EnhancedVolcano package**: For advanced volcano plot customization
- **ggrepel package**: For non-overlapping text labels
- **DESeq2 documentation**: For generating log fold change and p-value data
- **Test script**: See `scripts/create_test_graphs.R` for a complete working example with synthetic data generation

```


Verifying the updates and running validation:



Checking for the old diary file and removing it if present:



Removing the old diary file without the numeric prefix:



Verifying the troubleshooting section was added:



## Updates Made

1. **Data Structure Requirements** — Added a section documenting required columns for individual and merged datasets, including notes on `padj` vs `FDR`.

2. **Data Preparation Section** — Added code for merging two datasets with:
   - Proper column naming conventions
   - FDR combination strategies
   - Convenience columns for plotting

3. **Clarifications**:
   - Added a note explaining that `abs(logFC) > 2` means a 4-fold change (not 2-fold)
   - Fixed ggpubr example to use `gene_id` instead of `name`
   - Clarified threshold meanings throughout

4. **Troubleshooting Section** — Added common issues and solutions:
   - No significant genes showing up
   - Labels overlapping
   - Colors not displaying correctly
   - Reference lines not appearing
   - ggpubr installation issues
   - Performance tips for large datasets

5. **Related Files** — Updated to link to:
   - The test script (`scripts/create_test_graphs.R`)
   - The implementation diary (`log/01-2025-12-04-test-datasets-and-graphs.md`)

6. **Usage Examples** — Updated to reference the data preparation section and clarify the workflow.

All docmgr validation checks pass. The reference document now reflects the implementation experience and includes corrections and improvements.

