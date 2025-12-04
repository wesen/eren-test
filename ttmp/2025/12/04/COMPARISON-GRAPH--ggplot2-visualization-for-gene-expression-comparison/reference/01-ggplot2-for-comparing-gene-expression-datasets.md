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
RelatedFiles:
    - Path: 2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/reference/app.R
      Note: Interactive Shiny application for exploring gene expression visualizations
    - Path: 2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/reference/lib/gene_expr_data.R
      Note: Refactored data generation functions library
    - Path: 2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/reference/lib/gene_expr_plots.R
      Note: Refactored plotting functions library
    - Path: 2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/scripts/create_test_graphs.R
      Note: Implementation script that generates test gene expression datasets and creates ggplot2 visualizations
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
