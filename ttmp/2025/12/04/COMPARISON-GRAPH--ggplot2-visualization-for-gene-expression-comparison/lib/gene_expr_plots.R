# Gene Expression Plotting Library
# Functions for creating ggplot2 visualizations

library(ggplot2)
library(ggrepel)

#' Create a basic scatter plot comparing two datasets
#'
#' @param merged_data Merged dataset from merge_datasets()
#' @param fdr_threshold FDR threshold for significance (default: 0.05)
#' @param point_size Point size
#' @param point_alpha Point transparency
#' @return ggplot object
plot_scatter_basic <- function(merged_data, 
                                fdr_threshold = 0.05,
                                point_size = 1.2,
                                point_alpha = 0.7) {
  ggplot(merged_data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
    geom_point(aes(color = FDR < fdr_threshold), alpha = point_alpha, size = point_size) +
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "gray60"),
      labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
      name = "Significance"
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
    labs(
      x = "Log2 Fold Change (Dataset 1)",
      y = "Log2 Fold Change (Dataset 2)",
      title = "Comparison of Log2 Fold Changes Between Two Datasets"
    ) +
    theme_minimal()
}

#' Create a scatter plot with gene labels
#'
#' @param merged_data Merged dataset from merge_datasets()
#' @param fdr_threshold FDR threshold for significance
#' @param fc_threshold Fold change threshold for labeling
#' @param max_labels Maximum number of labels to show
#' @param point_size Point size
#' @param point_alpha Point transparency
#' @return ggplot object
plot_scatter_with_labels <- function(merged_data,
                                     fdr_threshold = 0.05,
                                     fc_threshold = 2,
                                     max_labels = 20,
                                     point_size = 1.2,
                                     point_alpha = 0.7) {
  ggplot(merged_data, aes(x = logFC_dataset1, y = logFC_dataset2)) +
    geom_point(aes(color = FDR < fdr_threshold), alpha = point_alpha, size = point_size) +
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "gray60"),
      labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
      name = paste("Significant (FDR <", fdr_threshold, ")")
    ) +
    geom_text_repel(
      data = subset(merged_data, 
                    FDR < fdr_threshold & 
                    abs(logFC_dataset1) > fc_threshold & 
                    abs(logFC_dataset2) > fc_threshold),
      aes(label = gene_id),
      size = 3,
      max.overlaps = max_labels
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
    labs(
      x = "Log2 Fold Change (Dataset 1)",
      y = "Log2 Fold Change (Dataset 2)"
    ) +
    theme_minimal()
}

#' Create an MA plot
#'
#' @param dataset Dataset from generate_gene_dataset()
#' @param padj_threshold Adjusted p-value threshold for significance
#' @param point_size Point size
#' @param point_alpha Point transparency
#' @param title Plot title
#' @return ggplot object
plot_ma <- function(dataset,
                    padj_threshold = 0.05,
                    point_size = 1.2,
                    point_alpha = 0.6,
                    title = "MA Plot") {
  ggplot(dataset, aes(x = log2(baseMean), y = log2FoldChange)) +
    geom_point(aes(color = padj < padj_threshold), alpha = point_alpha, size = point_size) +
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "gray60"),
      labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
      name = paste("Significant (padj <", padj_threshold, ")")
    ) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    labs(
      x = "Log2 Mean Expression",
      y = "Log2 Fold Change",
      title = title
    ) +
    theme_minimal()
}

#' Create a volcano plot
#'
#' @param dataset Dataset from generate_gene_dataset()
#' @param padj_threshold Adjusted p-value threshold for significance
#' @param fc_threshold Fold change threshold (for vertical lines)
#' @param point_size Point size
#' @param point_alpha Point transparency
#' @param title Plot title
#' @return ggplot object
plot_volcano <- function(dataset,
                         padj_threshold = 0.05,
                         fc_threshold = 1,
                         point_size = 1.2,
                         point_alpha = 0.6,
                         title = "Volcano Plot") {
  ggplot(dataset, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = padj < padj_threshold), alpha = point_alpha, size = point_size) +
    scale_color_manual(
      values = c("TRUE" = "red", "FALSE" = "gray60"),
      labels = c("TRUE" = "Significant", "FALSE" = "Not significant"),
      name = paste("Significant (padj <", padj_threshold, ")")
    ) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), 
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(padj_threshold), 
               linetype = "dashed", color = "black") +
    labs(
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      title = title
    ) +
    theme_minimal()
}

#' Create a multi-category significance scatter plot
#'
#' @param merged_data Merged dataset from merge_datasets()
#' @param padj_threshold Adjusted p-value threshold for significance
#' @param point_size Point size
#' @param point_alpha Point transparency
#' @param title Plot title
#' @return ggplot object
plot_scatter_multi_category <- function(merged_data,
                                        padj_threshold = 0.05,
                                        point_size = 1.2,
                                        point_alpha = 0.6,
                                        title = "Multi-Category Significance Comparison") {
  ggplot(merged_data, aes(x = log2FoldChange_ds1, y = log2FoldChange_ds2)) +
    geom_point(aes(color = sig_status), alpha = point_alpha, size = point_size) +
    scale_color_manual(
      values = c(
        "Both significant" = "red",
        "Dataset 1 only" = "orange",
        "Dataset 2 only" = "blue",
        "Neither" = "gray60"
      ),
      name = "Significance"
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(
      x = "Dataset 1 Log2FC",
      y = "Dataset 2 Log2FC",
      title = title
    ) +
    theme_minimal()
}

