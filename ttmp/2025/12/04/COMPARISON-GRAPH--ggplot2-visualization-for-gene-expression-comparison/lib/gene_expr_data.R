# Gene Expression Data Generation Library
# Functions for generating synthetic gene expression datasets

#' Generate a synthetic gene expression dataset
#'
#' @param n_genes Number of genes to generate
#' @param seed Random seed for reproducibility
#' @param logFC_mean Mean of log2 fold change distribution
#' @param logFC_sd Standard deviation of log2 fold change distribution
#' @param baseMean_mean Mean of base mean expression (log10 scale)
#' @param baseMean_sd Standard deviation of base mean expression (log10 scale)
#' @return A data.frame with columns: gene_id, log2FoldChange, baseMean, pvalue, padj, FDR
generate_gene_dataset <- function(n_genes = 1000, 
                                   seed = 42,
                                   logFC_mean = 0,
                                   logFC_sd = 1.5,
                                   baseMean_mean = 3,
                                   baseMean_sd = 1.5) {
  set.seed(seed)
  
  # Generate gene IDs
  gene_ids <- paste0("GENE_", sprintf("%04d", 1:n_genes))
  
  # Generate log2 fold changes (with some correlation structure)
  logFC <- rnorm(n_genes, mean = logFC_mean, sd = logFC_sd)
  # Add some genes with larger fold changes
  logFC[sample(1:n_genes, min(100, n_genes %/% 10))] <- rnorm(
    min(100, n_genes %/% 10), 
    mean = logFC_mean, 
    sd = logFC_sd * 2
  )
  
  # Generate base mean expression (log scale)
  baseMean <- 10^rnorm(n_genes, mean = baseMean_mean, sd = baseMean_sd)
  baseMean <- pmax(baseMean, 1)  # Ensure positive values
  
  # Generate p-values (smaller p-values for larger fold changes)
  large_fc_idx <- abs(logFC) > 1.5
  pvalue <- numeric(n_genes)
  pvalue[large_fc_idx] <- runif(sum(large_fc_idx), min = 0.0001, max = 0.01)
  pvalue[!large_fc_idx] <- runif(sum(!large_fc_idx), min = 0.01, max = 0.5)
  # Add some randomness
  pvalue <- pvalue * runif(n_genes, min = 0.5, max = 1.5)
  pvalue <- pmax(pvalue, 1e-10)  # Ensure minimum value
  
  # Generate adjusted p-values (padj)
  padj <- p.adjust(pvalue, method = "BH")
  
  # Create dataset
  data.frame(
    gene_id = gene_ids,
    log2FoldChange = logFC,
    baseMean = baseMean,
    pvalue = pvalue,
    padj = padj,
    FDR = padj  # FDR is same as padj for BH correction
  )
}

#' Generate a second dataset correlated with the first
#'
#' @param dataset1 First dataset (from generate_gene_dataset)
#' @param correlation Correlation coefficient between datasets (0-1)
#' @param seed Random seed for reproducibility
#' @return A data.frame with the same structure as dataset1
generate_correlated_dataset <- function(dataset1, 
                                         correlation = 0.7,
                                         seed = 42) {
  set.seed(seed)
  
  n_genes <- nrow(dataset1)
  
  # Create correlated log fold changes (with noise)
  logFC_ds1 <- dataset1$log2FoldChange
  logFC_ds2 <- correlation * logFC_ds1 + 
    sqrt(1 - correlation^2) * rnorm(n_genes, mean = 0, sd = 1.5)
  # Add some dataset-specific differences
  logFC_ds2[sample(1:n_genes, min(50, n_genes %/% 20))] <- rnorm(
    min(50, n_genes %/% 20), 
    mean = 0, 
    sd = 2
  )
  
  # Generate p-values for dataset 2
  large_fc_idx_ds2 <- abs(logFC_ds2) > 1.5
  pvalue_ds2 <- numeric(n_genes)
  pvalue_ds2[large_fc_idx_ds2] <- runif(sum(large_fc_idx_ds2), min = 0.0001, max = 0.01)
  pvalue_ds2[!large_fc_idx_ds2] <- runif(sum(!large_fc_idx_ds2), min = 0.01, max = 0.5)
  # Add some randomness
  pvalue_ds2 <- pvalue_ds2 * runif(n_genes, min = 0.5, max = 1.5)
  pvalue_ds2 <- pmax(pvalue_ds2, 1e-10)  # Ensure minimum value
  padj_ds2 <- p.adjust(pvalue_ds2, method = "BH")
  
  # Create dataset 2
  data.frame(
    gene_id = dataset1$gene_id,
    log2FoldChange = logFC_ds2,
    baseMean = dataset1$baseMean,  # Same base mean for comparison
    pvalue = pvalue_ds2,
    padj = padj_ds2,
    FDR = padj_ds2
  )
}

#' Merge two datasets for comparison
#'
#' @param dataset1 First dataset
#' @param dataset2 Second dataset
#' @return Merged data.frame with columns for both datasets
merge_datasets <- function(dataset1, dataset2) {
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
  
  # Add significance status
  merged_data$sig_status <- ifelse(
    merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05, "Both significant",
    ifelse(merged_data$padj_ds1 < 0.05, "Dataset 1 only",
    ifelse(merged_data$padj_ds2 < 0.05, "Dataset 2 only", "Neither"))
  )
  
  merged_data
}

#' Print summary statistics for datasets
#'
#' @param dataset1 First dataset
#' @param dataset2 Second dataset
#' @param merged_data Merged dataset
print_dataset_summary <- function(dataset1, dataset2, merged_data) {
  cat("Generated", nrow(dataset1), "genes\n")
  cat("Dataset 1: Significant genes (padj < 0.05):", sum(dataset1$padj < 0.05), "\n")
  cat("Dataset 2: Significant genes (padj < 0.05):", sum(dataset2$padj < 0.05), "\n")
  cat("Both datasets significant:", 
      sum(merged_data$padj_ds1 < 0.05 & merged_data$padj_ds2 < 0.05), "\n")
}

