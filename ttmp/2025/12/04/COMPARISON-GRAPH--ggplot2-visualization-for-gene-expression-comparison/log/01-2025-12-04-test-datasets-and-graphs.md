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
    - Path: 2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/log/2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/log/reference/01-ggplot2-for-comparing-gene-expression-datasets.md
      Note: |-
        Reference document containing ggplot2 code examples for gene expression visualization
        Reference document containing ggplot2 code examples used as basis for script
    - Path: 2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/log/2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/log/scripts/create_test_graphs.R
      Note: |-
        R script that generates synthetic gene expression datasets and creates all visualization types from the reference document
        R script that generates synthetic gene expression datasets and creates all visualization types
    - Path: 2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/log/2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/log/various/plots
      Note: Directory containing all generated PNG plots and CSV test datasets
    - Path: 2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/log/reference/01-ggplot2-for-comparing-gene-expression-datasets.md
      Note: Reference document containing ggplot2 code examples used as basis for script
    - Path: 2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/log/scripts/create_test_graphs.R
      Note: R script that generates synthetic gene expression datasets and creates all visualization types
    - Path: 2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/log/various/plots
      Note: Directory containing all generated PNG plots and CSV test datasets
ExternalSources: []
Summary: ""
LastUpdated: 2025-12-04T08:05:00-05:00
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

