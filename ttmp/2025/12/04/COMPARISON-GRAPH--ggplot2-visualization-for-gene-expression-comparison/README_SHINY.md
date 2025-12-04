# Shiny App Usage Guide

## Running the Shiny App

### Method 1: Using Rscript
```bash
cd /tmp/eren-test/ttmp/2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison
Rscript app.R
```

### Method 2: Using R console
```r
library(shiny)
runApp("app.R")
```

### Method 3: Using shiny::runApp()
```r
shiny::runApp("app.R", port = 3838, host = "0.0.0.0")
```

## Features

The Shiny app provides interactive exploration of gene expression data with:

### Data Generation Parameters
- **Number of genes**: Adjust the dataset size (100-10,000)
- **Dataset correlation**: Control correlation between datasets (0-1)
- **Random seed**: Set seed for reproducibility

### Plot Parameters
- **FDR threshold**: Adjust significance threshold (default: 0.05)
- **Fold change threshold**: Set fold change cutoff for volcano plots and labeling
- **Point size**: Adjust visualization point size
- **Point transparency**: Control point alpha/transparency

### Available Visualizations
1. **Scatter Plot (Basic)**: Compare log2 fold changes between datasets
2. **Scatter Plot (With Labels)**: Enhanced scatter with gene labels
3. **MA Plot - Dataset 1**: Mean expression vs. fold change for dataset 1
4. **MA Plot - Dataset 2**: Mean expression vs. fold change for dataset 2
5. **Volcano Plot - Dataset 1**: Fold change vs. significance for dataset 1
6. **Volcano Plot - Dataset 2**: Fold change vs. significance for dataset 2
7. **Multi-Category Scatter**: Color-coded by significance status
8. **Data Tables**: Browse the generated datasets

## Library Functions

The app uses refactored library functions:

- `lib/gene_expr_data.R`: Data generation functions
  - `generate_gene_dataset()`: Create synthetic gene expression data
  - `generate_correlated_dataset()`: Create correlated second dataset
  - `merge_datasets()`: Merge datasets for comparison

- `lib/gene_expr_plots.R`: Plotting functions
  - `plot_scatter_basic()`: Basic scatter plot
  - `plot_scatter_with_labels()`: Scatter with gene labels
  - `plot_ma()`: MA plot
  - `plot_volcano()`: Volcano plot
  - `plot_scatter_multi_category()`: Multi-category scatter

## Dependencies

Required R packages:
- `shiny` (for the app)
- `ggplot2` (for plotting)
- `ggrepel` (for label positioning)
- `DT` (optional, for enhanced data tables)

Install missing packages:
```r
install.packages(c("shiny", "ggplot2", "ggrepel", "DT"))
```

