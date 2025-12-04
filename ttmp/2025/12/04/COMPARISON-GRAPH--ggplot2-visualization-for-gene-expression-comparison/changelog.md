# Changelog

## 2025-12-04

- Initial workspace created


## 2025-12-04

Created reference document on ggplot2 for comparing gene expression datasets with significant/insignificant point visualization


## 2025-12-04

Created comprehensive diary documenting test dataset generation and graph creation process. Documented what worked, what didn't work, lessons learned, and future improvements.

### Related Files

- log/2025-12-04-test-datasets-and-graphs.md — Implementation diary with full session notes


## 2025-12-04

Renamed diary file to include numeric prefix (01-) to follow docmgr conventions

### Related Files

- log/01-2025-12-04-test-datasets-and-graphs.md — Renamed diary file with proper numeric prefix


## 2025-12-04

Updated reference document with data structure requirements, data preparation section, troubleshooting guide, and corrections based on implementation experience. Fixed ggpubr example to use gene_id, clarified threshold meanings, and added common pitfalls section.

### Related Files

- reference/01-ggplot2-for-comparing-gene-expression-datasets.md — Updated with corrections and improvements


## 2025-12-04

Successfully ran create_test_graphs.R script and generated 7 visualization plots: scatter plots comparing datasets, MA plots for each dataset, volcano plots for each dataset, and multi-category significance comparison plot. Script generated 1000 genes with 364 significant in dataset 1, 376 significant in dataset 2, and 206 significant in both.

### Related Files

- scripts/create_test_graphs.R — Main script that generates test datasets and creates all visualization plots


## 2025-12-04

Refactored code into library modules: lib/gene_expr_data.R for data generation functions and lib/gene_expr_plots.R for plotting functions. Created Shiny app (app.R) for interactive exploration with adjustable parameters. Updated create_test_graphs.R script to use library functions.

### Related Files

- app.R — Shiny application for interactive gene expression visualization
- lib/gene_expr_data.R — Data generation functions (generate_gene_dataset
- lib/gene_expr_plots.R — Plotting functions for all visualization types
- scripts/create_test_graphs.R — Refactored to use library functions


## 2025-12-04

Created comprehensive onboarding guide for interns covering project context, key concepts, resources, and next steps

