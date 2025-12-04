# Read documentation and generate comparison graph
_Exported on 12/4/2025 at 08:38:19 EST from Cursor (2.1.39)_

---

**User**

Read `docmgr help how-to-use` (redirect to file, and then read the file in a second step in order to have it in full). Make sure you read the file in full full full so you get all the infromation (it's around 1200 lines).

Then work on ticket COMPARISON-GRAPH.



Run the script and show us the graph.

---

**Cursor**


Reading `docmgr help how-to-use`, then working on COMPARISON-GRAPH and running the script to show the graph.



Reading the docmgr help file:



Reading the docmgr help file in sections:



Running the script to generate the graphs:



Script ran successfully. Listing generated plots:



All plots generated. Checking ticket status and updating documentation per docmgr workflow:



Updating the changelog and relating the script file per docmgr workflow:



Finding the reference document path and relating the script correctly:



Relating the script to the reference document:



Script ran successfully. Summary:

## Script execution results

**Generated datasets:**
- 1,000 genes total
- Dataset 1: 364 significant genes (padj < 0.05)
- Dataset 2: 376 significant genes (padj < 0.05)
- Both datasets: 206 genes significant in both

**Generated visualizations (7 plots):**

1. `01_scatter_basic.png` — Basic scatter plot comparing log2 fold changes between datasets, colored by significance
2. `02_scatter_with_labels.png` — Enhanced scatter plot with gene labels for highly significant genes
3. `03_ma_plot_dataset1.png` — MA plot for Dataset 1 (mean expression vs. fold change)
4. `04_ma_plot_dataset2.png` — MA plot for Dataset 2
5. `05_volcano_plot_dataset1.png` — Volcano plot for Dataset 1 (fold change vs. significance)
6. `06_volcano_plot_dataset2.png` — Volcano plot for Dataset 2
7. `07_scatter_multi_category.png` — Multi-category scatter plot showing significance status (both/neither/dataset-specific)

All plots are saved in: `various/plots/`  
Test datasets saved as CSV: `dataset1.csv`, `dataset2.csv`, `merged_data.csv`

The plots are PNG images (3000x2400 pixels, 300 DPI) ready for publication or presentation. The script successfully demonstrates the visualization patterns from your reference document.

I've updated the changelog and related the script file to the reference document per docmgr workflow. The ticket is ready for review.

---

**User**

Create a shiny app to allow me to play around. Refactor into a lib in the main directory as well.



Create tasks with docmgr, keep changelog.

---

**Cursor**


Creating a Shiny app and refactoring the code into a library. Reviewing the current script structure, then setting up the library and app.



Creating tasks, then building the library and Shiny app:



Creating the library structure and Shiny app:


```r
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

```

```r
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

```


Creating the Shiny app:


```r
#!/usr/bin/env Rscript
# Shiny App for Interactive Gene Expression Visualization
# Run with: Rscript app.R
# Or: shiny::runApp()

library(shiny)
library(ggplot2)
library(ggrepel)

# Source library functions
source("lib/gene_expr_data.R")
source("lib/gene_expr_plots.R")

# UI
ui <- fluidPage(
  titlePanel("Gene Expression Comparison Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Data Generation Parameters"),
      
      numericInput("n_genes", 
                   "Number of genes:", 
                   value = 1000, 
                   min = 100, 
                   max = 10000, 
                   step = 100),
      
      numericInput("correlation", 
                   "Dataset correlation:", 
                   value = 0.7, 
                   min = 0, 
                   max = 1, 
                   step = 0.1),
      
      numericInput("seed", 
                   "Random seed:", 
                   value = 42, 
                   min = 1),
      
      br(),
      h3("Plot Parameters"),
      
      numericInput("fdr_threshold", 
                   "FDR threshold:", 
                   value = 0.05, 
                   min = 0.001, 
                   max = 0.5, 
                   step = 0.01),
      
      numericInput("fc_threshold", 
                   "Fold change threshold:", 
                   value = 1, 
                   min = 0.1, 
                   max = 5, 
                   step = 0.1),
      
      numericInput("point_size", 
                   "Point size:", 
                   value = 1.2, 
                   min = 0.5, 
                   max = 5, 
                   step = 0.1),
      
      numericInput("point_alpha", 
                   "Point transparency:", 
                   value = 0.7, 
                   min = 0.1, 
                   max = 1, 
                   step = 0.1),
      
      br(),
      actionButton("generate_data", "Generate New Data", class = "btn-primary"),
      
      br(), br(),
      h4("Dataset Summary"),
      verbatimTextOutput("summary")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Scatter Plot (Basic)",
                 plotOutput("scatter_basic", height = "600px")),
        
        tabPanel("Scatter Plot (With Labels)",
                 numericInput("max_labels", 
                              "Max labels:", 
                              value = 20, 
                              min = 5, 
                              max = 100, 
                              step = 5,
                              width = "150px"),
                 plotOutput("scatter_labels", height = "600px")),
        
        tabPanel("MA Plot - Dataset 1",
                 plotOutput("ma_plot1", height = "600px")),
        
        tabPanel("MA Plot - Dataset 2",
                 plotOutput("ma_plot2", height = "600px")),
        
        tabPanel("Volcano Plot - Dataset 1",
                 plotOutput("volcano_plot1", height = "600px")),
        
        tabPanel("Volcano Plot - Dataset 2",
                 plotOutput("volcano_plot2", height = "600px")),
        
        tabPanel("Multi-Category Scatter",
                 plotOutput("scatter_multi", height = "600px")),
        
        tabPanel("Data Tables",
                 tabsetPanel(
                   tabPanel("Dataset 1",
                            DT::dataTableOutput("table1")),
                   tabPanel("Dataset 2",
                            DT::dataTableOutput("table2")),
                   tabPanel("Merged",
                            DT::dataTableOutput("table_merged"))
                 ))
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values to store datasets
  values <- reactiveValues(
    dataset1 = NULL,
    dataset2 = NULL,
    merged_data = NULL
  )
  
  # Generate initial data
  observe({
    dataset1 <- generate_gene_dataset(
      n_genes = input$n_genes,
      seed = input$seed
    )
    dataset2 <- generate_correlated_dataset(
      dataset1,
      correlation = input$correlation,
      seed = input$seed
    )
    merged_data <- merge_datasets(dataset1, dataset2)
    
    values$dataset1 <- dataset1
    values$dataset2 <- dataset2
    values$merged_data <- merged_data
  })
  
  # Regenerate data when button is clicked
  observeEvent(input$generate_data, {
    dataset1 <- generate_gene_dataset(
      n_genes = input$n_genes,
      seed = input$seed
    )
    dataset2 <- generate_correlated_dataset(
      dataset1,
      correlation = input$correlation,
      seed = input$seed
    )
    merged_data <- merge_datasets(dataset1, dataset2)
    
    values$dataset1 <- dataset1
    values$dataset2 <- dataset2
    values$merged_data <- merged_data
  })
  
  # Summary output
  output$summary <- renderText({
    if (!is.null(values$merged_data)) {
      paste(
        "Genes:", nrow(values$dataset1), "\n",
        "Dataset 1 significant:", sum(values$dataset1$padj < input$fdr_threshold), "\n",
        "Dataset 2 significant:", sum(values$dataset2$padj < input$fdr_threshold), "\n",
        "Both significant:", 
        sum(values$merged_data$padj_ds1 < input$fdr_threshold & 
            values$merged_data$padj_ds2 < input$fdr_threshold)
      )
    }
  })
  
  # Plots
  output$scatter_basic <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_basic(
        values$merged_data,
        fdr_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  output$scatter_labels <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_with_labels(
        values$merged_data,
        fdr_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        max_labels = input$max_labels,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  output$ma_plot1 <- renderPlot({
    if (!is.null(values$dataset1)) {
      plot_ma(
        values$dataset1,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "MA Plot - Dataset 1"
      )
    }
  })
  
  output$ma_plot2 <- renderPlot({
    if (!is.null(values$dataset2)) {
      plot_ma(
        values$dataset2,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "MA Plot - Dataset 2"
      )
    }
  })
  
  output$volcano_plot1 <- renderPlot({
    if (!is.null(values$dataset1)) {
      plot_volcano(
        values$dataset1,
        padj_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "Volcano Plot - Dataset 1"
      )
    }
  })
  
  output$volcano_plot2 <- renderPlot({
    if (!is.null(values$dataset2)) {
      plot_volcano(
        values$dataset2,
        padj_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "Volcano Plot - Dataset 2"
      )
    }
  })
  
  output$scatter_multi <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_multi_category(
        values$merged_data,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  # Data tables
  output$table1 <- DT::renderDataTable({
    if (!is.null(values$dataset1)) {
      DT::datatable(values$dataset1, options = list(pageLength = 25))
    }
  })
  
  output$table2 <- DT::renderDataTable({
    if (!is.null(values$dataset2)) {
      DT::datatable(values$dataset2, options = list(pageLength = 25))
    }
  })
  
  output$table_merged <- DT::renderDataTable({
    if (!is.null(values$merged_data)) {
      DT::datatable(values$merged_data, options = list(pageLength = 25))
    }
  })
}

# Run the application
if (!interactive()) {
  shinyApp(ui = ui, server = server, options = list(port = 3838, host = "0.0.0.0"))
} else {
  shinyApp(ui = ui, server = server)
}

```


Updating the script to use the library functions:


```r
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

```


Fixing the Shiny app to handle the DT package properly:


```r
#!/usr/bin/env Rscript
# Shiny App for Interactive Gene Expression Visualization
# Run with: Rscript app.R
# Or: shiny::runApp()

library(shiny)
library(ggplot2)
library(ggrepel)

# Load DT if available, otherwise use basic tables
if (requireNamespace("DT", quietly = TRUE)) {
  library(DT)
  use_dt <- TRUE
} else {
  use_dt <- FALSE
}

# Source library functions
source("lib/gene_expr_data.R")
source("lib/gene_expr_plots.R")

# UI
ui <- fluidPage(
  titlePanel("Gene Expression Comparison Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Data Generation Parameters"),
      
      numericInput("n_genes", 
                   "Number of genes:", 
                   value = 1000, 
                   min = 100, 
                   max = 10000, 
                   step = 100),
      
      numericInput("correlation", 
                   "Dataset correlation:", 
                   value = 0.7, 
                   min = 0, 
                   max = 1, 
                   step = 0.1),
      
      numericInput("seed", 
                   "Random seed:", 
                   value = 42, 
                   min = 1),
      
      br(),
      h3("Plot Parameters"),
      
      numericInput("fdr_threshold", 
                   "FDR threshold:", 
                   value = 0.05, 
                   min = 0.001, 
                   max = 0.5, 
                   step = 0.01),
      
      numericInput("fc_threshold", 
                   "Fold change threshold:", 
                   value = 1, 
                   min = 0.1, 
                   max = 5, 
                   step = 0.1),
      
      numericInput("point_size", 
                   "Point size:", 
                   value = 1.2, 
                   min = 0.5, 
                   max = 5, 
                   step = 0.1),
      
      numericInput("point_alpha", 
                   "Point transparency:", 
                   value = 0.7, 
                   min = 0.1, 
                   max = 1, 
                   step = 0.1),
      
      br(),
      actionButton("generate_data", "Generate New Data", class = "btn-primary"),
      
      br(), br(),
      h4("Dataset Summary"),
      verbatimTextOutput("summary")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Scatter Plot (Basic)",
                 plotOutput("scatter_basic", height = "600px")),
        
        tabPanel("Scatter Plot (With Labels)",
                 numericInput("max_labels", 
                              "Max labels:", 
                              value = 20, 
                              min = 5, 
                              max = 100, 
                              step = 5,
                              width = "150px"),
                 plotOutput("scatter_labels", height = "600px")),
        
        tabPanel("MA Plot - Dataset 1",
                 plotOutput("ma_plot1", height = "600px")),
        
        tabPanel("MA Plot - Dataset 2",
                 plotOutput("ma_plot2", height = "600px")),
        
        tabPanel("Volcano Plot - Dataset 1",
                 plotOutput("volcano_plot1", height = "600px")),
        
        tabPanel("Volcano Plot - Dataset 2",
                 plotOutput("volcano_plot2", height = "600px")),
        
        tabPanel("Multi-Category Scatter",
                 plotOutput("scatter_multi", height = "600px")),
        
        tabPanel("Data Tables",
                 tabsetPanel(
                   tabPanel("Dataset 1",
                            DT::dataTableOutput("table1")),
                   tabPanel("Dataset 2",
                            DT::dataTableOutput("table2")),
                   tabPanel("Merged",
                            DT::dataTableOutput("table_merged"))
                 ))
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values to store datasets
  values <- reactiveValues(
    dataset1 = NULL,
    dataset2 = NULL,
    merged_data = NULL
  )
  
  # Generate initial data
  observe({
    dataset1 <- generate_gene_dataset(
      n_genes = input$n_genes,
      seed = input$seed
    )
    dataset2 <- generate_correlated_dataset(
      dataset1,
      correlation = input$correlation,
      seed = input$seed
    )
    merged_data <- merge_datasets(dataset1, dataset2)
    
    values$dataset1 <- dataset1
    values$dataset2 <- dataset2
    values$merged_data <- merged_data
  })
  
  # Regenerate data when button is clicked
  observeEvent(input$generate_data, {
    dataset1 <- generate_gene_dataset(
      n_genes = input$n_genes,
      seed = input$seed
    )
    dataset2 <- generate_correlated_dataset(
      dataset1,
      correlation = input$correlation,
      seed = input$seed
    )
    merged_data <- merge_datasets(dataset1, dataset2)
    
    values$dataset1 <- dataset1
    values$dataset2 <- dataset2
    values$merged_data <- merged_data
  })
  
  # Summary output
  output$summary <- renderText({
    if (!is.null(values$merged_data)) {
      paste(
        "Genes:", nrow(values$dataset1), "\n",
        "Dataset 1 significant:", sum(values$dataset1$padj < input$fdr_threshold), "\n",
        "Dataset 2 significant:", sum(values$dataset2$padj < input$fdr_threshold), "\n",
        "Both significant:", 
        sum(values$merged_data$padj_ds1 < input$fdr_threshold & 
            values$merged_data$padj_ds2 < input$fdr_threshold)
      )
    }
  })
  
  # Plots
  output$scatter_basic <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_basic(
        values$merged_data,
        fdr_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  output$scatter_labels <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_with_labels(
        values$merged_data,
        fdr_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        max_labels = input$max_labels,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  output$ma_plot1 <- renderPlot({
    if (!is.null(values$dataset1)) {
      plot_ma(
        values$dataset1,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "MA Plot - Dataset 1"
      )
    }
  })
  
  output$ma_plot2 <- renderPlot({
    if (!is.null(values$dataset2)) {
      plot_ma(
        values$dataset2,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "MA Plot - Dataset 2"
      )
    }
  })
  
  output$volcano_plot1 <- renderPlot({
    if (!is.null(values$dataset1)) {
      plot_volcano(
        values$dataset1,
        padj_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "Volcano Plot - Dataset 1"
      )
    }
  })
  
  output$volcano_plot2 <- renderPlot({
    if (!is.null(values$dataset2)) {
      plot_volcano(
        values$dataset2,
        padj_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "Volcano Plot - Dataset 2"
      )
    }
  })
  
  output$scatter_multi <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_multi_category(
        values$merged_data,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  # Data tables
  output$table1 <- DT::renderDataTable({
    if (!is.null(values$dataset1)) {
      DT::datatable(values$dataset1, options = list(pageLength = 25))
    }
  })
  
  output$table2 <- DT::renderDataTable({
    if (!is.null(values$dataset2)) {
      DT::datatable(values$dataset2, options = list(pageLength = 25))
    }
  })
  
  output$table_merged <- DT::renderDataTable({
    if (!is.null(values$merged_data)) {
      DT::datatable(values$merged_data, options = list(pageLength = 25))
    }
  })
}

# Run the application
if (!interactive()) {
  shinyApp(ui = ui, server = server, options = list(port = 3838, host = "0.0.0.0"))
} else {
  shinyApp(ui = ui, server = server)
}


```

```r
#!/usr/bin/env Rscript
# Shiny App for Interactive Gene Expression Visualization
# Run with: Rscript app.R
# Or: shiny::runApp()

library(shiny)
library(ggplot2)
library(ggrepel)

# Load DT if available, otherwise use basic tables
if (requireNamespace("DT", quietly = TRUE)) {
  library(DT)
  use_dt <- TRUE
} else {
  use_dt <- FALSE
}

# Source library functions
source("lib/gene_expr_data.R")
source("lib/gene_expr_plots.R")

# UI
ui <- fluidPage(
  titlePanel("Gene Expression Comparison Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Data Generation Parameters"),
      
      numericInput("n_genes", 
                   "Number of genes:", 
                   value = 1000, 
                   min = 100, 
                   max = 10000, 
                   step = 100),
      
      numericInput("correlation", 
                   "Dataset correlation:", 
                   value = 0.7, 
                   min = 0, 
                   max = 1, 
                   step = 0.1),
      
      numericInput("seed", 
                   "Random seed:", 
                   value = 42, 
                   min = 1),
      
      br(),
      h3("Plot Parameters"),
      
      numericInput("fdr_threshold", 
                   "FDR threshold:", 
                   value = 0.05, 
                   min = 0.001, 
                   max = 0.5, 
                   step = 0.01),
      
      numericInput("fc_threshold", 
                   "Fold change threshold:", 
                   value = 1, 
                   min = 0.1, 
                   max = 5, 
                   step = 0.1),
      
      numericInput("point_size", 
                   "Point size:", 
                   value = 1.2, 
                   min = 0.5, 
                   max = 5, 
                   step = 0.1),
      
      numericInput("point_alpha", 
                   "Point transparency:", 
                   value = 0.7, 
                   min = 0.1, 
                   max = 1, 
                   step = 0.1),
      
      br(),
      actionButton("generate_data", "Generate New Data", class = "btn-primary"),
      
      br(), br(),
      h4("Dataset Summary"),
      verbatimTextOutput("summary")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Scatter Plot (Basic)",
                 plotOutput("scatter_basic", height = "600px")),
        
        tabPanel("Scatter Plot (With Labels)",
                 numericInput("max_labels", 
                              "Max labels:", 
                              value = 20, 
                              min = 5, 
                              max = 100, 
                              step = 5,
                              width = "150px"),
                 plotOutput("scatter_labels", height = "600px")),
        
        tabPanel("MA Plot - Dataset 1",
                 plotOutput("ma_plot1", height = "600px")),
        
        tabPanel("MA Plot - Dataset 2",
                 plotOutput("ma_plot2", height = "600px")),
        
        tabPanel("Volcano Plot - Dataset 1",
                 plotOutput("volcano_plot1", height = "600px")),
        
        tabPanel("Volcano Plot - Dataset 2",
                 plotOutput("volcano_plot2", height = "600px")),
        
        tabPanel("Multi-Category Scatter",
                 plotOutput("scatter_multi", height = "600px")),
        
        tabPanel("Data Tables",
                 tabsetPanel(
                   tabPanel("Dataset 1",
                            conditionalPanel(
                              condition = "output.use_dt",
                              DT::dataTableOutput("table1")
                            ),
                            conditionalPanel(
                              condition = "!output.use_dt",
                              tableOutput("table1_simple")
                            )),
                   tabPanel("Dataset 2",
                            conditionalPanel(
                              condition = "output.use_dt",
                              DT::dataTableOutput("table2")
                            ),
                            conditionalPanel(
                              condition = "!output.use_dt",
                              tableOutput("table2_simple")
                            )),
                   tabPanel("Merged",
                            conditionalPanel(
                              condition = "output.use_dt",
                              DT::dataTableOutput("table_merged")
                            ),
                            conditionalPanel(
                              condition = "!output.use_dt",
                              tableOutput("table_merged_simple")
                            ))
                 ))
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values to store datasets
  values <- reactiveValues(
    dataset1 = NULL,
    dataset2 = NULL,
    merged_data = NULL
  )
  
  # Generate initial data
  observe({
    dataset1 <- generate_gene_dataset(
      n_genes = input$n_genes,
      seed = input$seed
    )
    dataset2 <- generate_correlated_dataset(
      dataset1,
      correlation = input$correlation,
      seed = input$seed
    )
    merged_data <- merge_datasets(dataset1, dataset2)
    
    values$dataset1 <- dataset1
    values$dataset2 <- dataset2
    values$merged_data <- merged_data
  })
  
  # Regenerate data when button is clicked
  observeEvent(input$generate_data, {
    dataset1 <- generate_gene_dataset(
      n_genes = input$n_genes,
      seed = input$seed
    )
    dataset2 <- generate_correlated_dataset(
      dataset1,
      correlation = input$correlation,
      seed = input$seed
    )
    merged_data <- merge_datasets(dataset1, dataset2)
    
    values$dataset1 <- dataset1
    values$dataset2 <- dataset2
    values$merged_data <- merged_data
  })
  
  # Summary output
  output$summary <- renderText({
    if (!is.null(values$merged_data)) {
      paste(
        "Genes:", nrow(values$dataset1), "\n",
        "Dataset 1 significant:", sum(values$dataset1$padj < input$fdr_threshold), "\n",
        "Dataset 2 significant:", sum(values$dataset2$padj < input$fdr_threshold), "\n",
        "Both significant:", 
        sum(values$merged_data$padj_ds1 < input$fdr_threshold & 
            values$merged_data$padj_ds2 < input$fdr_threshold)
      )
    }
  })
  
  # Plots
  output$scatter_basic <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_basic(
        values$merged_data,
        fdr_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  output$scatter_labels <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_with_labels(
        values$merged_data,
        fdr_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        max_labels = input$max_labels,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  output$ma_plot1 <- renderPlot({
    if (!is.null(values$dataset1)) {
      plot_ma(
        values$dataset1,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "MA Plot - Dataset 1"
      )
    }
  })
  
  output$ma_plot2 <- renderPlot({
    if (!is.null(values$dataset2)) {
      plot_ma(
        values$dataset2,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "MA Plot - Dataset 2"
      )
    }
  })
  
  output$volcano_plot1 <- renderPlot({
    if (!is.null(values$dataset1)) {
      plot_volcano(
        values$dataset1,
        padj_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "Volcano Plot - Dataset 1"
      )
    }
  })
  
  output$volcano_plot2 <- renderPlot({
    if (!is.null(values$dataset2)) {
      plot_volcano(
        values$dataset2,
        padj_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "Volcano Plot - Dataset 2"
      )
    }
  })
  
  output$scatter_multi <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_multi_category(
        values$merged_data,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  # Data tables
  output$table1 <- DT::renderDataTable({
    if (!is.null(values$dataset1)) {
      DT::datatable(values$dataset1, options = list(pageLength = 25))
    }
  })
  
  output$table2 <- DT::renderDataTable({
    if (!is.null(values$dataset2)) {
      DT::datatable(values$dataset2, options = list(pageLength = 25))
    }
  })
  
  output$table_merged <- DT::renderDataTable({
    if (!is.null(values$merged_data)) {
      DT::datatable(values$merged_data, options = list(pageLength = 25))
    }
  })
}

# Run the application
if (!interactive()) {
  shinyApp(ui = ui, server = server, options = list(port = 3838, host = "0.0.0.0"))
} else {
  shinyApp(ui = ui, server = server)
}


```

```r
#!/usr/bin/env Rscript
# Shiny App for Interactive Gene Expression Visualization
# Run with: Rscript app.R
# Or: shiny::runApp()

library(shiny)
library(ggplot2)
library(ggrepel)

# Load DT if available, otherwise use basic tables
if (requireNamespace("DT", quietly = TRUE)) {
  library(DT)
  use_dt <- TRUE
} else {
  use_dt <- FALSE
}

# Source library functions
source("lib/gene_expr_data.R")
source("lib/gene_expr_plots.R")

# UI
ui <- fluidPage(
  titlePanel("Gene Expression Comparison Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Data Generation Parameters"),
      
      numericInput("n_genes", 
                   "Number of genes:", 
                   value = 1000, 
                   min = 100, 
                   max = 10000, 
                   step = 100),
      
      numericInput("correlation", 
                   "Dataset correlation:", 
                   value = 0.7, 
                   min = 0, 
                   max = 1, 
                   step = 0.1),
      
      numericInput("seed", 
                   "Random seed:", 
                   value = 42, 
                   min = 1),
      
      br(),
      h3("Plot Parameters"),
      
      numericInput("fdr_threshold", 
                   "FDR threshold:", 
                   value = 0.05, 
                   min = 0.001, 
                   max = 0.5, 
                   step = 0.01),
      
      numericInput("fc_threshold", 
                   "Fold change threshold:", 
                   value = 1, 
                   min = 0.1, 
                   max = 5, 
                   step = 0.1),
      
      numericInput("point_size", 
                   "Point size:", 
                   value = 1.2, 
                   min = 0.5, 
                   max = 5, 
                   step = 0.1),
      
      numericInput("point_alpha", 
                   "Point transparency:", 
                   value = 0.7, 
                   min = 0.1, 
                   max = 1, 
                   step = 0.1),
      
      br(),
      actionButton("generate_data", "Generate New Data", class = "btn-primary"),
      
      br(), br(),
      h4("Dataset Summary"),
      verbatimTextOutput("summary")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Scatter Plot (Basic)",
                 plotOutput("scatter_basic", height = "600px")),
        
        tabPanel("Scatter Plot (With Labels)",
                 numericInput("max_labels", 
                              "Max labels:", 
                              value = 20, 
                              min = 5, 
                              max = 100, 
                              step = 5,
                              width = "150px"),
                 plotOutput("scatter_labels", height = "600px")),
        
        tabPanel("MA Plot - Dataset 1",
                 plotOutput("ma_plot1", height = "600px")),
        
        tabPanel("MA Plot - Dataset 2",
                 plotOutput("ma_plot2", height = "600px")),
        
        tabPanel("Volcano Plot - Dataset 1",
                 plotOutput("volcano_plot1", height = "600px")),
        
        tabPanel("Volcano Plot - Dataset 2",
                 plotOutput("volcano_plot2", height = "600px")),
        
        tabPanel("Multi-Category Scatter",
                 plotOutput("scatter_multi", height = "600px")),
        
        tabPanel("Data Tables",
                 tabsetPanel(
                   tabPanel("Dataset 1",
                            conditionalPanel(
                              condition = "output.use_dt",
                              DT::dataTableOutput("table1")
                            ),
                            conditionalPanel(
                              condition = "!output.use_dt",
                              tableOutput("table1_simple")
                            )),
                   tabPanel("Dataset 2",
                            conditionalPanel(
                              condition = "output.use_dt",
                              DT::dataTableOutput("table2")
                            ),
                            conditionalPanel(
                              condition = "!output.use_dt",
                              tableOutput("table2_simple")
                            )),
                   tabPanel("Merged",
                            conditionalPanel(
                              condition = "output.use_dt",
                              DT::dataTableOutput("table_merged")
                            ),
                            conditionalPanel(
                              condition = "!output.use_dt",
                              tableOutput("table_merged_simple")
                            ))
                 ))
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values to store datasets
  values <- reactiveValues(
    dataset1 = NULL,
    dataset2 = NULL,
    merged_data = NULL
  )
  
  # Generate initial data
  observe({
    dataset1 <- generate_gene_dataset(
      n_genes = input$n_genes,
      seed = input$seed
    )
    dataset2 <- generate_correlated_dataset(
      dataset1,
      correlation = input$correlation,
      seed = input$seed
    )
    merged_data <- merge_datasets(dataset1, dataset2)
    
    values$dataset1 <- dataset1
    values$dataset2 <- dataset2
    values$merged_data <- merged_data
  })
  
  # Regenerate data when button is clicked
  observeEvent(input$generate_data, {
    dataset1 <- generate_gene_dataset(
      n_genes = input$n_genes,
      seed = input$seed
    )
    dataset2 <- generate_correlated_dataset(
      dataset1,
      correlation = input$correlation,
      seed = input$seed
    )
    merged_data <- merge_datasets(dataset1, dataset2)
    
    values$dataset1 <- dataset1
    values$dataset2 <- dataset2
    values$merged_data <- merged_data
  })
  
  # Summary output
  output$summary <- renderText({
    if (!is.null(values$merged_data)) {
      paste(
        "Genes:", nrow(values$dataset1), "\n",
        "Dataset 1 significant:", sum(values$dataset1$padj < input$fdr_threshold), "\n",
        "Dataset 2 significant:", sum(values$dataset2$padj < input$fdr_threshold), "\n",
        "Both significant:", 
        sum(values$merged_data$padj_ds1 < input$fdr_threshold & 
            values$merged_data$padj_ds2 < input$fdr_threshold)
      )
    }
  })
  
  # Plots
  output$scatter_basic <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_basic(
        values$merged_data,
        fdr_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  output$scatter_labels <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_with_labels(
        values$merged_data,
        fdr_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        max_labels = input$max_labels,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  output$ma_plot1 <- renderPlot({
    if (!is.null(values$dataset1)) {
      plot_ma(
        values$dataset1,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "MA Plot - Dataset 1"
      )
    }
  })
  
  output$ma_plot2 <- renderPlot({
    if (!is.null(values$dataset2)) {
      plot_ma(
        values$dataset2,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "MA Plot - Dataset 2"
      )
    }
  })
  
  output$volcano_plot1 <- renderPlot({
    if (!is.null(values$dataset1)) {
      plot_volcano(
        values$dataset1,
        padj_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "Volcano Plot - Dataset 1"
      )
    }
  })
  
  output$volcano_plot2 <- renderPlot({
    if (!is.null(values$dataset2)) {
      plot_volcano(
        values$dataset2,
        padj_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "Volcano Plot - Dataset 2"
      )
    }
  })
  
  output$scatter_multi <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_multi_category(
        values$merged_data,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  # Output flag for DT availability
  output$use_dt <- reactive(use_dt)
  outputOptions(output, "use_dt", suspendWhenHidden = FALSE)
  
  # Data tables
  if (use_dt) {
    output$table1 <- DT::renderDataTable({
      if (!is.null(values$dataset1)) {
        DT::datatable(values$dataset1, options = list(pageLength = 25))
      }
    })
    
    output$table2 <- DT::renderDataTable({
      if (!is.null(values$dataset2)) {
        DT::datatable(values$dataset2, options = list(pageLength = 25))
      }
    })
    
    output$table_merged <- DT::renderDataTable({
      if (!is.null(values$merged_data)) {
        DT::datatable(values$merged_data, options = list(pageLength = 25))
      }
    })
  } else {
    output$table1_simple <- renderTable({
      if (!is.null(values$dataset1)) {
        head(values$dataset1, 100)
      }
    })
    
    output$table2_simple <- renderTable({
      if (!is.null(values$dataset2)) {
        head(values$dataset2, 100)
      }
    })
    
    output$table_merged_simple <- renderTable({
      if (!is.null(values$merged_data)) {
        head(values$merged_data, 100)
      }
    })
  }
}

# Run the application
if (!interactive()) {
  shinyApp(ui = ui, server = server, options = list(port = 3838, host = "0.0.0.0"))
} else {
  shinyApp(ui = ui, server = server)
}


```


Testing the refactored script:



Checking if Shiny is installed, then updating tasks and changelog:



Relating the new files to the documentation and completing the final task:



Checking task status:



Checking off the script update task:



Creating a quick usage guide for the Shiny app:


```markdown
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

```


## Completed tasks

1. Refactored code into library modules:
   - `lib/gene_expr_data.R` — Data generation functions (`generate_gene_dataset()`, `generate_correlated_dataset()`, `merge_datasets()`)
   - `lib/gene_expr_plots.R` — Plotting functions for all visualization types

2. Created Shiny app (`app.R`):
   - Interactive interface with adjustable parameters
   - 7 visualization tabs (scatter plots, MA plots, volcano plots, multi-category)
   - Data tables for browsing datasets
   - Real-time parameter adjustment

3. Updated script (`scripts/create_test_graphs.R`):
   - Refactored to use library functions
   - Cleaner and more maintainable

4. Documentation:
   - Created `README_SHINY.md` with usage instructions
   - Updated changelog with all changes
   - Related all new files to documentation

## How to use

**Run the Shiny app:**
```bash
cd /tmp/eren-test/ttmp/2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison
Rscript app.R
```

Then open your browser to `http://localhost:3838` (or the URL shown in the console).

**Generate static plots:**
```bash
Rscript scripts/create_test_graphs.R
```

The app lets you:
- Adjust number of genes, correlation, and seed
- Change FDR thresholds, fold change cutoffs, and visualization parameters
- Explore all plot types interactively
- Browse the generated data tables

All tasks are complete and documented. The code is modular and reusable.

```shellscript
cd /tmp/eren-test/ttmp/2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison
Rscript app.R
```

```shellscript
Rscript scripts/create_test_graphs.R
```

---

**User**

run it in the background.

---

**Cursor**


Starting the Shiny app in the background:



Checking if the Shiny app is running and which port it's using:



Checking if the app started and reviewing the app.R configuration:



Fixing the DT package handling in the app:


```r
#!/usr/bin/env Rscript
# Shiny App for Interactive Gene Expression Visualization
# Run with: Rscript app.R
# Or: shiny::runApp()

library(shiny)
library(ggplot2)
library(ggrepel)

# Load DT if available, otherwise use basic tables
if (requireNamespace("DT", quietly = TRUE)) {
  library(DT)
  use_dt <- TRUE
} else {
  use_dt <- FALSE
}

# Source library functions
source("lib/gene_expr_data.R")
source("lib/gene_expr_plots.R")

# UI
ui <- fluidPage(
  titlePanel("Gene Expression Comparison Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Data Generation Parameters"),
      
      numericInput("n_genes", 
                   "Number of genes:", 
                   value = 1000, 
                   min = 100, 
                   max = 10000, 
                   step = 100),
      
      numericInput("correlation", 
                   "Dataset correlation:", 
                   value = 0.7, 
                   min = 0, 
                   max = 1, 
                   step = 0.1),
      
      numericInput("seed", 
                   "Random seed:", 
                   value = 42, 
                   min = 1),
      
      br(),
      h3("Plot Parameters"),
      
      numericInput("fdr_threshold", 
                   "FDR threshold:", 
                   value = 0.05, 
                   min = 0.001, 
                   max = 0.5, 
                   step = 0.01),
      
      numericInput("fc_threshold", 
                   "Fold change threshold:", 
                   value = 1, 
                   min = 0.1, 
                   max = 5, 
                   step = 0.1),
      
      numericInput("point_size", 
                   "Point size:", 
                   value = 1.2, 
                   min = 0.5, 
                   max = 5, 
                   step = 0.1),
      
      numericInput("point_alpha", 
                   "Point transparency:", 
                   value = 0.7, 
                   min = 0.1, 
                   max = 1, 
                   step = 0.1),
      
      br(),
      actionButton("generate_data", "Generate New Data", class = "btn-primary"),
      
      br(), br(),
      h4("Dataset Summary"),
      verbatimTextOutput("summary")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Scatter Plot (Basic)",
                 plotOutput("scatter_basic", height = "600px")),
        
        tabPanel("Scatter Plot (With Labels)",
                 numericInput("max_labels", 
                              "Max labels:", 
                              value = 20, 
                              min = 5, 
                              max = 100, 
                              step = 5,
                              width = "150px"),
                 plotOutput("scatter_labels", height = "600px")),
        
        tabPanel("MA Plot - Dataset 1",
                 plotOutput("ma_plot1", height = "600px")),
        
        tabPanel("MA Plot - Dataset 2",
                 plotOutput("ma_plot2", height = "600px")),
        
        tabPanel("Volcano Plot - Dataset 1",
                 plotOutput("volcano_plot1", height = "600px")),
        
        tabPanel("Volcano Plot - Dataset 2",
                 plotOutput("volcano_plot2", height = "600px")),
        
        tabPanel("Multi-Category Scatter",
                 plotOutput("scatter_multi", height = "600px")),
        
        tabPanel("Data Tables",
                 tabsetPanel(
                   tabPanel("Dataset 1",
                            tableOutput("table1_simple")),
                   tabPanel("Dataset 2",
                            tableOutput("table2_simple")),
                   tabPanel("Merged",
                            tableOutput("table_merged_simple"))
                 ))
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values to store datasets
  values <- reactiveValues(
    dataset1 = NULL,
    dataset2 = NULL,
    merged_data = NULL
  )
  
  # Generate initial data
  observe({
    dataset1 <- generate_gene_dataset(
      n_genes = input$n_genes,
      seed = input$seed
    )
    dataset2 <- generate_correlated_dataset(
      dataset1,
      correlation = input$correlation,
      seed = input$seed
    )
    merged_data <- merge_datasets(dataset1, dataset2)
    
    values$dataset1 <- dataset1
    values$dataset2 <- dataset2
    values$merged_data <- merged_data
  })
  
  # Regenerate data when button is clicked
  observeEvent(input$generate_data, {
    dataset1 <- generate_gene_dataset(
      n_genes = input$n_genes,
      seed = input$seed
    )
    dataset2 <- generate_correlated_dataset(
      dataset1,
      correlation = input$correlation,
      seed = input$seed
    )
    merged_data <- merge_datasets(dataset1, dataset2)
    
    values$dataset1 <- dataset1
    values$dataset2 <- dataset2
    values$merged_data <- merged_data
  })
  
  # Summary output
  output$summary <- renderText({
    if (!is.null(values$merged_data)) {
      paste(
        "Genes:", nrow(values$dataset1), "\n",
        "Dataset 1 significant:", sum(values$dataset1$padj < input$fdr_threshold), "\n",
        "Dataset 2 significant:", sum(values$dataset2$padj < input$fdr_threshold), "\n",
        "Both significant:", 
        sum(values$merged_data$padj_ds1 < input$fdr_threshold & 
            values$merged_data$padj_ds2 < input$fdr_threshold)
      )
    }
  })
  
  # Plots
  output$scatter_basic <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_basic(
        values$merged_data,
        fdr_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  output$scatter_labels <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_with_labels(
        values$merged_data,
        fdr_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        max_labels = input$max_labels,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  output$ma_plot1 <- renderPlot({
    if (!is.null(values$dataset1)) {
      plot_ma(
        values$dataset1,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "MA Plot - Dataset 1"
      )
    }
  })
  
  output$ma_plot2 <- renderPlot({
    if (!is.null(values$dataset2)) {
      plot_ma(
        values$dataset2,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "MA Plot - Dataset 2"
      )
    }
  })
  
  output$volcano_plot1 <- renderPlot({
    if (!is.null(values$dataset1)) {
      plot_volcano(
        values$dataset1,
        padj_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "Volcano Plot - Dataset 1"
      )
    }
  })
  
  output$volcano_plot2 <- renderPlot({
    if (!is.null(values$dataset2)) {
      plot_volcano(
        values$dataset2,
        padj_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "Volcano Plot - Dataset 2"
      )
    }
  })
  
  output$scatter_multi <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_multi_category(
        values$merged_data,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  # Output flag for DT availability
  output$use_dt <- reactive(use_dt)
  outputOptions(output, "use_dt", suspendWhenHidden = FALSE)
  
  # Data tables
  if (use_dt) {
    output$table1 <- DT::renderDataTable({
      if (!is.null(values$dataset1)) {
        DT::datatable(values$dataset1, options = list(pageLength = 25))
      }
    })
    
    output$table2 <- DT::renderDataTable({
      if (!is.null(values$dataset2)) {
        DT::datatable(values$dataset2, options = list(pageLength = 25))
      }
    })
    
    output$table_merged <- DT::renderDataTable({
      if (!is.null(values$merged_data)) {
        DT::datatable(values$merged_data, options = list(pageLength = 25))
      }
    })
  } else {
    output$table1_simple <- renderTable({
      if (!is.null(values$dataset1)) {
        head(values$dataset1, 100)
      }
    })
    
    output$table2_simple <- renderTable({
      if (!is.null(values$dataset2)) {
        head(values$dataset2, 100)
      }
    })
    
    output$table_merged_simple <- renderTable({
      if (!is.null(values$merged_data)) {
        head(values$merged_data, 100)
      }
    })
  }
}

# Run the application
if (!interactive()) {
  shinyApp(ui = ui, server = server, options = list(port = 3838, host = "0.0.0.0"))
} else {
  shinyApp(ui = ui, server = server)
}


```

```r
#!/usr/bin/env Rscript
# Shiny App for Interactive Gene Expression Visualization
# Run with: Rscript app.R
# Or: shiny::runApp()

library(shiny)
library(ggplot2)
library(ggrepel)

# Load DT if available, otherwise use basic tables
if (requireNamespace("DT", quietly = TRUE)) {
  library(DT)
  use_dt <- TRUE
} else {
  use_dt <- FALSE
}

# Source library functions
source("lib/gene_expr_data.R")
source("lib/gene_expr_plots.R")

# UI
ui <- fluidPage(
  titlePanel("Gene Expression Comparison Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Data Generation Parameters"),
      
      numericInput("n_genes", 
                   "Number of genes:", 
                   value = 1000, 
                   min = 100, 
                   max = 10000, 
                   step = 100),
      
      numericInput("correlation", 
                   "Dataset correlation:", 
                   value = 0.7, 
                   min = 0, 
                   max = 1, 
                   step = 0.1),
      
      numericInput("seed", 
                   "Random seed:", 
                   value = 42, 
                   min = 1),
      
      br(),
      h3("Plot Parameters"),
      
      numericInput("fdr_threshold", 
                   "FDR threshold:", 
                   value = 0.05, 
                   min = 0.001, 
                   max = 0.5, 
                   step = 0.01),
      
      numericInput("fc_threshold", 
                   "Fold change threshold:", 
                   value = 1, 
                   min = 0.1, 
                   max = 5, 
                   step = 0.1),
      
      numericInput("point_size", 
                   "Point size:", 
                   value = 1.2, 
                   min = 0.5, 
                   max = 5, 
                   step = 0.1),
      
      numericInput("point_alpha", 
                   "Point transparency:", 
                   value = 0.7, 
                   min = 0.1, 
                   max = 1, 
                   step = 0.1),
      
      br(),
      actionButton("generate_data", "Generate New Data", class = "btn-primary"),
      
      br(), br(),
      h4("Dataset Summary"),
      verbatimTextOutput("summary")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Scatter Plot (Basic)",
                 plotOutput("scatter_basic", height = "600px")),
        
        tabPanel("Scatter Plot (With Labels)",
                 numericInput("max_labels", 
                              "Max labels:", 
                              value = 20, 
                              min = 5, 
                              max = 100, 
                              step = 5,
                              width = "150px"),
                 plotOutput("scatter_labels", height = "600px")),
        
        tabPanel("MA Plot - Dataset 1",
                 plotOutput("ma_plot1", height = "600px")),
        
        tabPanel("MA Plot - Dataset 2",
                 plotOutput("ma_plot2", height = "600px")),
        
        tabPanel("Volcano Plot - Dataset 1",
                 plotOutput("volcano_plot1", height = "600px")),
        
        tabPanel("Volcano Plot - Dataset 2",
                 plotOutput("volcano_plot2", height = "600px")),
        
        tabPanel("Multi-Category Scatter",
                 plotOutput("scatter_multi", height = "600px")),
        
        tabPanel("Data Tables",
                 tabsetPanel(
                   tabPanel("Dataset 1",
                            tableOutput("table1_simple")),
                   tabPanel("Dataset 2",
                            tableOutput("table2_simple")),
                   tabPanel("Merged",
                            tableOutput("table_merged_simple"))
                 ))
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values to store datasets
  values <- reactiveValues(
    dataset1 = NULL,
    dataset2 = NULL,
    merged_data = NULL
  )
  
  # Generate initial data
  observe({
    dataset1 <- generate_gene_dataset(
      n_genes = input$n_genes,
      seed = input$seed
    )
    dataset2 <- generate_correlated_dataset(
      dataset1,
      correlation = input$correlation,
      seed = input$seed
    )
    merged_data <- merge_datasets(dataset1, dataset2)
    
    values$dataset1 <- dataset1
    values$dataset2 <- dataset2
    values$merged_data <- merged_data
  })
  
  # Regenerate data when button is clicked
  observeEvent(input$generate_data, {
    dataset1 <- generate_gene_dataset(
      n_genes = input$n_genes,
      seed = input$seed
    )
    dataset2 <- generate_correlated_dataset(
      dataset1,
      correlation = input$correlation,
      seed = input$seed
    )
    merged_data <- merge_datasets(dataset1, dataset2)
    
    values$dataset1 <- dataset1
    values$dataset2 <- dataset2
    values$merged_data <- merged_data
  })
  
  # Summary output
  output$summary <- renderText({
    if (!is.null(values$merged_data)) {
      paste(
        "Genes:", nrow(values$dataset1), "\n",
        "Dataset 1 significant:", sum(values$dataset1$padj < input$fdr_threshold), "\n",
        "Dataset 2 significant:", sum(values$dataset2$padj < input$fdr_threshold), "\n",
        "Both significant:", 
        sum(values$merged_data$padj_ds1 < input$fdr_threshold & 
            values$merged_data$padj_ds2 < input$fdr_threshold)
      )
    }
  })
  
  # Plots
  output$scatter_basic <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_basic(
        values$merged_data,
        fdr_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  output$scatter_labels <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_with_labels(
        values$merged_data,
        fdr_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        max_labels = input$max_labels,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  output$ma_plot1 <- renderPlot({
    if (!is.null(values$dataset1)) {
      plot_ma(
        values$dataset1,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "MA Plot - Dataset 1"
      )
    }
  })
  
  output$ma_plot2 <- renderPlot({
    if (!is.null(values$dataset2)) {
      plot_ma(
        values$dataset2,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "MA Plot - Dataset 2"
      )
    }
  })
  
  output$volcano_plot1 <- renderPlot({
    if (!is.null(values$dataset1)) {
      plot_volcano(
        values$dataset1,
        padj_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "Volcano Plot - Dataset 1"
      )
    }
  })
  
  output$volcano_plot2 <- renderPlot({
    if (!is.null(values$dataset2)) {
      plot_volcano(
        values$dataset2,
        padj_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "Volcano Plot - Dataset 2"
      )
    }
  })
  
  output$scatter_multi <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_multi_category(
        values$merged_data,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  # Data tables
  output$table1_simple <- renderTable({
    if (!is.null(values$dataset1)) {
      head(values$dataset1, 100)
    }
  })
  
  output$table2_simple <- renderTable({
    if (!is.null(values$dataset2)) {
      head(values$dataset2, 100)
    }
  })
  
  output$table_merged_simple <- renderTable({
    if (!is.null(values$merged_data)) {
      head(values$merged_data, 100)
    }
  })
}

# Run the application
if (!interactive()) {
  shinyApp(ui = ui, server = server, options = list(port = 3838, host = "0.0.0.0"))
} else {
  shinyApp(ui = ui, server = server)
}


```

```r
#!/usr/bin/env Rscript
# Shiny App for Interactive Gene Expression Visualization
# Run with: Rscript app.R
# Or: shiny::runApp()

library(shiny)
library(ggplot2)
library(ggrepel)

# Note: Using basic table output (DT package optional)

# Source library functions
source("lib/gene_expr_data.R")
source("lib/gene_expr_plots.R")

# UI
ui <- fluidPage(
  titlePanel("Gene Expression Comparison Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Data Generation Parameters"),
      
      numericInput("n_genes", 
                   "Number of genes:", 
                   value = 1000, 
                   min = 100, 
                   max = 10000, 
                   step = 100),
      
      numericInput("correlation", 
                   "Dataset correlation:", 
                   value = 0.7, 
                   min = 0, 
                   max = 1, 
                   step = 0.1),
      
      numericInput("seed", 
                   "Random seed:", 
                   value = 42, 
                   min = 1),
      
      br(),
      h3("Plot Parameters"),
      
      numericInput("fdr_threshold", 
                   "FDR threshold:", 
                   value = 0.05, 
                   min = 0.001, 
                   max = 0.5, 
                   step = 0.01),
      
      numericInput("fc_threshold", 
                   "Fold change threshold:", 
                   value = 1, 
                   min = 0.1, 
                   max = 5, 
                   step = 0.1),
      
      numericInput("point_size", 
                   "Point size:", 
                   value = 1.2, 
                   min = 0.5, 
                   max = 5, 
                   step = 0.1),
      
      numericInput("point_alpha", 
                   "Point transparency:", 
                   value = 0.7, 
                   min = 0.1, 
                   max = 1, 
                   step = 0.1),
      
      br(),
      actionButton("generate_data", "Generate New Data", class = "btn-primary"),
      
      br(), br(),
      h4("Dataset Summary"),
      verbatimTextOutput("summary")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Scatter Plot (Basic)",
                 plotOutput("scatter_basic", height = "600px")),
        
        tabPanel("Scatter Plot (With Labels)",
                 numericInput("max_labels", 
                              "Max labels:", 
                              value = 20, 
                              min = 5, 
                              max = 100, 
                              step = 5,
                              width = "150px"),
                 plotOutput("scatter_labels", height = "600px")),
        
        tabPanel("MA Plot - Dataset 1",
                 plotOutput("ma_plot1", height = "600px")),
        
        tabPanel("MA Plot - Dataset 2",
                 plotOutput("ma_plot2", height = "600px")),
        
        tabPanel("Volcano Plot - Dataset 1",
                 plotOutput("volcano_plot1", height = "600px")),
        
        tabPanel("Volcano Plot - Dataset 2",
                 plotOutput("volcano_plot2", height = "600px")),
        
        tabPanel("Multi-Category Scatter",
                 plotOutput("scatter_multi", height = "600px")),
        
        tabPanel("Data Tables",
                 tabsetPanel(
                   tabPanel("Dataset 1",
                            tableOutput("table1_simple")),
                   tabPanel("Dataset 2",
                            tableOutput("table2_simple")),
                   tabPanel("Merged",
                            tableOutput("table_merged_simple"))
                 ))
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values to store datasets
  values <- reactiveValues(
    dataset1 = NULL,
    dataset2 = NULL,
    merged_data = NULL
  )
  
  # Generate initial data
  observe({
    dataset1 <- generate_gene_dataset(
      n_genes = input$n_genes,
      seed = input$seed
    )
    dataset2 <- generate_correlated_dataset(
      dataset1,
      correlation = input$correlation,
      seed = input$seed
    )
    merged_data <- merge_datasets(dataset1, dataset2)
    
    values$dataset1 <- dataset1
    values$dataset2 <- dataset2
    values$merged_data <- merged_data
  })
  
  # Regenerate data when button is clicked
  observeEvent(input$generate_data, {
    dataset1 <- generate_gene_dataset(
      n_genes = input$n_genes,
      seed = input$seed
    )
    dataset2 <- generate_correlated_dataset(
      dataset1,
      correlation = input$correlation,
      seed = input$seed
    )
    merged_data <- merge_datasets(dataset1, dataset2)
    
    values$dataset1 <- dataset1
    values$dataset2 <- dataset2
    values$merged_data <- merged_data
  })
  
  # Summary output
  output$summary <- renderText({
    if (!is.null(values$merged_data)) {
      paste(
        "Genes:", nrow(values$dataset1), "\n",
        "Dataset 1 significant:", sum(values$dataset1$padj < input$fdr_threshold), "\n",
        "Dataset 2 significant:", sum(values$dataset2$padj < input$fdr_threshold), "\n",
        "Both significant:", 
        sum(values$merged_data$padj_ds1 < input$fdr_threshold & 
            values$merged_data$padj_ds2 < input$fdr_threshold)
      )
    }
  })
  
  # Plots
  output$scatter_basic <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_basic(
        values$merged_data,
        fdr_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  output$scatter_labels <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_with_labels(
        values$merged_data,
        fdr_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        max_labels = input$max_labels,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  output$ma_plot1 <- renderPlot({
    if (!is.null(values$dataset1)) {
      plot_ma(
        values$dataset1,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "MA Plot - Dataset 1"
      )
    }
  })
  
  output$ma_plot2 <- renderPlot({
    if (!is.null(values$dataset2)) {
      plot_ma(
        values$dataset2,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "MA Plot - Dataset 2"
      )
    }
  })
  
  output$volcano_plot1 <- renderPlot({
    if (!is.null(values$dataset1)) {
      plot_volcano(
        values$dataset1,
        padj_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "Volcano Plot - Dataset 1"
      )
    }
  })
  
  output$volcano_plot2 <- renderPlot({
    if (!is.null(values$dataset2)) {
      plot_volcano(
        values$dataset2,
        padj_threshold = input$fdr_threshold,
        fc_threshold = input$fc_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha,
        title = "Volcano Plot - Dataset 2"
      )
    }
  })
  
  output$scatter_multi <- renderPlot({
    if (!is.null(values$merged_data)) {
      plot_scatter_multi_category(
        values$merged_data,
        padj_threshold = input$fdr_threshold,
        point_size = input$point_size,
        point_alpha = input$point_alpha
      )
    }
  })
  
  # Data tables
  output$table1_simple <- renderTable({
    if (!is.null(values$dataset1)) {
      head(values$dataset1, 100)
    }
  })
  
  output$table2_simple <- renderTable({
    if (!is.null(values$dataset2)) {
      head(values$dataset2, 100)
    }
  })
  
  output$table_merged_simple <- renderTable({
    if (!is.null(values$merged_data)) {
      head(values$merged_data, 100)
    }
  })
}

# Run the application
if (!interactive()) {
  shinyApp(ui = ui, server = server, options = list(port = 3838, host = "0.0.0.0"))
} else {
  shinyApp(ui = ui, server = server)
}


```


Restarting the app in the background:



Checking for remaining DT references:



Testing the app directly to see the error:



Starting the app in the background:



Shiny app is running in the background.



**Shiny app is running in the background.**

- **Status**: Running (PID 2932041)
- **URL**: `http://localhost:3838` (or `http://0.0.0.0:3838`)
- **Log file**: `shiny_app.log` in the ticket directory

The app is listening on port 3838. Open `http://localhost:3838` in your browser to use it.

To stop the app later:
```bash
pkill -f "Rscript app.R"
```

Or find the process ID and kill it:
```bash
ps aux | grep "Rscript app.R"
kill <PID>
```

```shellscript
pkill -f "Rscript app.R"
```

```shellscript
ps aux | grep "Rscript app.R"
kill <PID>
```

