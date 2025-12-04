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

