---
Title: Onboarding Guide for Gene Expression Visualization Project
Ticket: COMPARISON-GRAPH
Status: active
Topics:
    - analysis
    - visualization
DocType: tutorial
Intent: long-term
Owners: []
RelatedFiles:
    - Path: 2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/tutorial/app.R
      Note: Interactive Shiny application for exploring visualizations
    - Path: 2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/tutorial/lib/gene_expr_data.R
      Note: Data generation functions for creating test datasets
    - Path: 2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/tutorial/lib/gene_expr_plots.R
      Note: Plotting function library with reusable visualization code
    - Path: 2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/tutorial/reference/01-ggplot2-for-comparing-gene-expression-datasets.md
      Note: Primary reference guide with complete code examples
    - Path: 2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/tutorial/scripts/create_test_graphs.R
      Note: Main test script demonstrating all visualization types
ExternalSources: []
Summary: Comprehensive onboarding guide for interns joining the gene expression visualization project, covering project context, key concepts, available resources, and next steps for contributing effectively
LastUpdated: 2025-12-04T08:14:44.039035551-05:00
---



# Onboarding Guide for Gene Expression Visualization Project

## Overview

This guide helps new team members quickly understand the COMPARISON-GRAPH project and start contributing effectively. You'll learn what this project does, why it exists, what's already been built, and how to navigate the codebase and documentation. By the end, you'll know where to find key resources, understand the technical approach, and have a clear path forward for your contributions.

**Learning objectives:**
- Understand the project's purpose and scope
- Navigate the ticket workspace and locate key documents
- Grasp core concepts in gene expression visualization
- Identify available tools, scripts, and reference materials
- Know how to continue research and implementation work

## Prerequisites

Before diving in, ensure you have:

- **Basic R knowledge**: Familiarity with R syntax, data frames, and basic plotting
- **ggplot2 basics**: Understanding of ggplot2's grammar of graphics (aes, geom_*, labs, theme)
- **R environment**: R installed with ggplot2, ggrepel, and optionally ggpubr packages
- **Git familiarity**: Ability to navigate repositories and understand file structures
- **Documentation tools**: Access to docmgr (for updating documentation) and glaze (for help)

**Optional but helpful:**
- Experience with gene expression analysis (DESeq2, edgeR, or similar)
- Understanding of statistical concepts (p-values, FDR, log fold changes)
- Familiarity with Shiny applications (for the interactive app)

## Project Context

### Purpose and Scope

The COMPARISON-GRAPH project creates a comprehensive reference guide and toolkit for visualizing comparisons between two gene expression datasets. Researchers often need to compare results from different experiments, conditions, or analysis pipelines, and this project provides copy-paste-ready ggplot2 code for common visualization patterns.

**Core problem solved:** When comparing two gene expression datasets, researchers need to visualize:
- Correlation between log fold changes from both datasets
- Which genes are statistically significant (FDR < 0.05)
- Magnitude of changes and their statistical significance
- Agreement or disagreement between datasets

**Solution provided:** Ready-to-use R code for three visualization types:
- **Scatter plots**: Direct comparison of log fold changes between datasets
- **MA plots**: Log fold change (M) vs average expression (A) for individual datasets
- **Volcano plots**: Statistical significance (-log10 p-value) vs fold change magnitude

### What's Been Accomplished

The project has established a solid foundation with both documentation and working code:

**Documentation:**
- Complete reference guide with code examples for all three plot types
- Data structure requirements and merging workflows
- Troubleshooting section for common issues
- Usage examples for different scenarios

**Implementation:**
- Test data generation functions (`lib/gene_expr_data.R`)
- Plotting function library (`lib/gene_expr_plots.R`)
- Complete test script (`scripts/create_test_graphs.R`)
- Interactive Shiny application (`app.R`)
- Sample output plots in `various/plots/`

**Current status:** The project is **active** and ready for enhancements, additional plot types, or improvements to existing visualizations.

## Key Concepts

Understanding these concepts will help you work effectively with the codebase and contribute meaningfully.

### Gene Expression Data Structure

Gene expression analysis produces datasets with specific columns that visualization code expects:

**Required columns for individual datasets:**
- `gene_id`: Unique identifier (character or factor)
- `log2FoldChange`: Log2 fold change value (numeric)
- `padj`: Adjusted p-value (FDR) after multiple testing correction (0-1)
- `baseMean`: Mean expression level across conditions (positive numeric)

**For comparison datasets (merged):**
- `logFC_dataset1` / `logFC_dataset2`: Log fold changes from each dataset
- `padj_ds1`, `padj_ds2`: Adjusted p-values from each dataset
- `FDR`: Combined FDR (typically `pmax(padj_ds1, padj_ds2)`)

**Key insight:** `padj` and `FDR` represent the same concept (False Discovery Rate after Benjamini-Hochberg correction). The codebase uses both terms depending on context.

### Visualization Types and Their Use Cases

Each plot type answers different research questions:

**Scatter plots** answer: "Do the two datasets agree? Which genes show consistent changes?"
- X-axis: log fold change from dataset 1
- Y-axis: log fold change from dataset 2
- Reference line (slope=1): Perfect agreement between datasets
- Color coding: Significant (FDR < 0.05) vs non-significant genes

**MA plots** answer: "Are fold changes consistent across expression levels?"
- X-axis: Log2 mean expression (A = average)
- Y-axis: Log2 fold change (M = minus, the difference)
- Reveals: Whether fold changes depend on expression level
- Common pattern: Higher variance at low expression levels

**Volcano plots** answer: "Which genes are both statistically significant and biologically meaningful?"
- X-axis: Log2 fold change (biological significance)
- Y-axis: -Log10 adjusted p-value (statistical significance)
- Thresholds: Vertical lines at ±1 (2-fold change), horizontal at -log10(0.05)
- Quadrants: Up-regulated significant, down-regulated significant, non-significant

### Significance Thresholds

The codebase uses standard thresholds that you should understand:

- **FDR < 0.05**: Standard threshold for statistical significance (5% false discovery rate)
- **|log2FC| > 1**: Represents a 2-fold change (2^1 = 2)
- **|log2FC| > 2**: Represents a 4-fold change (2^2 = 4), used for labeling extreme cases

**Design decision:** The reference uses `abs(logFC) > 2` for gene labeling to avoid overcrowding plots. Adjust thresholds based on your specific research needs.

## Navigating the Workspace

The ticket workspace follows a standard structure that organizes different types of content:

```
COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison/
├── index.md                    # Ticket overview and entry point
├── tasks.md                    # Current task list
├── changelog.md                # History of changes
├── reference/                  # Reference documentation
│   └── 01-ggplot2-for-comparing-gene-expression-datasets.md
├── lib/                        # Reusable R functions
│   ├── gene_expr_data.R       # Data generation functions
│   └── gene_expr_plots.R      # Plotting functions
├── scripts/                    # Implementation scripts
│   └── create_test_graphs.R   # Main test script
├── app.R                       # Shiny interactive application
├── various/                    # Working notes and outputs
│   └── plots/                 # Generated plot images
└── log/                        # Implementation diary
    └── 01-2025-12-04-test-datasets-and-graphs.md
```

**Where to start:**
1. **index.md**: Read first for project overview
2. **reference/01-ggplot2-for-comparing-gene-expression-datasets.md**: Complete code reference
3. **scripts/create_test_graphs.R**: Working implementation example
4. **lib/**: Reusable functions you can extend

## Available Resources

### Documentation

**Primary reference:**
- `reference/01-ggplot2-for-comparing-gene-expression-datasets.md`: Complete guide with all code examples, data preparation steps, troubleshooting, and usage patterns

**Implementation notes:**
- `log/01-2025-12-04-test-datasets-and-graphs.md`: Development diary showing what was tried and learned

**External resources:**
- ggpubr documentation: https://rpkgs.datanovia.com/ggpubr/reference/ggmaplot.html
- ggplot2 documentation: https://ggplot2.tidyverse.org/

### Code Libraries

**Data generation (`lib/gene_expr_data.R`):**
- `generate_gene_dataset()`: Creates synthetic gene expression data
- `generate_correlated_dataset()`: Creates correlated dataset for comparison
- `merge_datasets()`: Combines two datasets for comparison plots
- `print_dataset_summary()`: Diagnostic summary output

**Plotting functions (`lib/gene_expr_plots.R`):**
- `plot_scatter_basic()`: Basic scatter plot comparison
- `plot_scatter_with_labels()`: Scatter plot with gene labels
- `plot_ma()`: MA plot for individual datasets
- `plot_volcano()`: Volcano plot visualization

**Usage pattern:**
```r
# Source the libraries
source("lib/gene_expr_data.R")
source("lib/gene_expr_plots.R")

# Generate test data
dataset1 <- generate_gene_dataset(n_genes = 1000, seed = 42)
dataset2 <- generate_correlated_dataset(dataset1, correlation = 0.7, seed = 42)
merged_data <- merge_datasets(dataset1, dataset2)

# Create plots
p1 <- plot_scatter_basic(merged_data)
p2 <- plot_ma(dataset1, title = "MA Plot - Dataset 1")
```

### Interactive Tools

**Shiny application (`app.R`):**
- Interactive exploration of gene expression visualizations
- Upload your own data or use generated test data
- Adjust parameters and see results immediately
- Useful for: Exploring data, demonstrating to researchers, rapid prototyping

**To run:**
```r
shiny::runApp("app.R")
```

### Test Data and Outputs

**Sample outputs (`various/plots/`):**
- Pre-generated plots showing expected output
- CSV files with test datasets (`dataset1.csv`, `dataset2.csv`, `merged_data.csv`)
- Use these to verify your environment works correctly

**To regenerate:**
```r
# Run the test script
source("scripts/create_test_graphs.R")
```

## Next Steps for Contributing

Now that you understand the project, here's how to contribute effectively:

### Immediate Actions

1. **Set up your environment:**
   ```r
   # Install required packages
   install.packages(c("ggplot2", "ggrepel", "ggpubr"))
   
   # Verify installation
   library(ggplot2)
   library(ggrepel)
   ```

2. **Run the test script:**
   ```r
   # Navigate to the ticket directory
   setwd("ttmp/2025/12/04/COMPARISON-GRAPH--ggplot2-visualization-for-gene-expression-comparison")
   
   # Run the test script
   source("scripts/create_test_graphs.R")
   
   # Verify plots were created
   list.files("various/plots/")
   ```

3. **Explore the Shiny app:**
   ```r
   shiny::runApp("app.R")
   ```

### Potential Enhancement Areas

Based on the current state, consider these directions:

**New visualization types:**
- Heatmaps comparing multiple datasets
- Correlation matrices between datasets
- Venn diagrams for significant gene overlap
- Density plots for fold change distributions

**Code improvements:**
- Add more customization options to plotting functions
- Improve error handling and input validation
- Add unit tests for data generation functions
- Create vignettes or additional tutorials

**Documentation enhancements:**
- Add more real-world usage examples
- Create video tutorials or interactive demos
- Expand troubleshooting section with more edge cases
- Add performance tips for large datasets (>100k genes)

**Tooling:**
- Command-line interface for batch processing
- Integration with common analysis pipelines (DESeq2, edgeR)
- Export functions for publication-quality figures

### Working with Documentation

This project uses docmgr for documentation management. When making changes:

**Updating documentation:**
```bash
# Update a document's metadata
docmgr meta update --ticket COMPARISON-GRAPH --doc-type reference \
  --field Summary --value "Updated description"

# Relate files to documentation
docmgr doc relate --ticket COMPARISON-GRAPH --doc-type reference \
  --file-note "scripts/create_test_graphs.R:Main test implementation"

# Add changelog entry
docmgr changelog update --ticket COMPARISON-GRAPH \
  --entry "Added new visualization type" \
  --file-note "lib/gene_expr_plots.R:New plotting function"
```

**Best practices:**
- Always relate files when creating or modifying docs
- Update changelog after significant changes
- Keep index.md concise (~50 lines)
- Use docmgr commands for metadata, write markdown content in your editor

### Research Workflow

When researching new topics or approaches:

1. **Create working notes:**
   - Store research in `log/` directory with date prefixes
   - Document what you tried, what worked, what didn't
   - Include code snippets and references

2. **Update reference docs:**
   - Add new findings to the reference document
   - Include code examples that are minimal and focused
   - Link to external sources in frontmatter

3. **Relate implementation files:**
   - Use `docmgr doc relate` to connect code to documentation
   - Always include notes explaining why files matter
   - Keep bidirectional links between code and docs

4. **Validate your work:**
   ```bash
   # Check for issues
   docmgr doctor --ticket COMPARISON-GRAPH
   
   # Fix any warnings
   docmgr validate frontmatter --doc <file> --suggest-fixes
   ```

## Verification

Verify you're ready to contribute by completing these checks:

**Environment check:**
- [ ] R is installed and accessible
- [ ] Required packages (ggplot2, ggrepel) install successfully
- [ ] Can run `source("scripts/create_test_graphs.R")` without errors
- [ ] Plots appear in `various/plots/` directory

**Understanding check:**
- [ ] Can explain what each plot type shows (scatter, MA, volcano)
- [ ] Understand the data structure requirements
- [ ] Know where to find code examples and reference docs
- [ ] Can navigate the workspace structure

**Tooling check:**
- [ ] docmgr is installed and working (`docmgr status`)
- [ ] Can search for docs (`docmgr doc search --query "scatter"`)
- [ ] Understand how to relate files to documentation
- [ ] Know how to update changelog and tasks

**Next steps:**
- [ ] Identified an area you want to work on
- [ ] Created a plan for your contribution
- [ ] Know who to ask questions (check Owners in index.md)

## Troubleshooting

### Common Setup Issues

**Problem: R packages won't install**
- **Solution:** Ensure you have write permissions and internet access. Try installing from CRAN mirror: `install.packages("ggplot2", repos = "https://cran.rstudio.com/")`

**Problem: Scripts can't find library files**
- **Solution:** Ensure you're in the correct directory. Use `setwd()` or run scripts with full paths. Check that `lib/` directory exists relative to your working directory.

**Problem: Plots aren't generating**
- **Solution:** Check that output directory exists (`various/plots/`). Verify data generation succeeded before plotting. Check R console for error messages.

### Documentation Issues

**Problem: docmgr commands fail**
- **Solution:** Verify you're in the workspace root (where `.ttmp.yaml` exists). Check that ticket exists: `docmgr ticket list --ticket COMPARISON-GRAPH`

**Problem: Can't relate files**
- **Solution:** Use absolute paths or paths relative to workspace root (not relative to `ttmp/`). Always include notes: `--file-note "path:reason"`

**Problem: Validation warnings**
- **Solution:** Run `docmgr doctor --ticket COMPARISON-GRAPH` to see all issues. Fix frontmatter with `docmgr meta update`. Check vocabulary: `docmgr vocab list`

### Code Issues

**Problem: Plot looks wrong**
- **Solution:** Verify data structure matches requirements (check column names). Ensure significance thresholds match your data (FDR vs padj). Check that data is numeric, not character.

**Problem: Labels overlapping**
- **Solution:** Adjust `max.overlaps` parameter in `geom_text_repel()`. Increase threshold for labeling (e.g., `abs(logFC) > 2` instead of `> 1`). Reduce label size.

**Problem: Colors not showing**
- **Solution:** Verify logical condition evaluates correctly (`FDR < 0.05` should return TRUE/FALSE). Check `scale_color_manual()` values match condition results. Ensure data doesn't have NA values breaking the condition.

## Related Resources

### Internal Documentation

- **Reference guide**: `reference/01-ggplot2-for-comparing-gene-expression-datasets.md` - Complete code reference
- **Ticket index**: `index.md` - Project overview and entry point
- **Implementation log**: `log/01-2025-12-04-test-datasets-and-graphs.md` - Development history

### External Resources

- **ggplot2 documentation**: https://ggplot2.tidyverse.org/ - Core plotting library
- **ggpubr package**: https://rpkgs.datanovia.com/ggpubr/ - Simplified plotting functions
- **ggrepel documentation**: https://ggrepel.slowkow.com/ - Non-overlapping labels
- **DESeq2 vignette**: https://bioconductor.org/packages/DESeq2/ - Common source of gene expression data

### Getting Help

- **Documentation guidelines**: `glaze help how-to-write-good-documentation-pages` - Writing style guide
- **docmgr usage**: `docmgr help how-to-use` - Complete workflow guide
- **Ticket owners**: Check `index.md` frontmatter for Owners field
- **Vocabulary**: `docmgr vocab list` - See available topics and doc types

### Quick Command Reference

```bash
# Search documentation
docmgr doc search --query "scatter plot"

# List all docs in ticket
docmgr doc list --ticket COMPARISON-GRAPH

# Update metadata
docmgr meta update --ticket COMPARISON-GRAPH --field Status --value active

# Relate files
docmgr doc relate --ticket COMPARISON-GRAPH --doc-type reference \
  --file-note "path/to/file.R:Description"

# Validate
docmgr doctor --ticket COMPARISON-GRAPH
```

---

**Remember:** This is a living document. As you learn and contribute, consider updating this guide to help the next person who joins the project. Use docmgr to keep documentation current and well-connected to the codebase.
