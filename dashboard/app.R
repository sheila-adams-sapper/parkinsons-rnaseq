# dashboard/app.R

# Load libraries and resolve conflicts
library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(dplyr, warn.conflicts = FALSE)
library(rlang)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

# Explicitly assign the functions we want to use
select <- dplyr::select
filter <- dplyr::filter
arrange <- dplyr::arrange
mutate <- dplyr::mutate
# Explicitly resolve the filter conflict
filter <- dplyr::filter

# Load dashboard data with proper error handling
tryCatch({
    dashboard_data <- readRDS("data/dashboard_ready.rds")
    cat("Dashboard data loaded successfully\n")
    
    # Verify the data structure matches what we expect
    if (is.null(dashboard_data$results)) {
        stop("Results data is missing")
    }
    
    if (is.null(dashboard_data$metadata)) {
        stop("Metadata is missing")
    }
    
    cat("Data loaded with", nrow(dashboard_data$results), "results and", nrow(dashboard_data$metadata), "samples\n")
    
}, error = function(e) {
    cat("Error loading dashboard data:", e$message, "\n")
    stop("Cannot load dashboard data. Please regenerate with: snakemake --forcerun prepare_dashboard_data")
})

# Extract data for easier access
results_data <- dashboard_data$results
metadata_data <- dashboard_data$metadata
summary_data <- dashboard_data$summary

# Ensure required columns exist
if (!"contrast" %in% colnames(results_data)) {
    stop("Results data missing 'contrast' column")
}

if (!"padj" %in% colnames(results_data)) {
    stop("Results data missing 'padj' column")
}

if (!"log2FoldChange" %in% colnames(results_data)) {
    stop("Results data missing 'log2FoldChange' column")
}

# Get available contrasts
available_contrasts <- unique(results_data$contrast)
cat("Available contrasts:", paste(available_contrasts, collapse = ", "), "\n")

# Define the UI
ui <- dashboardPage(
    dashboardHeader(title = "Parkinson's Disease RNA-seq Analysis"),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Overview", tabName = "overview", icon = icon("home")),
            menuItem("Volcano Plot", tabName = "volcano", icon = icon("chart-line")),
            menuItem("PCA", tabName = "pca", icon = icon("project-diagram")),
            menuItem("Heatmap", tabName = "heatmap", icon = icon("th")),
            menuItem("DEG Table", tabName = "deg_table", icon = icon("table")),
            menuItem("Pathway Analysis", tabName = "pathway", icon = icon("sitemap")),
            menuItem("Gene Expression", tabName = "gene_expr", icon = icon("dna"))
        ),
        
        # Filters in sidebar
        hr(),
        h4("Filters"),
        
        selectInput("contrast_filter", 
                   "Contrast:",
                   choices = available_contrasts,
                   selected = available_contrasts[1]),
        
        numericInput("padj_threshold", 
                    "Adjusted P-value:",
                    value = 0.05, 
                    min = 0.001, 
                    max = 1, 
                    step = 0.01),
        
        numericInput("lfc_threshold", 
                    "Log2 Fold Change:",
                    value = 1, 
                    min = 0, 
                    max = 5, 
                    step = 0.1)
    ),
    
    dashboardBody(
        tabItems(
            # Overview tab
            tabItem(tabName = "overview",
                fluidRow(
                    valueBoxOutput("total_genes"),
                    valueBoxOutput("total_samples"),
                    valueBoxOutput("analysis_date")
                ),
                fluidRow(
                    box(title = "Study Overview", status = "primary", solidHeader = TRUE,
                        width = 12,
                        p(strong("This dashboard presents results from transcriptional profiling of skeletal muscle 
                    in Parkinson's disease patients compared to age-matched controls. A small subset of Parkinson's
                    disease patients completed a 16 week high intensity training program. These participants
                    submitted muscle tissue samples for pre and post training time points.")),
                        p(strong("Analysis includes", length(available_contrasts), "contrasts:",
                               paste(available_contrasts, collapse = ", "), "which can be explored using the
                               dropdown filter on the left.")),
                        p(strong("Total of", nrow(results_data), "gene results across all contrasts.")),
                        p(strong("Publication:"), "doi: 10.3389/fphys.2020.00653"),
                        p(strong("Samples:"), "12 PD patients, 12 older controls, 12 younger controls, with 5 PD patients 
                        completing pre and post training timepoints."),
                        p(strong("Analysis:"), "Differential expression analysis using DESeq2")
                    )
                )
            ),
            
            # Volcano plot tab
            tabItem(tabName = "volcano",
                fluidRow(
                    box(title = "Volcano Plot", status = "primary", solidHeader = TRUE,
                        width = 12, height = 600,
                        plotlyOutput("volcano_plot", height = "500px")
                    )
                )
            ),
            
            # PCA tab
            tabItem(tabName = "pca",
                fluidRow(
                    box(title = "Principal Component Analysis", status = "primary", solidHeader = TRUE,
                        width = 12, height = 600,
                        p("PCA analysis of the most variable genes across samples."),
                        plotlyOutput("pca_plot", height = "500px")
                    )
                )
            ),
            
            # Heatmap tab
            tabItem(tabName = "heatmap",
                fluidRow(
                    box(title = "Expression Heatmap", status = "primary", solidHeader = TRUE,
                        width = 12, height = 600,
                        p("Heatmap of top differentially expressed genes (scaled expression values)."),
                        plotlyOutput("heatmap_plot", height = "500px")
                    )
                )
            ),
            
            # DEG Table tab
            tabItem(tabName = "deg_table",
                fluidRow(
                    box(title = "Differentially Expressed Genes", status = "primary", solidHeader = TRUE,
                        width = 12,
                        DT::dataTableOutput("deg_table")
                    )
                )
            ),
            
            # Pathway Analysis tab
                        # Pathway Analysis tab - ENHANCED VERSION
            tabItem(tabName = "pathway",
                fluidRow(
                    box(title = "KEGG Pathway Enrichment Analysis", status = "primary", solidHeader = TRUE,
                        width = 12,
                        tags$div(
                            p(strong("Directional Analysis:"), 
                              "🔴 Red background = Up-regulated pathways (enhanced in treatment)", br(),
                              "🔵 Blue background = Down-regulated pathways (reduced in treatment)", br(),
                              "Yellow/Gray = Significance levels (p.adjust < 0.01, < 0.05, > 0.05)")
                        ),
                        hr(),
                        DT::dataTableOutput("pathway_table")
                    )
                )
            ),
            
            # Gene Expression tab
            tabItem(tabName = "gene_expr",
                fluidRow(
                    box(title = "Gene Expression Search", status = "primary", solidHeader = TRUE,
                        width = 12,
                        textInput("gene_search", "Enter Gene ID:", value = ""),
                        br(),
                        textOutput("gene_expr_result")
                    )
                )
            )
        ) # Close tabItems
    ) # Close dashboardBody
) # Close dashboardPage

# Define server logic
server <- function(input, output, session) {
    
    # Reactive data filtering
    filtered_data <- reactive({
        req(input$contrast_filter)
        
        data <- results_data %>%
            filter(contrast == input$contrast_filter) %>%
            filter(!is.na(padj), !is.na(log2FoldChange))
        
        return(data)
    })
    
    # Value boxes
    output$total_genes <- renderValueBox({
        valueBox(
            value = length(unique(results_data$gene)),
            subtitle = "Total Genes",
            icon = icon("dna"),
            color = "blue"
        )
    })
    
    output$total_samples <- renderValueBox({
        valueBox(
            value = nrow(metadata_data),
            subtitle = "Total Samples", 
            icon = icon("vials"),
            color = "green"
        )
    })
    
    output$analysis_date <- renderValueBox({
        valueBox(
            value = as.character(if(is.null(summary_data$analysis_date)) Sys.Date() else summary_data$analysis_date),
            subtitle = "Analysis Date",
            icon = icon("calendar"),
            color = "yellow"
        )
    })
    
    # Optimized volcano plot
    output$volcano_plot <- renderPlotly({
        req(filtered_data())
        
        data <- filtered_data()
        
        if (nrow(data) == 0) {
            p <- ggplot() + 
                annotate("text", x = 0, y = 0, label = "No data available for selected contrast") +
                theme_minimal()
            return(ggplotly(p))
        }
        
        # OPTIMIZATION: Sample data if too large
        if (nrow(data) > 10000) {
            # Keep all significant genes, sample the rest
            sig_genes <- data[data$padj < 0.05 & abs(data$log2FoldChange) > 0.5, ]
            non_sig_sample <- data[!(data$padj < 0.05 & abs(data$log2FoldChange) > 0.5), ]
            
            if (nrow(non_sig_sample) > 5000) {
                non_sig_sample <- non_sig_sample[sample(nrow(non_sig_sample), 5000), ]
            }
            
            data <- rbind(sig_genes, non_sig_sample)
        }
        
        # Create significance categories
        data$significance <- "Not Significant"
        data$significance[data$padj < input$padj_threshold & data$log2FoldChange > input$lfc_threshold] <- "Up-regulated"
        data$significance[data$padj < input$padj_threshold & data$log2FoldChange < -input$lfc_threshold] <- "Down-regulated"
        
        p <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj), color = significance, text = gene)) +
            geom_point(alpha = 0.6, size = 1) +
            scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not Significant" = "grey")) +
            geom_hline(yintercept = -log10(input$padj_threshold), linetype = "dashed") +
            geom_vline(xintercept = c(-input$lfc_threshold, input$lfc_threshold), linetype = "dashed") +
            labs(x = "Log2 Fold Change", y = "-log10(Adjusted P-value)",
                 title = paste("Volcano Plot -", input$contrast_filter)) +
            theme_minimal()
        
        ggplotly(p, tooltip = c("text", "x", "y"))
    })
    
        # Enhanced DEG table with gene symbol conversion - CORRECTED VERSION
    output$deg_table <- DT::renderDataTable({
        req(filtered_data())
        
        tryCatch({
            # Get filtered data
            data <- filtered_data()
            
            if (nrow(data) == 0) {
                return(DT::datatable(data.frame(Message = "No data available for selected contrast")))
            }
            
            # Filter for significant genes - use explicit dplyr:: prefix
            sig_data <- data %>%
                dplyr::filter(padj < input$padj_threshold, abs(log2FoldChange) > input$lfc_threshold) %>%
                dplyr::arrange(padj)
            
            if (nrow(sig_data) == 0) {
                return(DT::datatable(data.frame(Message = paste("No significant genes found with current thresholds:",
                                                              "padj <", input$padj_threshold, 
                                                              "and |log2FC| >", input$lfc_threshold))))
            }
            
            # Add gene symbol conversion function (same as heatmap)
            convert_gene_ids_to_symbols <- function(gene_ids) {
                tryCatch({
                    conversion <- clusterProfiler::bitr(gene_ids, 
                                                  fromType = "ENSEMBL", 
                                                  toType = "SYMBOL", 
                                                  OrgDb = org.Hs.eg.db)
                    
                    # Create a named vector for easy lookup
                    symbol_map <- setNames(conversion$SYMBOL, conversion$ENSEMBL)
                    
                    # Replace IDs with symbols, keep original if no symbol found
                    result <- ifelse(gene_ids %in% names(symbol_map), 
                                    symbol_map[gene_ids], 
                                    gene_ids)
                    
                    return(result)
                }, error = function(e) {
                    cat("Gene symbol conversion failed for DEG table, keeping original IDs\n")
                    return(gene_ids)
                })
            }
            
            # Convert gene IDs to symbols
            gene_symbols <- convert_gene_ids_to_symbols(sig_data$gene)
            
            # Create display table with symbols - keep original column names for formatting
            display_data <- data.frame(
                Gene_Symbol = gene_symbols,
                Gene_ID = sig_data$gene,
                log2FoldChange = round(sig_data$log2FoldChange, 3),
                padj = signif(sig_data$padj, 3),
                baseMean = round(sig_data$baseMean, 1),
                stringsAsFactors = FALSE
            )
            
            # Create the data table - apply formatting BEFORE changing column names
            dt <- DT::datatable(
                display_data, 
                options = list(
                    pageLength = 25, 
                    scrollX = TRUE,
                    columnDefs = list(
                        list(width = '150px', targets = 0),  # Gene Symbol column
                        list(width = '150px', targets = 1)   # Gene ID column
                    )
                ),
                rownames = FALSE
            ) %>%
            # Format using the original data frame column names
            DT::formatRound(c("log2FoldChange", "baseMean"), 3) %>%
            DT::formatSignif("padj", 3)
            
            # Now we can set the display column names without affecting formatting
            dt$x$data <- display_data
            dt$x$options$columnDefs <- append(dt$x$options$columnDefs, list(
                list(title = "Gene Symbol", targets = 0),
                list(title = "Gene ID", targets = 1),
                list(title = "Log2 Fold Change", targets = 2),
                list(title = "Adjusted P-value", targets = 3),
                list(title = "Base Mean", targets = 4)
            ))
            
            return(dt)
            
        }, error = function(e) {
            return(DT::datatable(data.frame(Error = paste("Error creating DEG table:", e$message))))
        })
    })

    # PCA implementation - CORRECTED VERSION
    output$pca_plot <- renderPlotly({
        tryCatch({
            # Load normalized counts for the selected contrast
            contrast_name <- input$contrast_filter
            norm_file <- paste0("../results/04_differential_expression/", contrast_name, "_normalized_counts.csv")
            
            if (file.exists(norm_file)) {
                norm_counts <- read.csv(norm_file, row.names = 1, check.names = FALSE)
                
                # Filter for most variable genes (top 1000)
                gene_vars <- apply(norm_counts, 1, var, na.rm = TRUE)
                top_var_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(1000, length(gene_vars))]
                
                # Subset to top variable genes
                pca_data <- norm_counts[top_var_genes, ]
                
                # Remove any genes with missing values
                pca_data <- pca_data[complete.cases(pca_data), ]
                
                if (nrow(pca_data) > 10 && ncol(pca_data) > 3) {
                    # Perform PCA (samples as rows, genes as columns)
                    pca_result <- prcomp(t(pca_data), scale. = TRUE)
                    
                    # Create PCA data frame
                    pca_df <- data.frame(
                        PC1 = pca_result$x[, 1],
                        PC2 = pca_result$x[, 2],
                        Sample = rownames(pca_result$x)
                    )
                    
                    # Add condition information based on the contrast
                    if (nrow(metadata_data) > 0) {
                        # Determine which condition column to use based on contrast
                        if (contrast_name == "post_vs_pre_training") {
                            # For post vs pre training contrast, use contrast2 column
                            if ("contrast2" %in% colnames(metadata_data)) {
                                pca_df <- merge(pca_df, metadata_data[, c("sample_id", "contrast2")], 
                                               by.x = "Sample", by.y = "sample_id", all.x = TRUE)
                                pca_df$condition <- pca_df$contrast2
                                pca_df$contrast2 <- NULL  # Remove the extra column
                            } else {
                                pca_df$condition <- "Unknown"
                            }
                        } else {
                            # For parkinsons vs control, use the regular condition column
                            pca_df <- merge(pca_df, metadata_data[, c("sample_id", "condition")], 
                                           by.x = "Sample", by.y = "sample_id", all.x = TRUE)
                        }
                        
                        # Handle missing values
                        pca_df$condition[is.na(pca_df$condition)] <- "Unknown"
                        
                        # Clean up condition names for better display
                        pca_df$condition <- gsub("_", " ", pca_df$condition)  # Replace underscores with spaces
                        pca_df$condition <- tools::toTitleCase(pca_df$condition)  # Capitalize properly
                        
                    } else {
                        pca_df$condition <- "Unknown"
                    }
                    
                    # Calculate variance explained
                    var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
                    
                    # Create PCA plot with contrast-specific colors
                    if (contrast_name == "post_vs_pre_training") {
                        color_values <- c("Pre Training" = "#1f77b4", "Post Training" = "#ff7f0e")
                    } else {
                        color_values <- c("Parkinsons" = "#d62728", "Control" = "#2ca02c")
                    }
                    
                    # Create PCA plot
                    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, text = Sample)) +
                        geom_point(size = 3, alpha = 0.7) +
                        scale_color_manual(values = color_values, name = "Group") +
                        labs(
                            x = paste0("PC1 (", var_explained[1], "% variance)"),
                            y = paste0("PC2 (", var_explained[2], "% variance)"),
                            title = paste("PCA Plot -", gsub("_", " ", contrast_name)),
                            color = "Group"
                        ) +
                        theme_minimal() +
                        theme(legend.position = "right")
                    
                    return(ggplotly(p, tooltip = c("text", "x", "y")))
                } else {
                    p <- ggplot() + 
                        annotate("text", x = 0, y = 0, label = "Insufficient data for PCA analysis") +
                        theme_minimal()
                    return(ggplotly(p))
                }
            } else {
                p <- ggplot() + 
                    annotate("text", x = 0, y = 0, label = paste("Normalized counts file not found:", norm_file)) +
                    theme_minimal()
                return(ggplotly(p))
            }
        }, error = function(e) {
            p <- ggplot() + 
                annotate("text", x = 0, y = 0, label = paste("Error creating PCA:", e$message)) +
                theme_minimal()
            return(ggplotly(p))
        })
    })
    
  #heatmap implementation:
  output$heatmap_plot <- renderPlotly({
      tryCatch({
          # Add gene symbol conversion function
          convert_gene_ids_to_symbols <- function(gene_ids) {
              tryCatch({
                  library(clusterProfiler)
                  library(org.Hs.eg.db)
                  conversion <- clusterProfiler::bitr(gene_ids, 
                                                fromType = "ENSEMBL", 
                                                toType = "SYMBOL", 
                                                OrgDb = org.Hs.eg.db)
                  
                  # Create a named vector for easy lookup
                  symbol_map <- setNames(conversion$SYMBOL, conversion$ENSEMBL)
                  
                  # Replace IDs with symbols, keep original if no symbol found
                  result <- ifelse(gene_ids %in% names(symbol_map), 
                                  symbol_map[gene_ids], 
                                  gene_ids)
                  
                  return(result)
              }, error = function(e) {
                  cat("Gene symbol conversion failed, keeping original IDs\n")
                  return(gene_ids)
              })
          }
          
          # Get top DE genes for the selected contrast
          data <- filtered_data()
          
          if (nrow(data) == 0) {
              p <- ggplot() + 
                  annotate("text", x = 0, y = 0, label = "No data available for selected contrast") +
                  theme_minimal()
              return(ggplotly(p))
          }
          
          # Get top 30 most significant genes
          top_genes <- data %>%
              filter(padj < input$padj_threshold, abs(log2FoldChange) > input$lfc_threshold) %>%
              arrange(padj) %>%
              head(30)
          
          if (nrow(top_genes) < 3) {
              p <- ggplot() + 
                  annotate("text", x = 0, y = 0, 
                          label = paste("Not enough significant genes for heatmap (found", nrow(top_genes), ")")) +
                  theme_minimal()
              return(ggplotly(p))
          }
          
          # Load normalized counts
          contrast_name <- input$contrast_filter
          norm_file <- paste0("../results/04_differential_expression/", contrast_name, "_normalized_counts.csv")
          
          if (file.exists(norm_file)) {
              norm_counts <- read.csv(norm_file, row.names = 1, check.names = FALSE)
              
              # Get counts for top genes
              gene_ids <- top_genes$gene
              heatmap_data <- norm_counts[intersect(gene_ids, rownames(norm_counts)), , drop = FALSE]
              
              if (nrow(heatmap_data) > 0) {
                  # Log transform (use log2(x+1) to match create_visualizations.R)
                  log_data <- log2(heatmap_data + 1)
                  
                  # Convert gene IDs to symbols BEFORE scaling
                  gene_symbols <- convert_gene_ids_to_symbols(rownames(log_data))
                  rownames(log_data) <- gene_symbols
                  
                  # Convert sample IDs to subject IDs if available
                  if ("subject_id" %in% colnames(metadata_data) && "sample_id" %in% colnames(metadata_data)) {
                      sample_to_subject <- setNames(metadata_data$subject_id, metadata_data$sample_id)
                      
                      # Update column names in heatmap data
                      valid_samples <- intersect(colnames(log_data), names(sample_to_subject))
                      if (length(valid_samples) > 0) {
                          colnames(log_data)[colnames(log_data) %in% valid_samples] <- 
                              sample_to_subject[valid_samples]
                          cat("Converted", length(valid_samples), "sample IDs to subject IDs\n")
                      }
                  }
                  
                  # Prepare data for ggplot (using log-transformed data, not scaled)
                  heatmap_long <- expand.grid(
                      Gene = rownames(log_data),
                      Sample = colnames(log_data)
                  )
                  heatmap_long$Value <- as.vector(as.matrix(log_data))
                  
                  # Create heatmap (matching create_visualizations.R style)
                  p <- ggplot(heatmap_long, aes(x = Sample, y = Gene, fill = Value)) +
                      geom_tile() +
                      scale_fill_gradient2(
                          low = "blue", mid = "white", high = "red", 
                          name = "log2(counts+1)"
                      ) +
                      theme_minimal() +
                      theme(
                          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                          axis.text.y = element_text(size = 8),
                          axis.title = element_text(size = 10)
                      ) +
                      labs(
                          x = "Subjects", 
                          y = "Gene Symbols",
                          title = paste("Top", nrow(heatmap_data), "DE Genes -", contrast_name)
                      )
                  
                  return(ggplotly(p))
              } else {
                  p <- ggplot() + 
                      annotate("text", x = 0, y = 0, label = "No matching genes found in normalized counts") +
                      theme_minimal()
                  return(ggplotly(p))
              }
          } else {
              p <- ggplot() + 
                  annotate("text", x = 0, y = 0, label = paste("Normalized counts file not found:", norm_file)) +
                  theme_minimal()
              return(ggplotly(p))
          }
      }, error = function(e) {
          p <- ggplot() + 
              annotate("text", x = 0, y = 0, label = paste("Error creating heatmap:", e$message)) +
              theme_minimal()
          return(ggplotly(p))
      })
  })
    
    # Pathway analysis table
        # Enhanced Pathway analysis table with directional results
        # Enhanced Pathway analysis table with directional results - CORRECTED VERSION
    output$pathway_table <- DT::renderDataTable({
        req(input$contrast_filter)
        
        tryCatch({
            # Try to load directional pathway results first
            directional_kegg_file <- paste0("../results/05_pathway_analysis/", input$contrast_filter, "_kegg_directional.csv")
            regular_kegg_file <- paste0("../results/05_pathway_analysis/", input$contrast_filter, "_kegg_enrichment.csv")
            
            pathway_data <- NULL
            analysis_type <- "Regular"
            
            # Prioritize directional results if available
            if (file.exists(directional_kegg_file)) {
                pathway_data <- read.csv(directional_kegg_file)
                analysis_type <- "Directional"
                cat("Loaded directional pathway data with", nrow(pathway_data), "pathways\n")
            } else if (file.exists(regular_kegg_file)) {
                pathway_data <- read.csv(regular_kegg_file)
                analysis_type <- "Regular"
                cat("Loaded regular pathway data with", nrow(pathway_data), "pathways\n")
            }
            
            if (!is.null(pathway_data) && nrow(pathway_data) > 0) {
                # Prepare display data based on analysis type
                if (analysis_type == "Directional" && "Direction" %in% colnames(pathway_data)) {
                    # Enhanced display with direction information - keep original column names
                    display_data <- pathway_data %>%
                        dplyr::select(Direction, Description, pvalue, p.adjust, Count) %>%
                        dplyr::arrange(Direction, p.adjust) %>%
                        head(25)  # Show top 25 pathways
                    
                    # Create proper column names for display
                    display_data_formatted <- data.frame(
                        `Direction` = display_data$Direction,
                        `Pathway Description` = display_data$Description,
                        `P-value` = round(display_data$pvalue, 4),
                        `Adjusted P-value` = round(display_data$p.adjust, 4),
                        `Gene Count` = display_data$Count,
                        stringsAsFactors = FALSE,
                        check.names = FALSE  # Keep spaces in column names
                    )
                    
                    # Create the data table with direction-based color coding
                    dt <- DT::datatable(
                        display_data_formatted,
                        options = list(
                            pageLength = 15, 
                            scrollX = TRUE,
                            columnDefs = list(
                                list(width = '120px', targets = 0),  # Direction column
                                list(width = '400px', targets = 1)   # Description column
                            )
                        ),
                        rownames = FALSE
                    ) %>%
                    DT::formatStyle(
                        "Direction",
                        backgroundColor = DT::styleEqual(
                            c("Up-regulated", "Down-regulated"), 
                            c("#ffcccc", "#ccccff")  # Light red for up, light blue for down
                        ),
                        fontWeight = "bold"
                    ) %>%
                    DT::formatStyle(
                        "Adjusted P-value",
                        backgroundColor = DT::styleInterval(
                            c(0.01, 0.05), 
                            c("lightcoral", "lightyellow", "lightgray")
                        )
                    )
                    
                    return(dt)
                    
                } else {
                    # Regular pathway display (fallback) - keep original column names
                    display_data <- pathway_data %>%
                        dplyr::select(Description, pvalue, p.adjust, Count) %>%
                        dplyr::arrange(p.adjust) %>%
                        head(20)
                    
                    # Create proper column names for display
                    display_data_formatted <- data.frame(
                        `Pathway Description` = display_data$Description,
                        `P-value` = round(display_data$pvalue, 4),
                        `Adjusted P-value` = round(display_data$p.adjust, 4),
                        `Gene Count` = display_data$Count,
                        stringsAsFactors = FALSE,
                        check.names = FALSE
                    )
                    
                    dt <- DT::datatable(
                        display_data_formatted,
                        options = list(pageLength = 15, scrollX = TRUE),
                        rownames = FALSE
                    ) %>%
                    DT::formatStyle(
                        "Adjusted P-value",
                        backgroundColor = DT::styleInterval(
                            c(0.01, 0.05), 
                            c("lightcoral", "lightyellow", "lightgray")
                        )
                    )
                    
                    return(dt)
                }
            } else {
                return(DT::datatable(data.frame(
                    Message = "No pathway enrichment results found",
                    Note = paste("Checked files:", directional_kegg_file, "and", regular_kegg_file)
                )))
            }
            
        }, error = function(e) {
            return(DT::datatable(data.frame(
                Error = paste("Error loading pathway data:", e$message)
            )))
        })
    })
    
    # Gene expression search
    output$gene_expr_result <- renderText({
        req(input$gene_search)
        
        if (input$gene_search == "") {
            return("Enter a gene ID to search for expression data.")
        }
        
        # Search for the gene in results
        gene_data <- results_data %>%
            filter(gene == input$gene_search | grepl(input$gene_search, gene, ignore.case = TRUE))
        
        if (nrow(gene_data) > 0) {
            result_text <- paste0("Found ", nrow(gene_data), " results for '", input$gene_search, "':\n")
            for (i in 1:min(5, nrow(gene_data))) {
                result_text <- paste0(result_text, 
                                     "Contrast: ", gene_data$contrast[i], 
                                     ", Log2FC: ", round(gene_data$log2FoldChange[i], 3),
                                     ", Adj.P: ", signif(gene_data$padj[i], 3), "\n")
            }
            return(result_text)
        } else {
            return(paste("No results found for gene:", input$gene_search))
        }
    })
}

# Run the app
shinyApp(ui = ui, server = server)