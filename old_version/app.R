# 

options(shiny.maxRequestSize = 100*1024^10)
options(repos = BiocManager::repositories())

library(shiny)
library(plotly)
library(data.table)
library(reshape2)
library(DT)
library(tidyverse)
library(shinydashboard)
library(DESeq2)
require(EnhancedVolcano)


# Start the dashboard ----

# Sidebar (inputs) ----
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItemOutput("load")),
  
  fileInput("file1", "Choose an R data-set",
            multiple = TRUE),
  
  fileInput("file2", "Choose CSV metadata",
            multiple = TRUE),
  
  # Input: Experimental desing to compare
  # Input an reactive selectiveInput later

  selectizeInput("selectedSample", label = "Samples to compare:",
                 choices = c('LOF_24DES','LOF_24POST',
                             'LOF_24PRE','LOF_30DES',
                             'LOF_30POST','LOF_30PRE'), 
                 multiple = TRUE, 
                 options = list(maxItems = 2)),
  
  actionButton("buttom", "Load!!"),
  
  helpText("Doble click to compare"),
  
  helpText("Play with the values above"), 
  
  # # 2. compare p-values distribution
  
  sliderInput('padj', 'Significance:', min = 0, max = 1,
              value = 0.05, step = NULL, round = 2),
  # 3. LogFC
  numericInput('logfc', 'log2(Fold Change)', min = 1, max = 8,
               value = 1, step = NA))


# The body of the dash ----
body <- dashboardBody(
  # Here the first box with tabs
  fluidRow(
    tabBox(width = 12,
      title = tagList(shiny::icon("gear"), "DE Visualization"),
      
      tabPanel("Pattern in expression",
               mainPanel(
                 plotly::plotlyOutput("areaPlot", width = "auto", height = "700px"))),
      
      tabPanel("Volcano Plot", 
               mainPanel(
                 plotOutput("volcano", width = "auto", height = "700px"))),
      
      tabPanel("P Values",  
               mainPanel(
                 plotly::plotlyOutput("phist", width = "auto", height = "700px"))),
      
      tabPanel("Significance genes",
               mainPanel(
                 plotly::plotlyOutput("sigdist", width = "auto", height = "700px")))
      
      ),
    fluidRow(
      infoBoxOutput("upgenesA"),
      infoBoxOutput("upgenesB"))
  )
)

# Put them together into a dashboardPage ----
ui <- dashboardPage(
  dashboardHeader(title = "Transcriptomic dashboard",
                  dropdownMenuOutput("messageMenu")),
  sidebar,
  body )


server <- function(input, output) {
  output$load <- renderMenu({
    menuItem("Load data", icon = icon("refresh"))
  })
  # dds table ----
  
  dds <-  eventReactive(input$buttom, {
    dds <- readRDS(file = input$file1$datapath)
  })
  
  contrast <- eventReactive(input$buttom, {
    
    sampleA <- input$selectedSample[1]
    sampleB <- input$selectedSample[2]
    
    contrast <- c("conditions", sampleA, sampleB)
    
    .res = results(dds(), contrast)
    
    baseMeanA <- rowMeans(counts(dds(), normalized=TRUE)[,colData(dds())$conditions == sampleA])
    baseMeanB <- rowMeans(counts(dds(),normalized=TRUE)[,colData(dds())$conditions == sampleB])
    
    
    res = cbind(baseMeanA, baseMeanB, as.data.frame(.res))
    
    res = cbind(round(select(res, -pvalue, -padj), 
                      digits = 2),
                select(res, pvalue, padj))
    
    res = cbind(sampleA=sampleA, sampleB=sampleB, as.data.frame(res))
    
    res$padj[is.na(res$padj)]  <- 1
    
    
    res <- data.table(res, ids = rownames(.res)) 
    
    # res <- res %>% rename(log2FC = log2FoldChange)
    
  })
  
  output$contrast <- DT::renderDataTable(
    DT::datatable(
      contrast()))
  
  # Metadata ----
  
  sam <- eventReactive(input$buttom, {    
    
    sam <- read.table(input$file2$datapath,
                      header = TRUE,
                      check.names = F, 
                      fill=T)
    sam = sam[sam[,2] != '',]
    # colnames(sam) = c('sample_name', 'replicate_name', 'f1','f2')
    
  })
  
  output$sam <- DT::renderDataTable(
    DT::datatable(sam(), 
                  options = list(searching = FALSE)))
  
  # histogram ----
  
  output$phist <- renderPlotly({
    
    col <- c("#899FE7", "#de2d26")
    
    hist <- contrast() %>%
      mutate(padjCol = ifelse(padj <= input$padj, "Used", "Dropped")) %>%
      filter(pvalue < 1) %>% 
      plot_ly(x = ~pvalue, 
              color = ~padjCol, colors = col,
              type = "histogram",
              cumulative = list(enabled=TRUE)) %>%
      layout(barmode = "overlay", 
             xaxis = list(title = "<b>P</b>values"),
             yaxis = list(title = "Cumulative Counts"))
    
    plotly_build(hist)
  })
  
  # Significanse gene distribution ---- 
  
  output$sigdist <- renderPlotly({
    
    signif <- contrast()[contrast()$padj <= input$padj]
    
    if(nrow(signif) > 500) {
      
      p <- contrast() %>%
        mutate(padjCol = ifelse(padj <= input$padj, "Used", "Dropped")) %>%
        mutate(logpvalue = -log10(pvalue)) %>%
        mutate(log2FoldChange = abs(log2FoldChange)) %>%
        sample_n(500) %>%
        plot_ly(x = ~ logpvalue, y = ~padj, color = ~log2FoldChange,
                type = 'scatter',
                mode = 'markers', 
                symbol = ~padjCol, symbols = c('circle','x'),
                text = ~ids,
                colors = 'YlGnBu',
                marker = list(size = 10, opacity = 0.7)) %>%
        layout(xaxis = list(title = "-Log<sub>10</sub><b>P</b>"),
               yaxis = list(title = "<b>P</b><i>adjust</i>")) %>%
        colorbar(title = "Log<sub>2</sub> Fold Change")
      
      plotly_build(p) 
    } else {
      
      p <- contrast() %>%
        mutate(padjCol = ifelse(padj <= input$padj, "Used", "Dropped")) %>%
        mutate(logpvalue = -log10(pvalue)) %>%
        mutate(log2FoldChange = abs(log2FoldChange)) %>%
        filter(padj <= input$padj) %>%
        plot_ly(x = ~ logpvalue, y = ~padj, color = ~log2FoldChange,
                type = 'scatter',
                mode = 'markers', 
                symbol = ~padjCol, symbols = c('circle','x'),
                text = ~ids,
                colors = 'YlGn',
                marker = list(size = 10, opacity = 0.7)) %>%
        layout(xaxis = list(title = "-Log<sub>10</sub><b>P</b>"),
               yaxis = list(title = "<b>P</b><i>adjust</i>")) %>%
        colorbar(title = "Log<sub>2</sub> Fold Change")
      
      plotly_build(p)
      
    }
    
  })
  
  # Area plot ----
  
  signif_ordered <- reactive({
    
    signif <- contrast()[contrast()$padj <= input$padj]
    
    genes_A <- signif[signif$log2FoldChange >= input$logfc]
    genes_A <- genes_A[order(abs(genes_A$log2FoldChange), decreasing = TRUE), ]
    
    genes_B <- signif[signif$log2FoldChange <= -input$logfc]
    genes_B <- genes_B[order(abs(genes_B$log2FoldChange), decreasing = FALSE), ]
    
    
    signif_ordered <- rbind(genes_A, genes_B)
    
  })
  
  normalized_sig_counts <- reactive({
    
    signifGenes <- unique(signif_ordered()$ids)
    
    counts(dds(), normalized = TRUE)[signifGenes, ] %>%
      reshape2::melt() %>%
      inner_join(sam(), by = c("Var2"= "names")) %>%
      data.frame() -> join_normalized_sig_counts
    
    # mean the replicates
    aggregate(join_normalized_sig_counts$value, 
              by = list(join_normalized_sig_counts$Var1, 
                        join_normalized_sig_counts$conditions), 
              mean) ->  x
    
    x <- data.table(genes = factor(x[,1], levels = signifGenes), 
                    conditions = x[,2], 
                    value = x[,3])
    
    
  })
  
  output$areaPlot <- renderPlotly({
    
    sample_color <- viridis::viridis(6)
    
    areaPlot <- normalized_sig_counts() %>%
      
      mutate(nsample = as.integer(factor(conditions))) %>% 
      
      plot_ly(x= ~genes, y = ~ value, 
              color = ~ conditions,
              fill = 'tonexty',
              mode = 'none',
              colors = sample_color,
              yaxis = ~ paste0("y", nsample)) %>%
      add_lines() %>%
      subplot(nrows = length(sample_color),
              shareX = TRUE)
    
    f1 <- list(
      family = "Arial, sans-serif",
      size = 7)
    
    z1 <- list(
      showticklabels = TRUE,
      tickangle = 45,
      tickfont = f1)
    
    areaPlot <- areaPlot %>% layout(xaxis = z1,
                                    xaxis = list(title = "cpm"))      
    
    plotly_build(areaPlot)
  })
  
  # volcano plot ---- 
  
  require(EnhancedVolcano) 
  
  output$volcano <- renderPlot({
    
    pCutoff <- round(max(contrast()[contrast()$padj <= 0.05, 'pvalue']), digits = 7)          
    
    # sampleA <- unique(contrast()$sampleA)
    # sampleB <- unique(contrast()$sampleB)
    # title <- paste0(sampleA, '_vs_', sampleB)
    
    
    v <- EnhancedVolcano(contrast(),
                         lab = contrast()$ids,
                         pCutoff = pCutoff,
                         FCcutoff = input$logfc,
                         boxedLabels = FALSE,
                         colAlpha = 3/5,
                         x = 'log2FoldChange',
                         y = 'pvalue',
                         title = '',
                         legendPosition = 'top',
                         subtitle = '',
                         caption = '')
    
    print(v)
  })
  
  # Boxes ---- 
  # output$upgenesA <- renderValueBox({
  #   sampleA <- unique(contrast()$sampleA)
  #   caption <- paste0("Up-genes in sample ", sampleA)
  #   
  #   signif <- contrast()[contrast()$padj <= input$padj]
  #   genes_A <- signif[signif$log2FoldChange >= input$logfc]
  #   valueBox(nrow(genes_A), icon = NULL, subtitle = caption)
  # })
  # 
  # output$upgenesB <- renderValueBox({
  #   sampleB <- unique(contrast()$sampleB)
  #   caption <- paste0("Up-genes in sample ", sampleB)
  #   
  #   signif <- contrast()[contrast()$padj <= input$padj]
  #   genes_B <- signif[signif$log2FoldChange <= -input$logfc]
  #   valueBox(nrow(genes_B), icon = NULL, subtitle = caption)
  # })
  
  # message menu istead of boxes
  output$messageMenu <- renderMenu({
    
    sampleA <- unique(contrast()$sampleA)
    sampleB <- unique(contrast()$sampleB)
    
    caption <- paste0("Up-genes in sample ", sampleA)
    caption <- paste0("Up-genes in sample ", sampleB)
    
    signif <- contrast()[contrast()$padj <= input$padj]
    
    n_genes_A <- nrow(signif[signif$log2FoldChange >= input$logfc])
    n_genes_B <- nrow(signif[signif$log2FoldChange <= -input$logfc])
    
    messageA <- paste0(n_genes_A ," ", caption)
    messageB <- paste0(n_genes_B ," ", caption)
    
    # pct_signif <- round(nrow(signif) / nrow(contrast()), 2)
    
    messageC <- paste0(nrow(signif) ," significance genes")
    
    messageA <- data.frame(from = sampleA, message = messageA)
    messageB <- data.frame(from = sampleB, message = messageB)
    messageC <- data.frame(from = 'Total', message = messageC)
    
    messageData <- rbind(messageA, messageB, messageC)

    msgs <- apply(messageData, 1, function(row) {
      messageItem(from = row[["from"]], message = row[["message"]])
    })
    dropdownMenu(type = "messages", .list = msgs)
  })
  
}

shinyApp(ui, server)




