
#
options(shiny.maxRequestSize = 100*1024^10, stringsAsFactors = FALSE)
#

library(shiny)
library(plotly)
library(DESeq2)
library(tidyverse)
library(DT)
library(data.table)
# UI SIDE-BAR ----

sidebarPanel <- sidebarPanel(
  titlePanel("trans-EDA"),

  fileInput("res", "Load the DE-result", multiple = TRUE),
  fileInput("count", "Load the matrix-count", multiple = TRUE),
  
  #selectInput('sample', 'Select Sample:', choices = c('s1', 's2')),
  
  textAreaInput("groupSelectViaText",
    "Input your group info",
    rows = 6, placeholder = paste(
      "Please input group information at here.",
                                  "-----", 
      "G1_rep1,Group1", "<enter-separated>G2_rep1,Group2", sep = "\n")),
  
  actionButton("buttom", "Load!!"),
  
  h2(),
  
  sliderInput('padj', 'Significance:', min = 0, max = 1,
              value = 0.05, step = NULL, round = 2),

  numericInput('logfc', 'log2(Fold Change)', min = 1, max = 8,
               value = 1, step = NA))

# UI - Render server to MainPanel ----
tabsetPanel <- tabsetPanel(
  tabPanel('Differential Expression significance', 
           plotly::plotlyOutput("phist", width = "auto", height = "400px"),
           plotly::plotlyOutput("pdist", width = "auto", height = "400px")),
  tabPanel('Tables',
           plotly::plotlyOutput('areaPlot', width = "auto", height = "400px"),
           DT::dataTableOutput('signif_ordered')),
  
  tabPanel('VolcanoPlot',
           plotly::plotlyOutput('volcano', width = "auto", height = "400px"))
  )

mainPanel <- mainPanel(
  tabsetPanel)

# UI-BUILD LAYOUT ----

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel,
    mainPanel
  )


)

# SERVER SIDE ----

server <- function(input, output, session) {
  
  help_text <- c("trans-EDA is a simple Exploratory Data Analysis dashboard. Pleas load your Differential Expresion results and count-matrix to start")
  
  showModal(modalDialog(help_text))
  
  # TAB-1 ----
  loadDE <- function(.sampleA, .sampleB) {
    
    read.csv(input$res$datapath,header=T, sep = '\t') %>%
      as_tibble(rownames = "ids") %>%
      mutate_at(vars(!matches("ids|sample|pvalue|padj")),
                round ,digits = 2)
      # mutate(pvalue = formatC(pvalue, format = "E", digits = 2)) %>%
      # mutate(padj = formatC(padj, format = "E", digits = 2))
    
    # %>% filter(sampleA == .sampleA & sampleB == .sampleB)
  }
  
  DE <- eventReactive(input$buttom, {
    
    # sampleA <- loadDE$selectedSample[1]
    # sampleB <- loadDE$selectedSample[2]
    
    DE <- loadDE(sampleA, sampleB)
  })
  
  output$DE <- DT::renderDataTable(
    DT::datatable(
      DE()))

  output$phist <- renderPlotly({
    
    col <- c("#899FE7", "#de2d26")
    
    hist <- DE() %>%
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
  
  output$pdist <- renderPlotly({
    
    signif <- DE() %>% 
     mutate(padjCol = ifelse(padj <= input$padj, "Used", "Dropped")) %>%
     mutate(logpvalue = -log10(pvalue)) %>%
     mutate(log2FoldChange = abs(log2FoldChange)) %>% 
     filter(padj <= input$padj)
    
    if(nrow(signif) > 1000) {
      
      p <- signif %>%
        sample_n(1000) %>%
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
      
      p <- signif %>%
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
  
  # TAB-2 ----
  
  sam <- eventReactive(input$buttom, {
    sam <- fread(input$groupSelectViaText, header = FALSE, col.names = c('name','group'))
  })
  
  signif_ordered <- reactive({
    signif <- DE() %>% filter(padj <= input$padj)
    
    genes_A <- signif %>% filter(log2FoldChange >= input$logfc) %>%
      arrange(desc(abs(log2FoldChange)))
    
    genes_B <- signif %>% filter(log2FoldChange <= -input$logfc) %>%
      arrange(desc(abs(log2FoldChange)))
    
    signif_ordered <- rbind(genes_A, genes_B)
    
  })
  
  dds <- eventReactive(input$buttom, {    
    
    x <- read.table(input$count$datapath, header=T, com='')
    
    countData <- round(x)
    
    colData <- names(countData)
    
    colData <- data.frame(conditions=factor(colData), row.names = colData)
    
    dds <- DESeqDataSetFromMatrix(
      countData = countData,
      colData = colData,
      design = ~ conditions)
    
    dds <- estimateSizeFactors(dds)
  })
  
  normalized_sig_counts <- reactive({
    
    signifGenes <- unique(signif_ordered()$ids)
    
    counts(dds(), normalized = TRUE)[signifGenes, ] %>%
      as_tibble(rownames = 'ids') %>% 
      pivot_longer(-ids) %>%
      inner_join(sam()) -> x
    
    n <- length(unique(sam()$group))
    if(n > 1) {
      # mean the replicates
      
      normalized_sig_counts <- x %>%
        group_by(ids, group) %>%
        summarise(value = mean(value)) %>%
        ungroup()


    } else { normalized_sig_counts <- x }
    
  })
  
  output$areaPlot <- renderPlotly({
    
    n <- length(unique(sam()$group))
    sample_color <- viridis::viridis(n)
    
    f1 <- list(
      family = "Arial, sans-serif",
      size = 7)
    
    z1 <- list(
      showticklabels = TRUE,
      tickangle = 45,
      tickfont = f1)
    
    
    areaPlot <- normalized_sig_counts() %>%
      mutate(nsample = as.integer(factor(group))) %>%
      plot_ly(x = ~ids, 
              y = ~value, 
              color = ~group,
              fill = 'tonexty',
              colors = sample_color,
              yaxis = ~paste0("y", sort(nsample, decreasing =F))) %>%
      layout(
        xaxis = z1,
        yaxis = f1,
        hoverlabel = list(font=list(size=20))
      ) %>% 
      add_lines() %>%
      subplot(nrows = length(sample_color), shareX = TRUE)
    
    plotly_build(areaPlot)
    
# 
#     panel <- . %>% 
#       plot_ly(x = ~ids, y = ~value,
#               color = ~ group,
#               fill = 'tonexty',
#               mode = 'none',
#               colors = sample_color) %>%
#       add_lines(color = ~ group) %>%
#       add_annotations(
#         text = ~unique(group),
#         x = 0.5,
#         y = 1,
#         yref = "paper",
#         xref = "paper",
#         yanchor = "bottom",
#         showarrow = FALSE,
#         font = list(size = 15)
#       )  %>%
#       layout(
#         showlegend = F,
#         xaxis = z1,
#         xaxis = list(title = "Count Per million")
#       )
#     # 
    # areaPlot <- normalized_sig_counts() %>%
    #   group_by(group) %>%
    #   do(p = panel(.)) %>%
    #   subplot(nrows = NROW(.), shareX = TRUE)
    # 
    #  
    # 

  })
  
  output$sam <- DT::renderDataTable(
    DT::datatable(
      sam()))
  
  output$signif_ordered <- DT::renderDataTable(
    DT::datatable(
      signif_ordered()))
  
  # volcano 
  # signif <- DE() %>% filter(padj <= input$padj)
  
  output$volcano <-  renderPlotly({
    
    x <- filter(DE(), !is.na(padj))
    # x$gene_name <- ifelse(is.na(x$ids), x$ids, x$gene_name )
    
    x$ids[x$padj > input$padj & abs(x$log2FoldChange) < input$logfc] <- NA
    
    x <- filter(x, padj > input$padj & log2FoldChange < input$logfc)
    
    vplot <- plot_ly(data = x, x = ~ log2FoldChange , y = ~ -log10(padj),
            type = "scatter", mode = "markers", 
            hoverinfo = "text", text = ~ ids,
            marker = list(size = 10, color = 'rgba(0, 0, 255, .3)')) %>%
      layout( yaxis = list(title = "-Log10 p-value", zeroline = FALSE),
              xaxis = list(title = "Log2 fold change", zeroline = FALSE, range=c(-6,6)))
    
    plotly_build(vplot)
    
  })
  
  
  
  

}

shinyApp(ui = ui, server = server)

# sam <- fread(c("LOF_24DES_1,LOF_24DES\nLOF_24DES_2,LOF_24DES\nLOF_24DES_3,LOF_24DES\nLOF_24POST_1,LOF_24POST\nLOF_24POST_2,LOF_24POST\nLOF_24POST_3,LOF_24POST"), header = FALSE, col.names = c('name','group'))


