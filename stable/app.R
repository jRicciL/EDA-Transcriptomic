
#
options(shiny.maxRequestSize = 100*1024^10, 
        stringsAsFactors = FALSE,
        repos = BiocManager::repositories())


library(shiny)
library(plotly)
library(DESeq2)
library(tidyverse)
library(DT)
library(data.table)
library(heatmaply)

# UI SIDE-BAR ----

sidebarPanel <- sidebarPanel(
  titlePanel("trans-EDA"),

  fileInput("res", "Load the DE-result", multiple = TRUE, 
            accept=c('.DE_results', '.rds')),
  fileInput("count", "Load the matrix-count", multiple = TRUE, 
            accept=c('.matrix', '.csv')),
  
  #selectInput('sample', 'Select Sample:', choices = c('s1', 's2')),
  
  textAreaInput("groupSelectViaText",
    "Input your group info",
    rows = 6, placeholder = paste(
      "Please input group information at here.",
                                  "-----", 
      "G1_rep1,Group1", "<enter-separated>G2_rep1,Group2", sep = "\n")),
  
  actionButton("buttom", "Load!!"),
  
  h2(),
  
  sliderInput('padj', 'Significance:', min = 0.01, max = 1,
              value = 0.95, step = NULL, round = 2),

  numericInput('logfc', 'log2(Fold Change)', min = 1, max = 8,
               value = 2, step = NA),
  )
  

# UI - Render server to MainPanel ----
tabsetPanel <- tabsetPanel(
  tabPanel('Differential Expression significance', 
           plotly::plotlyOutput("phist", width = "auto", height = "400px"),
           plotly::plotlyOutput("pdist", width = "auto", height = "400px")),
  
  tabPanel('Tables',
           plotly::plotlyOutput('areaPlot', width = "auto", height = "400px"),
           DT::dataTableOutput('signif_ordered')),
  
  tabPanel('VolcanoPlot',
           plotly::plotlyOutput('volcano', width = "auto", height = "400px")),
  
  tabPanel('Heatmap',
           plotly::plotlyOutput('heatmap', width = "auto", height = "400px"))
  )

mainPanel <- mainPanel(
  tabsetPanel)

# UI-BUILD LAYOUT ----
# uidb <- dashboardPage(
#   dashboardHeader(title = "Transcriptomic dashboard",
#                   dropdownMenuOutput("messageMenu")),
#   body = dashboardBody(tabsetPanel),
#   sidebar = dashboardSidebar(sidebarPanel, width = 350) )

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
  
  #****************************
  # Handle null or empty inputs
  #****************************
  # observe for not null inputs 
  observe({
    input$buttom # Take a dependence on button
    groupinfo_load <- input$groupSelectViaText == ''
    res_load <- is.null(input$res)
    count_load <- is.null(input$count)
    
    # If any imput is null  show info
    if(res_load){
      showNotification('Please upload your DE-result')
    } else if(count_load) {
      showNotification('Please upload your matrix coounts file.')
    } else if(count_load) {
      showNotification('Please paste your group information in the left-panel box.')
    }
  })

  # req to handle null inputs inside reactive functions
  check_empty_inputs <- function() { req(input$res, input$count, input$groupSelectViaText) }
  #****************************
  
  # TAB-1 ----
  sam <- eventReactive(input$buttom, {
    check_empty_inputs()
    group_info <- input$groupSelectViaText
    sam <- fread(group_info, header = FALSE, col.names = c('name','group'))
  })
  
  loadResults <- function(datapath) {
    
    read.csv(datapath, header=T, sep = '\t') %>%
      as_tibble(rownames = "ids") %>%
      mutate_at(vars(!matches("ids|sample|pvalue|padj")),
                round ,digits = 2)
  }
  
  DE <- eventReactive(input$buttom, {
    check_empty_inputs()
    samples <- names(table(sam()$group))
    sampleA <- samples[1]
    sampleB <- samples[2]
    
    DE <- loadResults(input$res$datapath) %>%
      filter(sampleA == sampleA & sampleB == sampleB)
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
  
  # not rendered yet
  output$sam <- DT::renderDataTable(
    DT::datatable(
      sam()))
  
  output$signif_ordered <- DT::renderDataTable(
    DT::datatable(
      signif_ordered(),
      extensions = 'Buttons', 
      options = list(
        pageLength = 25,  
        dom = 'Bfrtip',
        buttons = 
          list('copy', 'print', list(
            extend = 'collection',
            buttons = c('csv', 'excel'),
            text = 'Download'
          )))
      ))
  
  #********************************
  # VOLCANO PLOT
  #********************************
  
  # Helper funtions to draw vertical an horizontal lines with plotly
  vline <- function(x = 0, ... ){
    list(type='line', y0=0, y1=1, yref='paper',
         x0 = x, x1 = x,
         ...)}
  hline <- function(y = 0, ... ){
    list(type='line', x0=0, x1=1, xref='paper',
         y0=y, y1=y,
         ...)}
  # Default dashed line format
  dashed_line = list(color='black', dash='dot', widt = 0.8)
  
  #********************************
  # Reactive data filter
  get_X_volcano <- reactive({
    
    # Input values
    pCutoff <- input$padj
    nlog_pCutoff <- -log10(pCutoff)
    fcCutoff <- input$logfc
    
    # Get DE table
    X <- DE() %>% na.omit() %>% 
      # TODO: Only work with padj? Why combine both?
      mutate(
        neg_log10_p = -log10(pvalue)
      ) %>%
      mutate(signif_col  =
              ifelse(neg_log10_p >= nlog_pCutoff & abs(log2FoldChange) >= fcCutoff, sig_labels[1],
              ifelse(neg_log10_p <  nlog_pCutoff & abs(log2FoldChange) >= fcCutoff, sig_labels[2],
              ifelse(neg_log10_p >= nlog_pCutoff & abs(log2FoldChange) <  fcCutoff, sig_labels[3],
                  sig_labels[4]
          )))
      ) %>%
      mutate(signif_col = factor(signif_col, levels = sig_labels)) 
    # Subset some of the non-significative genes to ommit overload the plot
    non_sig <- filter(X, signif_col == 'NS')
    if (nrow(non_sig) > 1000){
      X <- rbind(
        # Subset non significative values
        filter(X, signif_col != 'NS'),
        # Keep a third of NS values
        non_sig %>% sample_n(nrow(non_sig) %/% 3)
      )
    }
    
    return(X)
  })
  
  #********************************
  output$volcano <-  renderPlotly({
    
    # Input values
    pCutoff <- input$padj
    nlog_pCutoff <- -log10(pCutoff)
    fcCutoff <- input$logfc
    
    # Get the filtered data
    X <- get_X_volcano()
    
    # Sig. Labels
    sig_labels = c('p-value and log<sub>2</sub>FC', 'p-value', 'log<sub>2</sub>FC', 'NS')
    volcano_colors = c('#DC0D0D',  '#27A871', '#0D8DB0', '#80847C')
    
    # Create the Volcano plot unisng plotly
    volPlotly <- X %>% 
      plot_ly() %>%
      add_trace(type='scatter', mode='markers',
                x = ~log2FoldChange, y = ~neg_log10_p, 
                color = ~signif_col, colors=volcano_colors, alpha=0.7,
                hoverinfo='none') %>%
      layout(shapes = list(vline(x = fcCutoff, line = dashed_line), 
                           vline(x = -fcCutoff, line = dashed_line), 
                           hline(y = nlog_pCutoff, line = dashed_line)))
    # Add Significative genes as a new trace 
    X_sig <- X %>% filter(signif_col == sig_labels[1])
    volPlotly <- volPlotly %>%  
      add_trace(type='scatter', mode='markers',
                x = X_sig$log2FoldChange, y = X_sig$neg_log10_p,
                hoverinfo = 'text', text = X_sig$ids, 
                showlegend=FALSE,
                marker=list(color=volcano_colors[1], 
                            size=7,
                            line=list(color = 'black',
                                      width=0.5))) %>%
      layout( yaxis = list(title = "-log<sub>10</sub>p"),
              xaxis = list(title = "log<sub>2</sub> Fold Change", zeroline=FALSE),
              legend = list(x = 0.45, y = -0.2,
                            orientation='h', xanchor='center'))
    
    volPlotly
    
  })
  #*********************************
  
  
  # message menu istead of boxes ----
  # output$messageMenu <- renderMenu({
    observeEvent(input$buttom, {  
    samples <- names(table(sam()$group))
    sampleA <- samples[1]
    sampleB <- samples[2]
    
    caption <- paste0("Up-genes in sample ", sampleA)
    caption <- paste0("Up-genes in sample ", sampleB)
    
    signif <- DE() %>%
      filter(padj <= input$padj)
    
    n_genes_A <- signif %>% filter(log2FoldChange >= input$logfc) %>% distinct(ids) %>% nrow
    n_genes_B <- signif %>% filter(log2FoldChange <= -input$logfc) %>% distinct(ids) %>% nrow
    
    messageA <- paste0(n_genes_A ," ", caption)
    messageB <- paste0(n_genes_B ," ", caption)
    

    messageC <- paste0(nrow(signif) ," significance genes")
    
    messageA <- data.frame(from = sampleA, message = messageA)
    messageB <- data.frame(from = sampleB, message = messageB)
    messageC <- data.frame(from = 'Total', message = messageC)
    
    messageData <- rbind(messageA, messageB, messageC)
    
    msgs <- apply(messageData, 1, function(row) {
      messageItem(from = row[["from"]], message = row[["message"]])
    })
    
    #dropdownMenu(type = "messages", .list = msgs)
    showNotification(msgs, type = "message")
  })
  
  # Heatmap ----
  
  output$heatmap <-  renderPlotly({
    
    minl <- min(normalized_sig_counts()$value)
    maxl <- max(normalized_sig_counts()$value)
    
    # row_side_colors <- colData(dds)
    # row_side_colors <- sam[, sam$group]
 
    normalized_sig_counts() %>%
      pivot_wider(names_from = group, values_from = value,  
                  values_fill = list(value = 0)) %>%
      data.frame(row.names = 1) %>% 
      heatmaply(
        seriate = "mean", 
        showticklabels = c(TRUE, FALSE),
        limits = c(minl, maxl),
        plot_method = "plotly",
        k_row = 3, k_col = 2
        # row_side_colors = row_side_colors
      )
  })
  
  
    

}

shinyApp(ui = ui, server = server)

# sam <- fread(c("LOF_24DES_1,LOF_24DES\nLOF_24DES_2,LOF_24DES\nLOF_24DES_3,LOF_24DES\nLOF_24POST_1,LOF_24POST\nLOF_24POST_2,LOF_24POST\nLOF_24POST_3,LOF_24POST"), header = FALSE, col.names = c('name','group'))


