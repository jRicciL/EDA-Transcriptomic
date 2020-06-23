
library(shiny)

# options(repos = BiocManager::repositories()

# Start the dashboard ----

# Sidebar (inputs) ----

sidebar <- dashboardSidebar(
  sidebarMenu(
    id = "sider",
    menuItemOutput("load")),
  
  fileInput("file1", "Choose a DE-result set",
            multiple = TRUE),
  # fileInput("file2", "Choose CSV metadata",
  #           multiple = TRUE),
  fileInput("file3", "Choose count-data",
            multiple = TRUE),
  
  textAreaInput(
    "groupSelectViaText",
    "Input your group info",
    rows = 6,
    placeholder = paste(
      "Please input group information at here. Here is a example format:",
      "-----",
      "G1_rep1,Group1",
      "G1_rep2,Group1",
      "G1_rep3,Group1",
      sep = '\n'
    )
  ),
  
  actionButton("buttom", "Load!!"),
  
  helpText("Doble click to compare"),
  
  # Input: Experimental desing to compare
  # Input an reactive selectiveInput later
  
  selectizeInput("selectedSample", label = "Samples to compare:",
                 choices = c('LOF_24DES','LOF_24POST',
                             'LOF_24PRE','LOF_30DES',
                             'LOF_30POST','LOF_30PRE'), 
                 multiple = TRUE, 
                 options = list(maxItems = 2)),
  
  
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


# test 
# in the ui
# Display this only if the density is shown
# conditionalPanel(condition = "input.density == true",
#                  sliderInput(inputId = "bw_adjust",
#                              label = "Bandwidth adjustment:",
#                              min = 0.2, max = 2, value = 1, step = 0.2)
# )
# 
# # in the server
# output$main_plot <- reactivePlot(width = 400, height = 300, function() {
#   
#   hist(faithful$eruptions,
#        probability = TRUE,
#        breaks = as.numeric(input$n_breaks),
#        xlab = "Duration (minutes)",
#        main = "Geyser eruption duration")
#   
#   if (input$individual_obs) {
#     rug(faithful$eruptions)
#   }
#   
#   if (input$density) {
#     dens <- density(faithful$eruptions, adjust = input$bw_adjust)
#     lines(dens, col = "blue")
#   }
#   
# })
# })



# Put them together into a dashboardPage ----
# 

ui <- dashboardPage(
  dashboardHeader(title = "Transcriptomic dashboard",
                  dropdownMenuOutput("messageMenu")),
  sidebar,
  body )


shinyApp(ui, server)



