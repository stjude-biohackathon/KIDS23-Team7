library(shiny)
library(shinyWidgets) # 
library(ggplot2) 
library(dplyr)  
library(data.table)
library(Seurat)
library(plotly)
source("~/KIDS2/GenericFunctions.R")
Sys.setenv(plotly_username = "plotly_name")


# Define UI for application that draws a histogram
ui <- fluidPage(

  # Application title
  titlePanel("KIDS23 TEAM7"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      
      ################
      # helper codes # 
      ################
      verbatimTextOutput("text_output"),   # instruction
      # to hide error ,,,
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"
      ),
      
      ##########
      # inputs #
      ##########
      radioButtons(inputId = "dataset",
                   label = "Sample numbers:",
                   choices = c("one sample", "multiple samples"),
                   selected = "''",
                   inline = TRUE,
      ),

      # uploading files 
      conditionalPanel(
        ## case1: one sample
        condition = "input.dataset == 'one sample'",
        fileInput("file1", "Choose matrix data (csv):",
                  multiple = TRUE,
                  buttonLabel = "Upload...",
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        fileInput("file2", "Choose meta data (csv):",
                  multiple = TRUE,
                  buttonLabel = "Upload...",
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv"))
      ),
      ## case2: multiple samples 
      conditionalPanel(
        condition = "input.dataset == 'multiple samples'",
        fileInput("file3", "Choose your File (zip):",
                  multiple = TRUE,
                  buttonLabel = "Upload...",
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv"))),
      
      ##############
      # assay type #
      ##############
      selectInput(inputId = "assaytype",
                  label = "Assay type:",
                  choices = c("Spatial(visium)", "Akoya(phenoCycler)", "Vizgen(MERSCOPE)", "Nanostring(CosMx)"),
                  selected = 'Akoya(phenoCycler)'),
      
      actionButton(inputId = "run_code", label = "Run your data...",)
      
    ),
    
    #############
    # mainPanel #
    #############
    mainPanel(
      textOutput("text"),
      tableOutput("table"),

      tabsetPanel(
        # tabPanel("Processing",
        #          verbatimTextOutput(c("obj","obj2")),
        #          # verbatimTextOutput("obj2"),
        #          verbatimTextOutput("obj3"),
        #          verbatimTextOutput("obj4"),
        #          ),
        tabPanel("QC",
          plotOutput(outputId = "qc_violin"),
          uiOutput("gene_dropdown"),
          plotOutput(outputId = "qc_histogram"),
          ),
        tabPanel("Analysis", plotOutput(outputId = "plot")),
        tabPanel("Visualization",  plotlyOutput(outputId = "cell_annotation",height = 800, width = 1200)),
        tabPanel("STAT", verbatimTextOutput(outputId = "stat"))
      )
    )
    
  )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
  options(shiny.maxRequestSize = 1000 * 1024^2) # increase max upload size to 1000 MB
  output$text_output <- renderText({"Instructions:\n1.choose your sample number\n2.upload your files\n3.run your data\n4.filtering cells..."})

  # 1. csv to Seurat object
  # seurat_object <- reactive({
  #   req(input$file1)
  #   req(input$file2)
  #   csv_to_seurat(ExpressionMarker = input$file1$datapath, MetaData = input$file2$datapath)
  # })
  # 
  # output$obj <- renderPrint({
  #   print(seurat_object())})

  # seurat_object <- reactive({
  #   if(input$dataset == "one sample") {
  #     if(!is.null(input$file1) & !is.null(input$file2)) {
  #       csv_to_seurat(ExpressionMarker = input$file1$datapath, MetaData = input$file2$datapath)
  #     }
  #   } else if(input$dataset == "multiple samples") {
  #     if(!is.null(input$file3)) {
  #       # code for multiple samples
  #     }
  #   }
  # }
  
  # output$obj <- renderPrint(
  #   print(seurat_object())
  # )
  
  # using RDS
  seurat_object <- reactive({
    # req(file.exists("/home/lead/Akoya/rds/Sample2_Group2_SeuratObj_sketch.rds"))
    readRDS("/home/ybae/KIDS2/test_annoated.rds")
  })
  

  # ui <- fluidPage(
  #   tags$head(tags$script(src = "message-handler.js")),
  #   actionButton("do", "Click Me")
  # )
  # 
  # server <- function(input, output, session) {
  #   observeEvent(input$do, {
  #     session$sendCustomMessage(type = 'testmessage',
  #                               message = 'Thank you for clicking')
  #   })
  # }
  
  # 2. QC # 
  observeEvent(input$run_code, 
               { output$qc_violin = renderPlot({
    AssayType <- 'Akoya'
    obj = seurat_object()
    Idents(obj) = "CelltypePrediction"
    VlnPlot(obj, features = c(paste0('nCount_', AssayType),
                                          paste0('nFeature_', AssayType)), pt.size = 0)
  })

  # dropdown genes
  output$gene_dropdown <- renderUI({
    req(seurat_object())
    genes <- rownames(seurat_object())
    selectInput(inputId = "gene", label = "Select Gene", choices = genes)
  })
  
  output$qc_histogram = renderPlot({
    obj = seurat_object()
    Idents(obj) = "CelltypePrediction"
    RidgePlot(obj, features = input$gene, ncol = 2, layer = 'counts') & xlab("")
  })})
  
  ### since we are using rds ### 
  # 3. normalization and sketch  # 
  # skected <- reactive({
  #   norm_and_sketch(seurat_object())
  # })
  
  # to checking
  # output$obj2 <- renderPrint({
  #   print(skected())
  # })
  
  # 4. seurat_workflow # 
  # processed_seurat <- reactive({
  #   (seurat_workflow(skected()))
  # })
  
  # to checking 
  # output$obj3 <- renderPrint({
  #   print(processed_seurat())
  # })
  
  # 5. annotation # 
  # annotate_seurat <- reactive({
  #   (annotate_cell_types())
  # })
  
  # to checking 
  # output$obj4 <- renderPrint({
  #   print(annotate_seurat())
  # })
  
  output$cell_annotation <- renderPlotly({
    plot <- ImageDimPlot(seurat_object(), fov = c('Sample2_Group2'), size = 0.3, group.by = 'CelltypePrediction', coord.fixed = T, dark.background = F)
    ggplotly(plot)
    
  })
  
  

}


# Run the application 
shinyApp(ui = ui, server = server)
