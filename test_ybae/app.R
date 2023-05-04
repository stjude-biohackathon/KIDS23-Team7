library(shiny)
library(shinyWidgets) # 
library(ggplot2) 
library(dplyr)  
library(data.table)
library(Seurat)
source("generic_functions.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("KIDS23 TEAM7"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      
      verbatimTextOutput("text_output"),  
      
      # Input
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
      
      # assay type 
      selectInput(inputId = "assaytype",
                  label = "Assay type:",
                  choices = c("Spatial(visium)", "Akoya(phenoCycler)", "Vizgen(MERSCOPE)", "Nanostring(CosMx)"),
                  selected = ''),
    ),
    
    # Push to analysis
    mainPanel(
      textOutput("text"),
      tableOutput("table"),
      verbatimTextOutput("obj"),
      verbatimTextOutput("text2"),
      verbatimTextOutput("text3"),
      
      tabsetPanel(
        tabPanel("QC", 
                 
                 tableOutput(outputId = "qcplot1")),
        tabPanel("Analysis", plotOutput(outputId = "plot")),
        tabPanel("Visualization", verbatimTextOutput(outputId = "visualization"))
      )
    )
    
  )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
  options(shiny.maxRequestSize = 1000 * 1024^2) # increase max upload size to 1000 MB
  output$text_output <- renderText({"Instructions:\n1.choose your sample number\n2.upload your files\n3.run your QC\n4.filtering cells"})
  
  # Read CSV file (testing purpose)
  # suppressWarnings(data <- reactive({
  #   req(input$file1)
  #   read.csv(input$file1$datapath, header = TRUE)
  # }))
  # output$table <- renderTable({
  #   head(data())
  # })
  
  # using RDS
  # output$obj <- renderPrint(
  #   readRDS("/home/lead/Akoya/rds/Sample2_Group1_SeuratObj_sketch.rds")
  # )

  # csv to Seurat object
  seurat_object <- reactive({
    csv_to_seurat(ExpressionMarker = input$file1$datapath, MetaData = input$file2$datapath)})
  
  output$obj <- renderPrint(
    print(seurat_object())
  )
  
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
  
  # QC # 
  # qc_histogram <- reactive({RidgePlot(seurat_object, features = rownames(test)[1:10], ncol = 2, layer = 'counts') & xlab("")})
  # output$qcplot1 <- renderPlot({
  # qc_histogram()
  # })
  
  
  # normalization and sketch  # 
  skected <- reactive({
    norm_and_sketch(seurat_object())
  })
  
  # to checking
  output$text2 <- renderPrint({
    print(class(skected()))
  })
  
  # seurat_workflow # 
  processed_seurat <- reactive({
    (seurat_workflow(skected()))
  })
  
  # to checking 
  output$text3 <- renderPrint({
    print(processed_seurat())
  })

  
}


# Run the application 
shinyApp(ui = ui, server = server)
