library(shiny)
library(shinyWidgets) # 
library(ggplot2) 
library(dplyr)  
library(data.table)

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
          ## one sample
          condition = "input.dataset == 'one sample'",
          fileInput("file1", "Choose matrix data (csv):",
                    multiple = TRUE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          fileInput("file2", "Choose meta data (csv):",
                    multiple = TRUE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv"))
        ),
        
        ##  multiple samples 
        conditionalPanel(
          condition = "input.dataset == 'multiple samples'",
          fileInput("file3", "Choose your File:",
                    multiple = TRUE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv"))),
        
        # assay type 
        selectInput(inputId = "assaytype",
                    label = "Assay type:",
                    choices = c("Spatial(visium)", "Akoya(phenoCycler)", "Vizgen(MERSCOPE)", "Nanostring(CosMx)"),
                    selected = '"'),
      ),
      
      # Push to analysis
      
      mainPanel(
        tableOutput("table"),
        tabsetPanel(
          tabPanel("QC", tableOutput(outputId = "qcplot")),
          tabPanel("Analysis", plotOutput(outputId = "plot")),
          tabPanel("Visualization", verbatimTextOutput(outputId = "visualization"))
        )
      )
      
    )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
  options(shiny.maxRequestSize = 1000 * 1024^2) # increase max upload size to 1000 MB
  
  # message("Upload files and run QC first") # this gives message in your console
  
  output$text_output <- renderText("Instructions:\n1.choose your sample number\n2. upload your files\n3. run your QC\n4. filtering cells ")
  # Read CSV file
  data <- reactive({
    req(input$file1)
    fread(input$file1$datapath, header = TRUE, data.table=FALSE)
  })
  
  # Render table
  output$table <- renderTable({
    data()
  })
  
  
  # output$table <- renderTable({
  #   # input$file1 will be NULL initially. After the user selects
  #   # and uploads a file, head of that data file by default,
  #   # or all rows if selected, will be shown.
  #   req(input$file1) # to check if an argument or input is present and not NUL
  #   fread(input$file1$datapath,
  #                  header = TRUE,
  #                  sep = ",",
  #                  quote = "",
  #               data.table=FALSE)
  # })
}


# Run the application 
shinyApp(ui = ui, server = server)
