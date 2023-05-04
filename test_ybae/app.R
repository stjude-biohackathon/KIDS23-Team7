#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2) 
library(dplyr)  

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("KIDS23 TEAM7"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(
        
        # Input: Checkbox if file has header ----
        #checkboxInput("sample number", "one sample", TRUE),
        #checkboxInput("sample number", "multiple samples", TRUE),
       radioButtons(inputId = "dataset",
                    label = "Sample numbers:",
                    choices = c("one sample", "multiple samples"),
                    selected = "''"),
        
        # uploading files 
        fileInput("file1", "Choose CSV File",
                  multiple = TRUE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        
        # assay type 
        selectInput(inputId = "dataset",
                    label = "Assay type:",
                    choices = c("akoya", "visium", "vizgen"),
                    selected = '"'),
        
        sliderInput(inputId = "obs",
                    label = "Number of observations to show:",
                    min = 1, max = 10, value = 5)
      ),
      
      
      
      mainPanel(
        tabsetPanel(
          tabPanel("Table", tableOutput(outputId = "table")),
          tabPanel("Plot", plotOutput(outputId = "plot")),
          tabPanel("Summary", verbatimTextOutput(outputId = "summary"))
        )
      )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white',
             xlab = 'Waiting time to next eruption (in mins)',
             main = 'Histogram of waiting times')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
