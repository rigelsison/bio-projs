library(shiny)
library(shinydashboard)
library(DT)
library(tidyverse)

source("scripts.R")

fun_choices <- c("promega", "qpcr_draft", "qpcr_final")

header <- dashboardHeader(title = "Mawi Data Transformation and Analysis Hub")

sidebar <- dashboardSidebar(
  
  width = 300,
  
  sidebarMenu(
    
    menuItem("Upload and transform", tabName = "Upload", icon = icon("upload"),
             
             # Input: Select a file ----
             fileInput("fileUploaded", "Choose CSV File",
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
             
             selectInput("script_fun",
                         "Choose output function",
                         choices = fun_choices),
             
             # Output: Download a file ----
             div(style="display:inline-block;width:100%;text-align: center;",
             downloadButton('downloadFile', 'Process & Download', class= "action")),
             tags$style(".skin-blue .sidebar .shiny-download-link { color: #444; }"),
             
             # CSS style for the download button ----
             div(style="display:inline-block;width:100%;text-align: center;",
             tags$style(type='text/css', "#downloadFile { width:90%; margin-top: 35px; margin-bottom: 10px}"))),
             
    menuItem("Charts", tabName = "tab1", icon = icon("chart-bar"),
             
             p("Coming soon!")
                      
                      
    )
  )
)

body <- dashboardBody(
  
  fluidRow(
    box(width = 15,
        title = "Output File Preview",
        DT::dataTableOutput("mytable"),
        style = "height:500px; overflow-x: scroll;")
  ),
  
  tags$head(tags$style(HTML('
                            .main-header .logo {
                            font-family: "Roboto", Arial, "Helvetica", sans-serif;
                            font-weight: normal;
                            font-size: 24px;
                            }
                            ')))
)

ui <- dashboardPage(header,
                    sidebar,
                    body,
                    skin = "blue")

server <- function(input, output) {
  
  # Download handler in Server
  output$downloadFile <- downloadHandler(
    filename = function() {
      paste(input$script_fun, "-", Sys.Date(), '.csv', sep='')
    },
    content = function(con) {
      originalData <- input$fileUploaded
      df <- originalData$datapath
      
      ##############################
      script_fun <- get(input$script_fun)
      final <- script_fun(df)
      
      ##############################
      dataProcessed <- final
      write.csv(dataProcessed, con)
      output$done <- renderText({
        ifelse(!is.null(dataProcessed),"Transformation done! :D","Error during the process :(")
      })
    }
  )
  
  output$mytable = DT::renderDataTable({
    get(input$script_fun)(input$fileUploaded$datapath)
  })
  
  
}
# Create Shiny app ----
shinyApp(ui, server)
