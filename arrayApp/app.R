library(shiny)

source("/Users/lennartkester/Documents/EPIC/UMC_arrayTools/EPIC_functions.R")

inputChoicesEPIC <- loadEPICFolders()

ui <- navbarPage("UMCU Array R tools", 
                 tabPanel(tags$b("EPIC array"),              
                          # row layout with input and output definitions ----
                          fluidRow(
                            column(4,titlePanel(h4("Make EPIC CNV plot")),selectInput("sampleFolderEPIC", "Choose a sample folder:",choices = inputChoicesEPIC)),
                           ),
                          fluidRow(
                            column(2,actionButton("CNVplot", tags$b("Plot CNVs"), icon("paper-plane"), style="color: #fff; background-color: #0088c7; border-color: #ffffff")),
                            column(2,actionButton("refreshFolders", "Refresh")),
                            column(2,checkboxInput("EPICpdf","make PDF")),
                            
                           ),
                          fluidRow(
                            br(),
                            br(),
                            mainPanel(width = 12,plotOutput("plotEPIC"))
                          )
                 ),
                 tags$style(HTML(".navbar-default .navbar-brand {color: #ffffff;}
                                 .navbar { background-color: #0088c7;}
                                 .navbar-default .navbar-nav > li > a {color:#ffffff;}
                                 "))
)


server <- function(input, output, session) {
  observeEvent(input$CNVplot,ignoreInit = T,{
    showModal(modalDialog("Loading data and making plot", footer=NULL))
    if(!is.null(refData)){
      refData <- loadAndProcessControls()    
    }
    patient <- loadAndProcessSamples(input$sampleFolderEPIC)
    x <- segmentData(patient,refData)
    if(input$EPICpdf){
      pdf(paste0(parameters$dataFolder,sampleNames(patient),"/",sampleNames(patient),".pdf"),width=12,height = 7)
      CNV.genomeplot(x)
      dev.off()
    }
    removeModal()
    output$plotEPIC <- renderPlot({
      CNV.genomeplot(x)
    })
  })
  
  observeEvent(input$refreshFolders,ignoreInit = T,{
    updateSelectInput(session,"sampleFolderEPIC",choices=loadEPICFolders() )
  })

  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)


