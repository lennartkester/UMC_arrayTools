library(shiny)

source("T:/pathologie/KMBP/EPIC/UMC_arrayTools/EPIC_functions.R")

inputChoicesEPIC <- "First select input directory"
refData <- NULL

ui <- navbarPage("UMCU Array R tools", 
                 tabPanel(tags$b("EPIC array"),              
                          # row layout with input and output definitions ----
                          fluidRow(
                            column(2,titlePanel(h4("Select folder")),actionButton("dir",tags$b("Input directory"), icon("folder-open"),style="color: #fff; background-color: #0088c7; border-color: #ffffff")),
                            column(4,titlePanel(h4("Select sample")),selectInput("samplesEPIC", NULL,choices = inputChoicesEPIC)),
                            column(2,titlePanel(h4("Plot CNVs")),actionButton("CNVplot", tags$b("Plot CNVs"), icon("paper-plane"), style="color: #fff; background-color: #0088c7; border-color: #ffffff"))
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

  observeEvent(input$dir,ignoreInit = T,{
    chosenDir <<- choose.dir(parameters$dataFolder)
    if(!is.na(chosenDir)){
      updateSelectInput(session,"samplesEPIC",choices=c("All",loadEPICsamples(chosenDir)) )
    }
    
  })
  
  observeEvent(input$CNVplot,ignoreInit = T,{
    if(is.null(refData)){
      showModal(modalDialog("Loading Controls", footer=NULL))
      refData <<- loadAndProcessControls()
      removeModal()
    }
    if(input$samplesEPIC == "All"){
      samples <- loadEPICsamples(chosenDir)
      for ( i in 1:length(samples)){
        showModal(modalDialog(paste0("Loading data and making plot for ",samples[i]), footer=NULL))
        patient <- loadAndProcessSamples(dataFolder = chosenDir,sampleName = samples[i])
        x <- segmentData(patient = patient,refData = refData,dataFolder = chosenDir)
        pdf(paste0(chosenDir,"/",sampleNames(patient),".pdf"),width=18,height = 7)
        CNV.genomeplot(x)
        dev.off()
        removeModal()
      }
    }
    if(input$samplesEPIC != "All"){
      showModal(modalDialog(paste0("Loading data and making plot for ",input$samplesEPIC), footer=NULL))
      patient <- loadAndProcessSamples(dataFolder = chosenDir,sampleName = input$samplesEPIC)
      x <- segmentData(patient = patient,refData = refData,dataFolder = chosenDir)
      pdf(paste0(chosenDir,"/",sampleNames(patient),".pdf"),width=12,height = 7)
      CNV.genomeplot(x)
      dev.off()
      output$plotEPIC <- renderPlot({
        CNV.genomeplot(x)
      })
      removeModal()
    }
  })
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)


