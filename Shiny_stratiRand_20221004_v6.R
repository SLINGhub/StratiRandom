
library(shiny)
library(OSAT)
library(ggplot2)

ui <- fluidPage(
  titlePanel("Stratified randomization"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      #File upload
      fileInput("file1", "Upload CSV File", multiple = TRUE, accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),

      # Display all or only the head of the table
      radioButtons("disp", "Display", choices = c(Head = "head", All = "all"), selected = "head"),
      
      tags$hr(),
      
      #Choose the covariates for the stratify randomization
      selectInput('Select1', '1st important factor', choices = NULL),
      selectInput('Select2', '2nd important factor', choices = NULL),
      selectInput('Select3', '3rd important factor', choices = NULL),
      selectInput('Select4', '4th important factor', choices = NULL),
      selectInput('Select5', '5th important factor', choices = NULL),
      selectInput('Select6', '6th important factor', choices = NULL),
      selectInput('Select7', '7th important factor', choices = NULL),
      
      #Choose number of batches
      numericInput('BatchNum', 'No. of Batches (no more than 96 samples per batch)', 3, min = 1, max = 200)
      
   ),
   
   mainPanel(
     
     img(src='SLING.png', width="40%"),
     
     tags$hr(),
     
     HTML("<b>R packages referred: "),
     tags$a(href="https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-689", "OSAT"),
     HTML(".<b>"),
     
     tags$hr(),
     
     HTML("<b>Sample size:<b>"),
     textOutput("sampleSize"),
     tags$hr(),
     
     HTML("<b>Data:<b>"),
     # Output: Raw data file
     tableOutput("contents"),
     
     tags$hr(),
     
     HTML("<br/><br/>
          <b>Please click the 'Stratify randomization' button below to stratify randomize the data and check the p-values between the batch the other covariates.  Be cautious if p-value < 0.1.</b>
          <br/><br/>"),
     
     # Click to perform stratify randomization
     inputPanel(actionButton("goButton", "Stratify randomization")),
     # Output: chi-square test of the stratify radomization data
     tableOutput("pResults"),
     
     tags$hr(),
     
     HTML("<br/><br/>
          <b>Please click the 'Download Results' button below to download the results of stratify randomization.</b>
          <br/><br/>"),
     
     # Click to download the stratify randomized data file
     inputPanel(downloadButton("downloadData","Download Data")),

     # Click to download the QC plot of the stratify randomized data file
     inputPanel(downloadButton("downloadQC","Download QC plot"))
     
   )
  )
)


server <- function(session, input, output) {
  
  # Fetch the appropriate data object
  rawData <- reactive({
    req(input$file1)
    df <- read.csv(input$file1$datapath, header=TRUE, check.names = FALSE)
    rownames(df) <- df[, 1]
    df
  })
  
  ####Filter the data, links to the selection input "Select1"
  observeEvent(input$file1, {
    updateSelectizeInput( session, 'Select1', choices = c( "NA", colnames(rawData()) )  )
  })

  observeEvent(input$Select1,{
    updateSelectInput(session,'Select2',
        choices=c( "NA", colnames(rawData())[!(colnames(rawData()) %in% input$Select1)] ) )
  })
  
  observeEvent(c(input$Select1, input$Select2),{
    updateSelectInput(session,'Select3',
        choices=c( "NA", colnames(rawData())[!(colnames(rawData()) %in% c(input$Select1, input$Select2))] ) )
  })
  
  observeEvent(c(input$Select1, input$Select2, input$Select3),{
    updateSelectInput(session,'Select4',
        choices=c( "NA", colnames(rawData())[!(colnames(rawData()) %in% c(input$Select1, input$Select2, input$Select3))] ) )
  })

  observeEvent(c(input$Select1, input$Select2, input$Select3, input$Select4),{
    updateSelectInput(session,'Select5',
                      choices=c( "NA", colnames(rawData())[!(colnames(rawData()) %in% c(input$Select1, input$Select2, input$Select3, input$Select4))] ) )
  })
  
  observeEvent(c(input$Select1, input$Select2, input$Select3, input$Select4, input$Select5),{
    updateSelectInput(session,'Select6',
                      choices=c( "NA", colnames(rawData())[!(colnames(rawData()) %in% c(input$Select1, input$Select2, input$Select3, input$Select4, input$Select5))] ) )
  })
  
  observeEvent(c(input$Select1, input$Select2, input$Select3, input$Select4, input$Select5, input$Select6),{
    updateSelectInput(session,'Select7',
                      choices=c( "NA", colnames(rawData())[!(colnames(rawData()) %in% c(input$Select1, input$Select2, input$Select3, input$Select4, input$Select5, input$Select6))] ) )
  })
  
  # Display the sample size
  output$sampleSize <- renderText({ 
    nrow(rawData())
  })
  
  # Display the selected raw data table
  output$contents <- renderTable({
    if(input$disp == "head") {return(head(rawData()))}
    else {return(rawData())}
  })
  
  # Perform stratify randomization 
  straRandCovar <- reactive({
    
    covariates <- c(input$Select1, input$Select2, input$Select3, input$Select4, input$Select5, input$Select6, input$Select7)
    covariates = covariates[!covariates %in% "NA"]

  })
  
  straRandResults <- reactive({
    gs <- setup.sample(rawData(), optimal=straRandCovar() ) ## create object that represents the used in the experiment
    gc <- setup.container(IlluminaBeadChip96Plate, input$BatchNum, batch='plates')
    ## create an optimized setup.
    # demonstration only. nSim=5000 or more are commonly used.
    gSetup <- create.optimized.setup(sample=gs, container=gc, nSim=1000)
  
    })
  
  straRandOut <- reactive({
    
    out <- get.experiment.setup(straRandResults())
    out$order <- paste(formatC(out$plate, width=3, flag="0"), "_", out$chips, "_", formatC(out$wells, width=2, flag="0"), sep = "")
    out <- out[order(out$order),]
    rownames(out) <- 1:nrow(out)
    out <- out[, -c((ncol(out)-6) : ncol(out))]
    
  })
  
  # Display the stratify randomized data table
  observeEvent(input$goButton,{
    output$pResults <- renderTable({
      return( multi.chisq.test(straRandOut(), grpVar='plates', 
            varList=straRandCovar() ) $stat ) } ) 
  })
  
  # Download the stratify randomized data file
  output$downloadData <- downloadHandler(
    filename = function() {
      fpath <- input$file1$datapath
      paste("StraRand_", input$file1, sep="")
    },
    content = function(file) {write.csv(straRandOut(), file, row.names = FALSE)}
  )

  # Download the QC plots of the stratify randomized data file
  output$downloadQC <- downloadHandler(
    filename = function() {
      fpath <- input$file1$datapath
      paste("QC_", substr(input$file1, 1, nchar(input$file1)-4), ".pdf", sep="")
    },
    content <- function(file) {
      pdf(file, width = 14, height = 5)
      
      
      ######Stratification of samples across the run order
      test_StratiRand <- get.experiment.setup(straRandResults())
      test_StratiRand$order <- paste(formatC(test_StratiRand$plate, width=3, flag="0"), "_", test_StratiRand$chips, "_", formatC(test_StratiRand$wells, width=2, flag="0"), sep = "")
      test_StratiRand <- test_StratiRand[order(test_StratiRand$order),]
      
      test_StratiRand$BatchChange <- 0
      for(i in 2:nrow(test_StratiRand))
      {
        test_StratiRand$BatchChange[i] <- as.integer(test_StratiRand$plates)[i] - as.integer(test_StratiRand$plates)[i-1]
      }
      BatchIndex <- which(test_StratiRand$BatchChange %in% 1)
      BatchLine <- BatchIndex - 0.5
      Order <- row(test_StratiRand)[,1]
      
      #colnames(test_StratiRand)
      #Convert categorical data to numberical data
      covariates <- c(input$Select1, input$Select2, input$Select3, input$Select4, input$Select5, input$Select6, input$Select7)
      covariates = covariates[!covariates %in% "NA"]
      
      cate_check <- test_StratiRand[, covariates]
      cate_check2 <- cate_check
      
      
      for( i in 1:ncol(cate_check))
      {
        levels(cate_check2[, i]) <-  1:length(levels(as.factor(cate_check2[, i])))
      }
      
      Sample_order <- data.frame(Order, cate_check2)
      
      #plot
      for( i in 2:ncol(Sample_order) )
      {
        p1 <- ggplot(Sample_order, aes(x = Order, Sample_order[, i], color = as.factor(Sample_order[, i]) )) +
          geom_segment(
            aes(
              x = Order, xend = Order,
              y = as.integer(as.factor(Sample_order[, i])) - 0.35,
              yend = as.integer(as.factor(Sample_order[, i])) + 0.35
            ),
            size = 0.5) +
          xlab("Order") +
          scale_y_discrete(
            name = "Category",
            limits = 1:length(levels(as.factor(Sample_order[, i]))),
            expand = c(0, 0.05) 
          ) +
          geom_vline(xintercept=BatchLine, colour="grey", linetype="dashed", size=0.5) + 
          theme_bw() + 
          ggtitle(colnames(Sample_order)[i]) +
          theme(legend.position="right", legend.title=element_blank()) + 
          theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank())
        print(p1)
      }
      
      
      dev.off()},
    contentType = 'image/pdf'
  )  

}


shinyApp(ui = ui, server = server)



