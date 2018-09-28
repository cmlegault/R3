#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

packages = c("shiny",       # interactive components
             "dplyr")       # data wrangling

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

# create set of random resids for demo plots using max combos of Year and Age
demo_resids_df <- data.frame(Year = rep(seq(1969, 2018), each = 20),
                             Age = rep(seq(1, 20), 50),
                             randresid = rnorm(50 * 20))

# settings for biases
bias <- list()
bias$add <- list(High = 4, Medium = 3, Low = 2)
bias$mult <- list(High = 3, Medium = 2.5, Low = 2)

#-------------------------------------------------------------------------
# plotting (and other) functions go here
# df is a data frame with at least columns Year, Age, plotresid, can have other columns
plotASAPstyle <- function(df){
  years <- sort(unique(df$Year))
  nyrs <- length(years)
  ages <- sort(unique(df$Age))
  nages <- length(ages)
  z1 <- matrix(NA, nrow = nyrs, ncol = nages)
  for (i in 1:length(df$Year)){
    z1[(df$Year[i] - min(years) + 1), df$Age[i]] <- df$plotresid[i]
  }
  scale.catch.bubble.resid <- 2 # ASAP default
  my.title <- "ASAP Age Comp Residuals"
  pos.resid.col <- "#ffffffaa"
  neg.resid.col <- "#ff1111aa"
  zr <- z1
  resid.col <- ifelse(zr > 0.0, pos.resid.col, neg.resid.col)
  plot(ages, rev(ages),  xlim = c(1, nages), ylim = c(years[nyrs],(years[1]-2)), 
       xlab = "Age", ylab = "Pearson Residuals (Obs-Pred)/SQRT(Pred*(1-Pred)/NESS)", 
       type = "n", axes=F)
  axis(1, at= ages, lab=ages)
  axis(2, at = rev(years), lab = rev(years), cex.axis=0.75, las=1)
  box()
  abline(h=years, col="lightgray")
  segments(x0=seq(ages[1], nages), y0=rep(years[1],nages),
           x1=seq(ages[1], nages), y1=rep(years[nyrs],nages), col = "lightgray", lty = 1)
  for (j in 1:nyrs){
    points(ages, rep((years[1]+j-1), nages), cex=abs(zr[j,])* scale.catch.bubble.resid,
           col="black", bg = resid.col[j,],  pch = 21 )
  }
  tmp.zr <- matrix((abs(zr)),nrow=1, ncol=length(zr[,1])*length(zr[1,])  )
  tmp.zr1 <- tmp.zr[1,]
  bubble.legend.pearson<- summary(tmp.zr1, na.rm=T) [c(2,3,5)]
  legend("topright", xpd=T, legend=round(bubble.legend.pearson,2), pch=rep(1, 3), 
         pt.cex=as.numeric(bubble.legend.pearson)*  scale.catch.bubble.resid,
         horiz=T , col='black'  )
  legend("topleft", xpd=T, legend=c("Neg.", "Pos."), pch=rep(21, 2), pt.cex=3,
         horiz=T , pt.bg=c(neg.resid.col, pos.resid.col), col="black"  )
  text(x= trunc(nages/2), y=(years[1]-1),   cex=0.8,
       label=paste("Max(resid)=",round(max(abs(zr), na.rm=T),2), sep="") )
  title (my.title, outer=T, line=-1 ) 
  title(sub=paste("Mean resid = ", round(mean(zr, na.rm=T ),2), "   SD(resid) = ", 
                  round(sd(as.vector(zr), na.rm=T ),2), sep=""), col.sub='blue', cex.sub=0.8)
  return()
}

plotBAMstyle <- function(df){
  plot(1:10, 1:10, type='n', axes = FALSE, xlab = "", ylab = "")
   text(5,5, "BAM plot under development")
  return()
}

plotr4ssstyle <- function(df){
  plot(1:10, 1:10, type='n', axes = FALSE, xlab = "", ylab = "")
   text(5,5, "r4ss plot under development")
  return()
}

plotSAMstyle <- function(df){
  # not quite the same as SAM due to location of legend bar, but close
  x <- df$Year
  y <- df$Age
  z <- df$plotresid
  SAMxlab <- "Year"
  SAMylab <- "Age"
  SAMxlim <- c(min(x) - 1, max(x) + 1)
  SAMylim <- c(min(y) - 1, max(y) + 1)
  zmax <- max(sqrt(abs(z)), na.rm = TRUE)
  SAMcex <- sqrt(abs(z))/zmax*5
  plot(x, y, xlab=SAMxlab, ylab=SAMylab, xlim=SAMxlim, ylim=SAMylim, type="n")
  neg <- z<0
  points(x[neg], y[neg], cex=SAMcex[neg], col=rgb(1, 0, 0, alpha=0.5), pch=19)
  points(x[!neg], y[!neg], cex=SAMcex[!neg], col=rgb(0, 0, 1, alpha=0.5), pch=19)
  legend("top", bty="n", legend="SAM Age Comp Residuals", text.col=gray(0.5))
  add_legend <- function(z){
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
                mar=c(0, 0, 0, 0), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    zscale <- pretty(z, min.n=4)
    uu<-par("usr")
    yy<-rep(uu[3]+.03*(uu[4]-uu[3]), length(zscale))
    xx<-seq(uu[1]+.10*(uu[2]-uu[1]),uu[1]+.4*(uu[2]-uu[1]), length=length(zscale))
    text(xx,yy,labels=zscale)
    colb <- ifelse(zscale<0, rgb(1, 0, 0, alpha=.5), rgb(0, 0, 1, alpha=.5))
    points(xx,yy,cex=sqrt(abs(zscale))/max(sqrt(abs(zscale)), na.rm=TRUE)*5, pch=19, col=colb)
  }
  add_legend(z)
  return()
}

#-------------------------------------------------------------------------
get_bias_df <- function(difficulty){
  bias_df <- list(Source = character(),
                  Value = integer())

  if (difficulty == "Easy"){
    nbiases <- 1
  } else if (difficulty == "Moderate"){
    nbiases <- sample(1:3, 1)
  } else if (difficulty == "Hard"){
    nbiases <- sample(2:5, 1)
  } else if (difficulty == "Wicked Hard"){
    nbiases <- sample(3:6, 1)
  } else {
    return("problem with difficulty switch")
  }
  Source = sample(c("Year", "Age", "Cohort"), nbiases)
  
  return(bias_df)
}

#-------------------------------------------------------------------------
# Define UI using tabs for different topics
ui <- navbarPage("Recognizing Random Residuals",
   
  tabPanel("Introduction",
    sidebarLayout(
      sidebarPanel(
        fluidRow(
          column(12,
                 h2("Welcome"),
                 br(),
                 p("This is the welcome text..."),
                 br(),
                 p("Have fun!")
                 )
        )
      ),
      mainPanel(
        h3("A title"),
        p("Some text..."),
        br(),
        h3("Another title, if needed"),
        p("Some more text...")
      )
    )
  ),
  
  tabPanel("Demo",
    sidebarLayout(
      sidebarPanel(
        selectInput("demoPlotType",
                    label = "Plot Type",
                    choices = list("ASAP", "r4ss", "BAM", "SAM", "other?"),
                    selected = "ASAP"),
        
        selectInput("demoBiasAmount",
                    label = "How much bias?",
                    choice = list("High", "Medium", "Low"),
                    selected = "High"),
        
        checkboxGroupInput("demoBiases",
                           label = "Which biases to show?",
                           choices = list("Year", "Age", "Cohort"),
                           selected = NULL),
        
        radioButtons("demoBiasType",
                     label = "Type of bias?",
                     choices = list("Additive", "Multiplicative"),
                     selected = "Additive",
                     inline = TRUE),
        
        radioButtons("demoBiasDir",
                     label = "Direction of bias?",
                     choices = list("Positive", "Negative"),
                     selected = "Positive",
                     inline = TRUE),
        
        checkboxInput("DemoBiasRescale",
                      label = "Rescale residuals?",
                      value = FALSE)
      ),
      mainPanel(
        plotOutput("demoBiasPlot")
      )
    )
  ),
  
  tabPanel("Settings",
    sidebarLayout(
      sidebarPanel(
        sliderInput("nyears",
                    label = "Number of Years",
                    min = 1,
                    max = 50,
                    step = 1,
                    value = 15),
        
        sliderInput("nages",
                    label = "Number of Ages",
                    min = 1,
                    max = 20,
                    step = 1,
                    value = 8),
        
        selectInput("plottype",
                    label = "Plot Type)",
                    choices = list("ASAP", "r4ss", "BAM", "SAM", "other?"),
                    selected = "ASAP"),
        
        selectInput("difficulty",
                    label = "How difficult?",
                    choices = list("Easy", "Moderate", "Hard", "Wicked Hard"),
                    selected = "Moderate")
        
      ),
      mainPanel(
        plotOutput("demoPlot"),
        dataTableOutput("randomTable")
      )
    )
  ),
  
  tabPanel("Random or Not?",
    sidebarLayout(
      sidebarPanel(
        actionButton("createCase",
                     label = "Create Test Case"),
        br(),
        br(),
        radioButtons("Response",
                    label = "Your response",
                    choices = list("Random", "Biased"),
                    selected = character(0),
                    inline = TRUE),
        actionButton("submit",
                     label = "Submit my response")
      ),
      mainPanel(
        textOutput("correctText"),
        plotOutput("testingPlot", height = "600px")
      )
    )
  ),
  
  tabPanel("Results so far",
    sidebarLayout(
      sidebarPanel(
        downloadButton("downloadResults", "Download")
        ),
      mainPanel(
        #plotOutput("resultsOverTimePlot"),
        #plotOutput("resultsCrossTablePlot"),
        dataTableOutput("resultsTable")
      )
    )
  )
  
) # close navbarPage parens


# Define server logic using reactive data frames
server <- function(input, output) {
  
  icount <- reactiveValues(icase = 0)
  values <- reactiveValues()
  
  # which button was last pressed?
  clickvalues <- reactiveValues(create = 0, respond = 0)
  observeEvent(input$createCase, {
    clickvalues$create <- 1
    clickvalues$respond <- 0
  })
  observeEvent(input$submit, {
    clickvalues$create <- 0
    clickvalues$respond <- 1
  })
  
  # demo values
  demoList <- reactive({
    mylist <- list()
    demonyears <- 15
    demonages <- 6
    demonvals <- demonyears * demonages
    demoendyear <- 2018
    demostartyear <- demoendyear - demonyears + 1
    demorand_df <- demo_resids_df %>%
      filter(Age <= demonages, Year >= demostartyear, Age <= demoendyear) %>%
      mutate(Cohort = Year - Age)
    mylist$resid_df <- demorand_df
    
    mylist$biased <- FALSE
    
    if (input$demoBiasType == "Additive"){
      if (input$demoBiasAmount == "High") biasamount <- bias$add$High
      if (input$demoBiasAmount == "Medium") biasamount <- bias$add$Medium
      if (input$demoBiasAmount == "Low") biasamount <- bias$add$Low
      if (input$demoBiasDir == "Negative") biasamount <- -1 * biasamount
    }
    
    if (input$demoBiasType == "Multiplicative"){
      if (input$demoBiasAmount == "High") biasamount <- bias$mult$High
      if (input$demoBiasAmount == "Medium") biasamount <- bias$mult$Medium
      if (input$demoBiasAmount == "Low") biasamount <- bias$mult$Low
      if (input$demoBiasDir == "Negative") biasamount <- 1 / biasamount
    }
    
    additivebias0 <- data.frame(Year = rep(seq(demostartyear, demoendyear), each = demonages),
                               Age = rep(seq(1, demonages), demonyears),
                               bias = 0) %>%
      mutate(Cohort = Year - Age)
    multiplicativebias1 <- data.frame(Year = rep(seq(demostartyear, demoendyear), each = demonages),
                                     Age = rep(seq(1, demonages), demonyears),
                                     bias = 1) %>%
      mutate(Cohort = Year - Age)
    additivebias <- additivebias0
    multiplicativebias <- multiplicativebias1
    
    # Checkbox group Year selected
    if (any(input$demoBiases == "Year")){
      mylist$biased <- TRUE
      if (input$demoBiasType == "Additive"){
        thisadditivebias <- additivebias0 %>%
          mutate(bias = ifelse(Year == 2014, bias + biasamount, bias))
        additivebias$bias <- additivebias$bias + thisadditivebias$bias
      }
      if (input$demoBiasType == "Multiplicative"){
        thismultiplicativebias <- multiplicativebias1 %>%
          mutate(bias = ifelse(Year == 2014, bias * biasamount, bias))
        multiplicativebias$bias <- multiplicativebias$bias * thismultiplicativebias$bias
      }
    }
    
    # Checkbox group Age selected
    if (any(input$demoBiases == "Age")){
      mylist$biased <- TRUE
      if (input$demoBiasType == "Additive"){
        thisadditivebias <- additivebias0 %>%
          mutate(bias = ifelse(Age == 3, bias + biasamount, bias))
        additivebias$bias <- additivebias$bias + thisadditivebias$bias
      }
      if (input$demoBiasType == "Multiplicative"){
        thismultiplicativebias <- multiplicativebias1 %>%
          mutate(bias = ifelse(Age == 3, bias * biasamount, bias))
        multiplicativebias$bias <- multiplicativebias$bias * thismultiplicativebias$bias
      }
    }
    
    # Checkbox group Cohort selected
    if (any(input$demoBiases == "Cohort")){
      mylist$biased <- TRUE
      if (input$demoBiasType == "Additive"){
        thisadditivebias <- additivebias0 %>%
          mutate(bias = ifelse(Cohort == 2006, bias + biasamount, bias))
        additivebias$bias <- additivebias$bias + thisadditivebias$bias
      }
      if (input$demoBiasType == "Multiplicative"){
        thismultiplicativebias <- multiplicativebias1 %>%
          mutate(bias = ifelse(Cohort == 2006, bias * biasamount, bias))
        multiplicativebias$bias <- multiplicativebias$bias * thismultiplicativebias$bias
      }
    }
    
    biasedresid_df <- mylist$randresid_df
    biasedresid_df$biasedresid <- biasedresid_df$randresid * multiplicativebias$bias + additivebias$bias
    mylist$biasedresid_df <- biasedresid_df
    
    mylist$resid_df <- mylist$resid_df %>%
      mutate(multiplicativebias = multiplicativebias$bias,
             additivebias = additivebias$bias,
             biasresid = randresid * multiplicativebias + additivebias,
             plotresid = randresid)
    if(mylist$biased == TRUE) mylist$resid_df$plotresid <- mylist$resid_df$biasresid
    
    # rescale
    if (input$DemoBiasRescale == TRUE){
      mymean <- round(mean(mylist$resid_df$plotresid), 2)
      mysd <- round(sd(mylist$resid_df$plotresid), 2)
      mylist$resid_df$plotresid <- (mylist$resid_df$plotresid - mymean) / mysd
    }
    mylist
  
  })
  
  # generate a case for evaluation
  caseList <- eventReactive(input$createCase, {
    icount$icase <- icount$icase + 1

    mylist <- list()
    mylist$icase <- icount$icase
    
    # random residuals
    nvals <- input$nyears * input$nages
    endyear <- 2018
    startyear <- endyear - input$nyears + 1
    rand_df <- expand.grid(Year = startyear:endyear, Age = 1:input$nages) %>%
      mutate(Cohort = Year - Age,
             randresid = rnorm(nvals))
    mylist$resid_df <- as.data.frame(rand_df)
    
    # determine whether to apply bias or not
    #applybias <- sample(c(TRUE, FALSE), 1)
    applybias <- FALSE
    mylist$Actual <- ifelse(applybias == TRUE, "Biased", "Random")
    
    # residuals with bias
    mylist$resid_df$biased <- mylist$resid_df$randresid
    if (applybias == TRUE){ 
      mylist$bias_df <- get_bias_df(input$difficulty)
    }

    # residuals to plot
    if (applybias == TRUE) mylist$resid_df$plotresid <- mylist$resid_df$biased
    if (applybias == FALSE) mylist$resid_df$plotresid <- mylist$resid_df$randresid

    # rescale
    applyrescale <- sample(c(TRUE, FALSE), 1)
    mylist$Rescale <- applyrescale
    if (applyrescale == TRUE){
      mymean <- mean(mylist$resid_df$plotresid)
      mysd <- sd(mylist$resid_df$plotresid)
      mylist$resid_df$plotresid <- (mylist$resid_df$plotresid - mymean) / mysd
    }
#    print(mylist)
    mylist
  })
  
  responseList <- eventReactive(input$submit,{
    if (is.null(input$Response)){
      stop("need error trap for no guess")
    }
    response.list <- list()
    response.list$icase <- icount$icase
    response.list$Response <- input$Response
    response.list$Correct <- ifelse(caseList()$Actual == input$Response, 1, 0)
    response.list
  })
  
  observeEvent(input$submit,{
    # if (is.null(input$Response)){
    #   stop("need error trap for no guess")
    # }
    thisresult_df <- data.frame(Case = icount$icase,
                                Difficult = input$difficulty,
                                Actual = caseList()$Actual,
                                Response = responseList()$Response,
                                Correct = responseList()$Correct)
    if (input$submit == 1){
      values$results_df <- thisresult_df
    }else{
      values$results_df <- rbind(values$results_df, thisresult_df)
    }
  })

  output$demoBiasPlot <- renderPlot({
    if (input$demoPlotType == "ASAP"){
      plotASAPstyle(demoList()$resid_df)
    }
    if (input$demoPlotType == "r4ss"){
      plotr4ssstyle(demoList()$resid_df)
    }
    if (input$demoPlotType == "BAM"){
      plotBAMstyle(demoList()$resid_df)
    }
    if (input$demoPlotType == "SAM"){
      plotSAMstyle(demoList()$resid_df)
    }
    if (input$demoPlotType == "other?"){
      plot(1:100,1:100, type = 'n', axes = FALSE)
      text(50,50,"need to define other") 
    }
  })
  
  output$demoPlot <- renderPlot({
    demo_df <- demo_resids_df %>%
      filter(Year >= 2018 - input$nyears + 1,
             Age <= input$nages) %>%
      mutate(plotresid = randresid)
    
    if (input$plottype == "ASAP"){
      plotASAPstyle(demo_df)
    }
    if (input$plottype == "r4ss"){
      plotr4ssstyle(demo_df)
    }
    if (input$plottype == "BAM"){
      plotBAMstyle(demo_df)
    }
    if (input$plottype == "SAM"){
      plotSAMstyle(demo_df)
    }
    if (input$plottype == "other?"){
      plot(1:100,1:100, type = 'n', axes = FALSE)
      text(50,50,"need to define other") 
    }
  })
  
  output$testingPlot <- renderPlot({
    if (clickvalues$create == 0) return(NULL)

    pr_df <- caseList()$resid_df

    if (is.null(pr_df)){
      return(NULL)
    }
    if (input$plottype == "ASAP"){
      plotASAPstyle(pr_df)
    }
    
    if (input$plottype == "r4ss"){
      plotr4ssstyle(pr_df)
    }
    
    if (input$plottype == "BAM"){
      plotBAMstyle(pr_df)
    }
    
    if (input$plottype == "SAM"){
      plotSAMstyle(pr_df)
    }
  })
   
  output$correctText <- renderText({
    if (clickvalues$respond == 0) return(NULL)
    
    ifelse (caseList()$Actual == responseList()$Response, "Correct!", "Sorry, wrong response")
  })
  
  output$randomTable <- renderDataTable(caseList()$plotresids_df)
  
  output$resultsTable <- renderDataTable(values$results_df)
   

  ## download buttons ##
  # note: if Run App in RStudio window, the filename will not default correctly (known RStudio bug),
  # but if Run App in External browser then filename will show up correctly however it will be
  # downloaded to directory C:\Users\your.name\AppData\Local\Temp
  output$downloadResults <- downloadHandler(
    filename = function() {
      paste0("R3_",Sys.Date(),".csv")
    },
    content = function(file) {
      write.csv(values$results_df, file, row.names = FALSE)
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)

