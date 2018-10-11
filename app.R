#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

packages = c("shiny",       # interactive components
             "dplyr",       # data wrangling
             "tidyr",       # data wrangling
             "ggplot2")     # nice plotting

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
bias <- data.frame(Type = c(rep("Additive", 3), rep("Multiplicative", 3)),
                   Amount = rep(c("High", "Medium", "Low"), 2),
                   Value = c(4, 2, 1, 3, 2, 1.5))

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
# Define UI using tabs for different topics
ui <- navbarPage(strong("Recognizing Random Residuals"), 
                 header = img(src="noaanefsclogo.png"), 
                 windowTitle = "R3",
                 
  tabPanel("Introduction",
    sidebarLayout(
      sidebarPanel(
        fluidRow(
          column(12,
                 h2("Welcome"),
                 br(),
                 p("Think you can recognize when residuals are random versus biased? Here's a game to let you see how good you really are. Start with the 'Demo' tab to see how the different biases look. Use the 'Settings' tab to create a situation like you are accustomed to seeing in your or someone else's assessment. The 'Random or Not?' tab is the fun part. A plot is provided and you guess whether the residuals are random or biased with immediate feedback. The 'Results so far' tab lets you see how you are doing. Can you do better than random? If so, bump up the difficulty a notch and try again."),
                 br(),
                 p("Have fun!")
                 )
        )
      ),
      mainPanel(
        h3("Demo"),
        p("In this tab, a set of residuals is created randomly from a normal distribution with mean zero and standard deviation one. This is the desired property of residuals in age composition plots. There are currently two types of plots available, ASAP and SAM. You can then add bias as a Year, Age, or Cohort effect. The 2006 Cohort is age 0 in 2006, age 1 in 2007, age 2 in 2008, etc. The higher the bias, the easier is it to see the effect of adding or multiplying the random values by the bias. The bias can be either positive or negative. In this demonstration, mixing of additive and multiplicative with positive and negative biases is not allowed, but it can occur in the test plots."),
        h3("Settings"),
        p("This is where you determine how the test plots will appear and how challenging it will be to detect the biases. The sliders allow you to change the number of years and ages to mimic situations that you encounter in your assessment. There are only two plot types available currently, ASAP and SAM, because the r4ss and BAM plots require the observed and expected values, not just a matrix of residuals. This will be developed in the future. Feedback is welcome on the difficulty settings. Is Easy too easy, or Hard too hard? The plot on this page shows an example of the settings selected with no bias."),
        h3("Random or Not?"),
        p("This is the fun part! Click the Create Test Case button to generate a test plot. Select either Random or Biased and click the Submit my response button. You will see the plot disappear and be replace by a message of either 'Correct!' or Sorry, wrong response'. If the residuals were biased, a small table will also appear showing all the biases that were present in the residuals (did you see them all?). A smaller version of the test plot is shown at the bottom of the page (you may have to scroll down to see it) so you can see where biases were or were not present. Just click the Create Test Case again to get your next test plot."),
        p(" Please note there are some known issues with this page due to my relative newness with Shiny. Specifically, if you click either button multiple times, it will advance the counter and create problems for the percent correct calculations. If you are feeling bad about how you are doing, you can simply click the submit my response button numerous times on a correct reponse to inflate your percent correct. <grin>"),
        h3("Results so far"),
        p("This is where you can see your progress. If you've responded fewer than five times, the top plot will just show a message that you haven't done enough yet. Once you've completed at least five test plot responses, the plot shows your percent correct as the solid line with the red area indicating the 95% confidence interval associated with random guessing. So if you're in the red area, you're only doing as well as flipping a coin! If you are below the red area, then you should spend some more time with the Demo tab to see how residuals respond to biases. If you are above the red area, then you are better than random and can say that yes, you can recognize random residuals! The red area appears jagged at small number of responses due to the confidence interval resulting in whole numbers. As the number of responses gets large, the red area appears much smoother. The first table compiles your correct and incorrect responses according to the level of difficulty and whether the test plot was random or biased. The next table contains the response-specific results and can be sorted by any of the columns by clicking on the up or down arrow to the right of the column header."),
        h3("Technical Details"),
        p("The additive biases just add or subtract the Value in the table below to the random residuals in that year, age, or cohort. The multiplicative biases multiply the residuals by the Value when the direction is positive and multiply the residuals by (1/Value) when the direction is negative (meaning the size of the residuals decrease)."),
        tableOutput("biasTable"),
        p("Cohort biases can be applied at age zero in years startyear - 3 to endyear - 3. Age and Year biases can be applied in any of the years or ages. Any given year, age, or cohort cannot have multiple biases applied, meaning there cannot be a case when year 2014 has both an additive and multiplicative bias. However, intersections of a year, age, or cohort can have multiple biases applied as seen in the 'Demo' tab."),
        p("The difficulty settings change the number of biases possible, the amount of bias, and the type of bias. When the difficulty setting is Easy, there is only one additive bias applied and it has a high amount of bias. The Moderate setting has either 3 or 4 additive biases applied with either medium or low amounts of bias. The Hard setting has from 1 to 3 biases applied that can be either additive or multiplicative with low amount of bias. The Easy, Moderate, and Hard settings have the direction of each bias applied randomly (meaning both positive and negative or all positive or all negative directions can occur). The Wicked Hard setting applies two multiplicative biases, one positive and one negative, with low amounts of bias."),
        p("Code for this R Shiny app is available at ",
        a("https://github.com/cmlegault/R3", href="https://github.com/cmlegault/R3", target="_blank"))
      )
    )
  ),
  
  tabPanel("Demo",
    sidebarLayout(
      sidebarPanel(
        selectInput("demoPlotType",
                    label = "Plot Type",
                    choices = list("ASAP", "r4ss", "BAM", "SAM"),
                    selected = "ASAP"),
        
        selectInput("demoBiasAmount",
                    label = "How much bias?",
                    choice = list("High", "Medium", "Low"),
                    selected = "High"),
        
        checkboxGroupInput("demoBiases",
                           label = "Which biases to show?",
                           choices = list("Year 2014", "Age 3", "Cohort 2006"),
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
                     inline = TRUE)
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
                    min = 6,
                    max = 50,
                    step = 1,
                    value = 15),
        
        sliderInput("nages",
                    label = "Number of Ages",
                    min = 4,
                    max = 20,
                    step = 1,
                    value = 8),
        
        selectInput("plottype",
                    label = "Plot Type)",
                    choices = list("ASAP", "SAM"),
                    selected = "ASAP"),
        
        selectInput("difficulty",
                    label = "How difficult?",
                    choices = list("Easy", "Moderate", "Hard", "Wicked Hard"),
                    selected = "Moderate")
        
      ),
      mainPanel(
        plotOutput("demoPlot")
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
                    selected = "Random",
                    inline = TRUE),
        actionButton("submit",
                     label = "Submit my response")
      ),
      mainPanel(
        textOutput("correctText"),
        br(),
        br(),
        tableOutput("biasesTable"),
        plotOutput("testingPlot", height = "600px"),
        plotOutput("testingPlotsmall", height = "400px")
      )
    )
  ),
  
  tabPanel("Results so far",
           plotOutput("resultsOverTimePlot"),
           tableOutput("crossTable"),
           dataTableOutput("resultsTable")
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
    
    biasamount <- bias %>%
      filter(Type == input$demoBiasType, Amount == input$demoBiasAmount) %>%
      select(Value) %>%
      as.numeric()
    
    if (input$demoBiasType == "Additive"){
      if (input$demoBiasDir == "Negative") biasamount <- -1 * biasamount
    }
    
    if (input$demoBiasType == "Multiplicative"){
      if (input$demoBiasDir == "Negative") biasamount <- 1 / biasamount
    }
    
    bias0 <- expand.grid(Age = 1:demonages, Year = demostartyear:demoendyear) %>%
      mutate(Cohort = Year - Age) 
    
    additivebias0 <- bias0 %>%
      mutate(bias = 0)
    multiplicativebias1 <- bias0 %>%
      mutate(bias = 1)

    additivebias <- additivebias0
    multiplicativebias <- multiplicativebias1
    
    # Checkbox group Year selected
    if (any(input$demoBiases == "Year 2014")){
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
    if (any(input$demoBiases == "Age 3")){
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
    if (any(input$demoBiases == "Cohort 2006")){
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
    
    mylist
  
  })
  
  #################################################################################
  #################################################################################
  #################################################################################
  # generate a case for evaluation
  caseList <- eventReactive(input$createCase, {
    icount$icase <- icount$icase + 1

    mylist <- list()
    mylist$icase <- icount$icase
    
    # random residuals
    nvals <- input$nyears * input$nages
    endyear <- 2018
    startyear <- endyear - input$nyears + 1
    rand_df <- expand.grid(Age = 1:input$nages, Year = startyear:endyear) %>%
      mutate(Cohort = Year - Age,
             randresid = rnorm(nvals))
    mylist$resid_df <- as.data.frame(rand_df)
    
    # determine whether to apply bias or not
    applybias <- sample(c(TRUE, FALSE), 1)
    mylist$Actual <- ifelse(applybias == TRUE, "Biased", "Random")
    
    # residuals with bias
    mylist$resid_df$biased <- mylist$resid_df$randresid
    if (applybias == TRUE){
      dirflag <- TRUE
      if (input$difficulty == "Easy"){
        nbiases <- 1
        biasType <- "Additive"
        biasAmount <- "High"
      } else if (input$difficulty == "Moderate"){
        nbiases <- sample(3:4, 1)
        biasType <- rep("Additive", nbiases)
        biasAmount <- sample(c("Medium", "Low"), nbiases, replace = TRUE)
      } else if (input$difficulty == "Hard"){
        nbiases <- sample(1:3, 1)
        biasType <- sample(c("Additive", "Multiplicative"), nbiases, replace = TRUE)
        biasAmount <- rep("Low", nbiases)
      } else if (input$difficulty == "Wicked Hard"){
        nbiases <- 2
        biasType <- rep("Multiplicative", nbiases)
        biasAmount <- rep("Low", nbiases)
        dirflag <- FALSE
      } else {
        return("problem with difficulty switch")
      }
      Source <- rep(NA, nbiases)
      SourceValue <- rep(NA, nbiases)
      biasvalue <- rep(NA, nbiases)
      nbiasyear <- 0
      nbiasage <- 0
      nbiascohort <- 0
      for (ibias in 1:nbiases){
        Source[ibias] <- sample(c("Year", "Age", "Cohort"), 1)
        if (Source[ibias] == "Year") nbiasyear <- nbiasyear + 1
        if (Source[ibias] == "Age") nbiasage <- nbiasage + 1
        if (Source[ibias] == "Cohort") nbiascohort <- nbiascohort + 1
      }
      # which years, ages, or cohororts are biased
      biasyears <- sample(seq(startyear, endyear), nbiasyear, replace = FALSE)
      biasages <- sample(seq(1, input$nages), nbiasage, replace = FALSE)
      biascohorts <- sample(seq(startyear - 3, endyear - 3), nbiascohort, replace = FALSE)

      biasdirection <- sample(c("Positive", "Negative"), nbiases, replace = dirflag)
      
      icounty <- 0
      icounta <- 0
      icountc <- 0
      for (ibias in 1:nbiases){
        biasvalue[ibias] <- bias %>%
          filter(Type == biasType[ibias], Amount == biasAmount[ibias]) %>%
          select(Value) %>%
          as.numeric()
        
        if (biasType[ibias] == "Additive"){
          if (biasdirection[ibias] == "Negative") biasvalue[ibias] <- -1 * biasvalue[ibias]
        }
        if (biasType[ibias] == "Multiplicative"){
          if (biasdirection[ibias] == "Negative") biasvalue[ibias] <- 1 / biasvalue[ibias]
        }
        if (Source[ibias] == "Year"){
          icounty <- icounty + 1
          SourceValue[ibias] <- biasyears[icounty]
        }
        if (Source[ibias] == "Age"){
          icounta <- icounta + 1
          SourceValue[ibias] <- biasages[icounta]
        }
        if (Source[ibias] == "Cohort"){
          icountc <- icountc + 1
          SourceValue[ibias] <- biascohorts[icountc]
        }
      }
      bias_df <- data.frame(Count = 1:nbiases,
                            Source = Source,
                            SourceValue = SourceValue,
                            Type = biasType,
                            Amount = biasAmount,
                            Direction = biasdirection,
                            BiasValue = biasvalue)
      
      mylist$bias_df <- bias_df
      
      # bias containers
      bias0 <- expand.grid(Age = 1:input$nages, Year = startyear:endyear) %>%
        mutate(Cohort = Year - Age) 
      
      additivebias0 <- bias0 %>%
        mutate(bias = 0)
      multiplicativebias1 <- bias0 %>%
        mutate(bias = 1)
      
      additivebias <- additivebias0
      multiplicativebias <- multiplicativebias1
      
      for (ibias in 1:nbiases){
        # Year bias
        if (bias_df$Source[ibias] == "Year"){
          if (bias_df$Type[ibias] == "Additive"){
            thisadditivebias <- additivebias0 %>%
              mutate(bias = ifelse(Year == bias_df$SourceValue[ibias], 
                                   bias + bias_df$BiasValue[ibias], bias))
            additivebias$bias <- additivebias$bias + thisadditivebias$bias
          }
          if (bias_df$Type[ibias] == "Multiplicative"){
            thismultiplicativebias <- multiplicativebias1 %>%
              mutate(bias = ifelse(Year == bias_df$SourceValue[ibias], 
                                   bias * bias_df$BiasValue[ibias], bias))
            multiplicativebias$bias <- multiplicativebias$bias * thismultiplicativebias$bias
          }
        }
        
        # Age bias
        if (bias_df$Source[ibias] == "Age"){
          if (bias_df$Type[ibias] == "Additive"){
            thisadditivebias <- additivebias0 %>%
              mutate(bias = ifelse(Age == bias_df$SourceValue[ibias], 
                                   bias + bias_df$BiasValue[ibias], bias))
            additivebias$bias <- additivebias$bias + thisadditivebias$bias
          }
          if (bias_df$Type[ibias] == "Multiplicative"){
            thismultiplicativebias <- multiplicativebias1 %>%
              mutate(bias = ifelse(Age == bias_df$SourceValue[ibias], 
                                   bias * bias_df$BiasValue[ibias], bias))
            multiplicativebias$bias <- multiplicativebias$bias * thismultiplicativebias$bias
          }
        }
        
        # Cohort bias
        if (bias_df$Source[ibias] == "Cohort"){
          if (bias_df$Type[ibias] == "Additive"){
            thisadditivebias <- additivebias0 %>%
              mutate(bias = ifelse(Cohort == bias_df$SourceValue[ibias], 
                                   bias + bias_df$BiasValue[ibias], bias))
            additivebias$bias <- additivebias$bias + thisadditivebias$bias
          }
          if (bias_df$Type[ibias] == "Multiplicative"){
            thismultiplicativebias <- multiplicativebias1 %>%
              mutate(bias = ifelse(Cohort == bias_df$SourceValue[ibias], 
                                   bias * bias_df$BiasValue[ibias], bias))
            multiplicativebias$bias <- multiplicativebias$bias * thismultiplicativebias$bias
          }
        }
      }      
      
      # apply the bias
      mylist$resid_df$biased <- mylist$resid_df$randresid * multiplicativebias$bias + additivebias$bias
    }

    # residuals to plot
    if (applybias == TRUE) mylist$resid_df$plotresid <- mylist$resid_df$biased
    if (applybias == FALSE) mylist$resid_df$plotresid <- mylist$resid_df$randresid

    mylist
  })
  
  responseList <- eventReactive(input$submit,{
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
                                Difficulty = input$difficulty,
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
   
  output$testingPlotsmall <- renderPlot({
    if (clickvalues$respond == 0) return(NULL)
    
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
  
  output$biasesTable <- renderTable({
    if (is.null(caseList()$bias_df) | clickvalues$respond == 0){
      return(NULL)
    }
    caseList()$bias_df
  })
  
  output$biasTable <- renderTable({
    bias
  })
  
  output$resultsOverTimePlot <- renderPlot({
    if (is.null(values$results_df)){
      return(NULL)
    }
    results_plot_df <- values$results_df %>%
      mutate(cdf = cumsum(Correct)) %>%
      mutate(percCorrect = 100 * cdf / Case) %>%
      mutate(lowerCI = 100 * qbinom(0.025, Case, 0.5) / Case) %>%
      mutate(upperCI = 100 * qbinom(0.975, Case, 0.5) / Case)
    d_ends <- results_plot_df %>%
      filter(Case == max(Case))
    if (length(results_plot_df$Case) < 5){
      plot(1:10,1:10,type='n',axes=FALSE,xlab="",ylab="")
        text(5,5,"Need at least 5 responses to plot time series")
    }else{
      ggplot(results_plot_df, aes(x=Case, y=percCorrect)) +
        geom_point() +
        geom_line() +
        geom_ribbon(aes(ymin=lowerCI, ymax=upperCI), fill = "red", alpha = 0.5) +
        xlab("Number of Responses") +
        ylab("Percent Correct") +
        expand_limits(y=c(0, 100)) +
        scale_y_continuous(sec.axis = sec_axis(~ ., breaks = round(d_ends$percCorrect, 1))) +
        scale_x_continuous(expand = c(0,0)) +
        theme_bw()
    }
  })
  
  output$crossTable <- renderTable({
    values$results_df %>%
      group_by(Difficulty, Actual, Response) %>%
      summarize(n = n()) %>%
      spread(key = Response, value = n, fill = 0) %>%
      mutate(PercentCorrect = ifelse(Actual == "Biased", 
                                     100 * Biased / (Biased + Random), 
                                     100 * Random / (Biased + Random))) 
  })
  
  output$resultsTable <- renderDataTable(values$results_df)
   
}

# Run the application 
shinyApp(ui = ui, server = server)

