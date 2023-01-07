library(shiny)
library(ellipse)


source('mvmixdist.R')


# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # css hack
  tags$style(HTML("
    thead:first-child > tr:first-child > th {
        border-top: 0;
        font-weight: normal;
    }
")),
  
  # Application title
  titlePanel("Bayesian fitting of two peaks and a constant background"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      HTML("<center>"),
      "Proportions",
      tableOutput("props"),
      "Uniform background (bg)",
      tableOutput("box"),
      "Multivariate normal distribution #1 (mvn1)",
      tags$table(align="center",
                 tags$td(align="center",
                         HTML("&mu;<sub>1</sub>"),
                         tableOutput("mu1")
                 ),
                 tags$td(align="center",
                         HTML("&Sigma;<sub>1</sub>"),
                         tableOutput("sigma1")
                 )
      ),
      "Multivariate normal distribution #2 (mvn2)",
      tags$table(align="center",
                 tags$td(align="center",
                         HTML("&mu;<sub>2</sub>"),
                         tableOutput("mu2")
                 ),
                 tags$td(align="center",
                         HTML("&Sigma;<sub>2</sub>"),
                         tableOutput("sigma2")
                 )
      ),
      sliderInput("numSamples", "Number samples",
                  min = 1, max = 500, value = 200, step=1),
      actionButton("drawSample", "Draw sample"),
      actionButton("getPosterior", "Get Posterior"),
      actionButton("getMLEstimate", "Get ML estimate"),
      HTML("</center>")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Introduction",
                 tagList(
                   tags$br(),
                   tags$p("In experimental nuclear physics during data analysis, one is 
                        facing the task to split the measured signal of a particle detector
                        in several contributing components in order to extract the relevant one.
                        The relevant contribution are the particles being produced in a
                        controlled way within the experiment whereas the
                        other contributions are associated with particles entering the
                        detector from the outside (e.g., from space) or particles that
                        are produced in the experiment but whose impact on the detector
                        output is undesired.
                        "),
                   tags$p("In statistics processes being made up of several independent subprocesses",
                          "can be modeled as ", tags$a("mixture models", 
                                                     href="https://en.wikipedia.org/wiki/Mixture_model"), ".",
                          "In the example presented here, we assume that the detector is capable ",
                          "to resolve the position where a particle has hit the detector, i.e. ",
                          "the position on a two dimensional plane. ",
                          "Further, we assume that there are two subprocesses that lead to distribution patterns ",
                          "which can be described by bivariate normal distributions and one process that ",
                          "produces a uniform distribution of particles on the detection surface. ",
                          "The parameters controlling the specifics of this scenario can be adjusted using the ",
                          "control elements on the left. The abbreviation ", tags$em("bg"), " denotes the background",
                          "contribution and ", tags$em("mvn1"), " and ", tags$em("mvn2"), " the first and second ",
                          "bivariate normal distribution, respectively.",
                          "A brief description of the various sections in the control panel follows:",
                          tags$dl(tags$dt("Proportions"),
                                  tags$dd("Determines the proportion of particles associated with the different",
                                          "subprocesses. For instance, 20% of particles due to background, ",
                                          "30% due to the first and 50% due to the second bivariate normal ",
                                          "distribution corresponds to the specifications 0.2, 0.3, and 0.5. ",
                                          "If the proportions do not sum to one, they will be normalized."),
                                  tags$dt("Uniform background"),
                                  tags$dd("The fields in this section specify the area in which particles ",
                                          "associated with the background are produced. "),
                                  tags$dt("Multivariate normal distribution #1"),
                                  tags$dd("The two fields below ", HTML("&mu;<sub>1</sub>"),
                                          "determine the x- and y-position of the center of ", tags$em("mvn1"), ". ",
                                          "The top-left and bottom-right field below ", 
                                          HTML("&Sigma;<sub>1</sub>"),
                                          " determine the horizontal and vertical extension, respectively, ",
                                          "of the distribution. The top-right and bottom-left fields define ",
                                          "the covariance between the x- and y-position, and should be the same.",
                                          "For the sake of illustration, they can be left at zero."),
                                  tags$dt("Multivariate normal distribution #2"),
                                  tags$dd("The fields controlling the characteristics of the second bivariate ",
                                          "normal distribution ", tags$em("mvn2"), " are analogously defined to ",
                                          "those for ", tags$em("mvn1")))),
                   tags$p("The above numbers completely define the data generation process.",
                          "A click on the button 'Draw sample' produces points according to the specification, ",
                          "which can be seen by clicking on the tab 'Scatter plot' above. ",
                          "The parameters in the control panel can be modified at any time and a click ",
                          "on 'Draw sample' updates the scatter plot."),
                   tags$p("These data points can now be used to estimate the underlying parameters of the ",
                          "data generation process, i.e. the values in the fields of the control panel. ",
                          "Two inference approaches are implemented in this demonstration. ",
                          "Clicking the button 'Get ML estimate' will try to recover the values for the ",
                          "parameters using the ",
                          tags$a("expectation-maximization algorithm",
                                 href="https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm"), ". ",
                          "A click on the button 'Get posterior' yields the estimate using ",
                          tags$a("Gibbs sampling", href="https://en.wikipedia.org/wiki/Gibbs_sampling"),
                          "(prior on proportions is a symmetric Dirichlet pdf and an inverse-Wishart pdf for each mvn component)",
                          "The results are shown in both cases in the 'Scatter plot' tab.",
                          "If Gibbs sampling has been used, the trace plots are available under the equally named tab.")
                 )
        ),
        tabPanel("Scatter plot", 
                 plotOutput("distPlot"),
                 tableOutput("postProps"),
                 tableOutput("postMean1"),
                 tableOutput("postMean2")
        ),
        tabPanel("Trace plots",
                 plotOutput("tracePlotMean1"),
                 plotOutput("tracePlotMean2"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  initDistPars <- list(
    props = rbind(c(0.1, 0.4, 0.5)),
    box = t(matrix(c(-20, -20, 40, 40), 2, 2)),
    mu1 = c(3,2),
    mu2 = c(15,15),
    sigma1 = matrix(c(2,0,0,4), nrow=2, ncol=2),
    sigma2 = matrix(c(5,0,0,5), nrow=2, ncol=2),
    numSamples = 200
  )
  colnames(initDistPars$props) <-  c("bg", "mvn1", "mvn2")
  
  rownames(initDistPars$box) <- c("lower-left",
                                  "upper-right")
  colnames(initDistPars$box) <- c("x", "y")  
  
  with(initDistPars, {
    output$props <- drawNumTable(
      "props", rbind(props),
      cn = TRUE, min=0, max=1, step=0.05) 
    output$box <- drawNumTable(
      "box", rbind(box),
      rn=TRUE, cn=TRUE, min=-100, max=100, step=5)
    output$mu1 <- drawNumTable("mu1", cbind(mu1),
                               min=-100, max=100, step=5)
    output$mu1 <- drawNumTable("mu1", cbind(mu1),
                               min=-100, max=100, step=5)
    output$sigma1 <- drawNumTable("sigma1", sigma1,
                                  min=-100, max=100, step=5)
    output$mu2 <- drawNumTable("mu2", cbind(mu2),
                               min=-100, max=100, step=5)   
    output$sigma2 <- drawNumTable("sigma2", sigma2,
                                  min=-100, max=100, step=5)
  })
  
  postStruc <- reactiveVal(value = NULL)
  
  distPars <- eventReactive(input$drawSample, {
    
    newDistPars <- list(
      props = sapply(1:3, function(i) input[[paste0("props",i)]]),
      box = sapply(c(1,3,2,4), function(i) input[[paste0("box",i)]]),
      mu1 = sapply(1:2, function(i) input[[paste0("mu1",i)]]),
      mu2 = sapply(1:2, function(i) input[[paste0("mu2",i)]]),
      sigma1 = matrix(
        sapply(1:4,  function(i) input[[paste0("sigma1",i)]]),
        nrow = 2, ncol = 2),
      sigma2 = matrix(
        sapply(1:4, function(i) input[[paste0("sigma2",i)]]),
        nrow = 2, ncol = 2),
      numSamples = input[["numSamples"]])
    newDistPars$props <- pmin(pmax(newDistPars$props,0),1)
    newDistPars$props <- with(newDistPars, props / sum(props))
    sapply(1:3, function(i)
      updateNumericInput(session, paste0("props",i), 
                         value = newDistPars$props[i]))
    # also append the distribution objects
    mvnPrior <- list(mu0 = c(5,5), kappa0 = 1,
                    nu0 = 2, phi = rbind(c(1,0), c(0,1)))
    mixPrior  <- list(alpha = 1)
    
    constDist <- with(newDistPars, createDist_Const(box))
    mvnDist1 <- with(newDistPars, createDist_MVN(mu1, sigma1,
                                                 prior = mvnPrior))
    mvnDist2 <- with(newDistPars, createDist_MVN(mu2, sigma2,
                                                 prior = mvnPrior))
    mixDist <- isolate(with(
      newDistPars, 
      createDist_Mix(prop=props,
                     list(constDist, mvnDist1, mvnDist2),
                     prior = mixPrior)))
    obsPoints <- getSample(mixDist, newDistPars$numSamples)
    newDistPars$obsPoints <- obsPoints
    newDistPars$mixDist <- mixDist
    # return everything
    postStruc(NULL)
    newDistPars
  })
  
  observeEvent(input$getPosterior, {
    
    x <- postStruc()
    propsObs <- x$propsObs
    mean1Obs <- x$mean1Obs
    mean2Obs <- x$mean2Obs
    sigma1Obs <- x$sigma1Obs
    sigma2Obs <- x$sigma2Obs
    
    x <- with(distPars(), getPosteriorSample(mixDist, obsPoints, 200))
    propsObs <- cbind(propsObs, sapply(x, function(x) x@prop))
    mean1Obs <- cbind(mean1Obs, sapply(x, function(x) x@comp[[2]]@mean))
    mean2Obs <- cbind(mean2Obs, sapply(x, function(x) x@comp[[3]]@mean))
    sigma1Obs <- cbind(sigma1Obs, sapply(x, function(x) as.vector(x@comp[[2]]@sigma)))
    sigma2Obs <- cbind(sigma2Obs, sapply(x, function(x) as.vector(x@comp[[3]]@sigma)))
    
    postStruc(list(
      numObs = ncol(propsObs),
      propsObs = propsObs,
      propsEst = rowMeans(propsObs),
      propsSd = apply(propsObs, 1, sd),
      mean1Obs = mean1Obs,
      mean1Est = rowMeans(mean1Obs),
      mean1Sd = apply(mean1Obs, 1, sd),
      sigma1Obs = sigma1Obs,
      sigma1Est = rowMeans(sigma1Obs),
      sigma1Sd = apply(sigma1Obs, 1, sd),
      mean2Obs = mean2Obs,
      mean2Est = rowMeans(mean2Obs),
      mean2Sd = apply(mean2Obs, 1, sd),
      sigma2Obs = sigma2Obs,
      sigma2Est = rowMeans(sigma2Obs),
      sigma2Sd = apply(sigma2Obs, 1, sd)
    ))
  })
  
  
  observeEvent(input$getMLEstimate, {
    
    x <- with(distPars(), list(getMLEstimate(mixDist, obsPoints)))
    propsObs <- sapply(x, function(x) x@prop)
    mean1Obs <- sapply(x, function(x) x@comp[[2]]@mean)
    mean2Obs <- sapply(x, function(x) x@comp[[3]]@mean)
    sigma1Obs <- sapply(x, function(x) as.vector(x@comp[[2]]@sigma))
    sigma2Obs <- sapply(x, function(x) as.vector(x@comp[[3]]@sigma))
    
    postStruc(list(
      propsObs = propsObs,
      propsEst = rowMeans(propsObs),
      propsSd = NA,
      mean1Obs = mean1Obs,
      mean1Est = rowMeans(mean1Obs),
      mean1Sd = NA,
      sigma1Obs = sigma1Obs,
      sigma1Est = rowMeans(sigma1Obs),
      sigma1Sd = NA,
      mean2Obs = mean2Obs,
      mean2Est = rowMeans(mean2Obs),
      mean2Sd = NA,
      sigma2Obs = sigma2Obs,
      sigma2Est = rowMeans(sigma2Obs),
      sigma2Sd = NA
    ))
  })
  
  numUncToStr <- function(est, unc) {
    paste0(sprintf("%.2f", est),
           ifelse(!is.na(unc), paste0(" ± ", sprintf("%.2f", unc)),  ""))
  }
  

  output$distPlot <- renderPlot({
    
    with(distPars(), {
         plot(t(obsPoints), xlim=box[c(1,3)], 
              ylim=box[c(2,4)], cex=0.5,
              xlab = "x", ylab="y")
         for (i in 1:2) {
           curMean <- distPars()[[paste0("mu", i)]]
           curCov <- distPars()[[paste0("sigma", i)]]
           dim(curCov) <- c(2,2)
           for (l in c(0.68,0.95)) {
             myEll <- ellipse(x = curCov, centre = curMean, level = l)
             lines(myEll, col="green")
           }
         }
    })
    if (!is.null(postStruc())) {
      
      for (i in 1:2) {
        curMean <- postStruc()[[paste0("mean", i, "Est")]]
        curCov <- postStruc()[[paste0("sigma", i, "Est")]]
        dim(curCov) <- c(2,2)
        for (l in c(0.68,0.95)) {
          myEll <- ellipse(x = curCov, centre = curMean, level = l)
          lines(myEll, col="red")
        }
      }
    }
    
  })
  
  output$tracePlotMean1 <- renderPlot({
    
    if (!is.null(postStruc()) && isTRUE(postStruc()$numObs > 1)) 
      with(postStruc(), {
        yrng <- range(c(mean1Obs[1,],mean1Obs[2,]))
        plot(seq_len(ncol(mean1Obs)), mean1Obs[1,], type="l", ylim=yrng, col="green",
             main = expression("Monte Carlo chain of" ~ mu[1*x] ~ "(green) and" ~ mu[1*y] ~ "(red)"),
             xlab = "iteration number", ylab = expression("mean of" ~ mu[1]))
        lines(seq_len(ncol(mean1Obs)), mean1Obs[2,], col="red")
      })
  })
  
  output$tracePlotMean2 <- renderPlot({
    
    if (!is.null(postStruc()) && isTRUE(postStruc()$numObs > 1)) 
      with(postStruc(), {
        yrng <- range(c(mean2Obs[1,],mean2Obs[2,]))
        plot(seq_len(ncol(mean2Obs)), mean2Obs[1,], type="l", ylim=yrng, col="green",
             main = expression("Monte Carlo chain of" ~ mu[2*x] ~ "(green) and" ~ mu[2*y] ~ "(red)"),
             xlab = "iteration number", ylab = expression("mean of" ~ mu[2]))
        lines(seq_len(ncol(mean2Obs)), mean2Obs[2,], col="red")
      })
  })
  
  
  output$postProps <- renderTable({
    
    if (!is.null(postStruc())) {
      with(postStruc(),
           data.frame(
             "Est props" = numUncToStr(propsEst[1], propsSd[1]),
             " " = numUncToStr(propsEst[2], propsSd[2]),
             " " = numUncToStr(propsEst[3], propsSd[3]),
             check.names=FALSE)
      )}
  }, sanitize.text.function = identity)
  
  output$postMean1 <- renderTable({
    
    if (!is.null(postStruc())) {
      with(postStruc(),
           data.frame(
             "Est &mu;<sub>1</sub>" = numUncToStr(mean1Est, mean1Sd),
             "Est &Sigma;<sub>1</sub>" = numUncToStr(sigma1Est[1:2], sigma1Sd[1:2]),
             " " = numUncToStr(sigma1Est[3:4], sigma1Sd[3:4]),
             check.names=FALSE)
      )}
  }, sanitize.text.function = identity)
  
  output$postMean2 <- renderTable({
    
    if (!is.null(postStruc())) {
      with(postStruc(),
           data.frame(
             "Est &mu;<sub>2</sub>" = numUncToStr(mean2Est, mean2Sd),
             "Est &Sigma;<sub>2</sub>" = numUncToStr(sigma2Est[1:2], sigma2Sd[1:2]),
             " " = numUncToStr(sigma2Est[3:4], sigma2Sd[3:4]),
             check.names=FALSE)
      )}
  }, sanitize.text.function = identity)
    
}

# custom functions

drawNumTable <- function(nm, A, rn = FALSE, cn = FALSE,
                         step = 0.1, min = 0, max=1) {
  
  renderTable({
    inputMat <- paste0("<input id='", nm, 1:length(A),
                       "' class='shiny-bound-input'",
                       " type='number'",
                       " value='", as.vector(A),"'",
                       " step='", step, "'",
                       " min='", min, "'",
                       " max='", max, "'",
                       " style='width: 60px;'>")
    dim(inputMat) <- dim(A)
    rownames(inputMat) <- rownames(A)
    colnames(inputMat) <- colnames(A)
    inputMat
  }, sanitize.text.function = function(x) x,
  colnames = cn, rownames = rn)
}



# Run the application 
shinyApp(ui = ui, server = server)

