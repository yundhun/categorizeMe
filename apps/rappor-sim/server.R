library(shiny)
library(dplyr)
source("../../analysis/R/decode.R")
source("../../analysis/R/simulation.R")
source("../../analysis/R/encode.R")
source("../../analysis/R/categorize.R")
options(scipen=999)
#set.seed(123)
Plot <- function(x, color = "grey") {
  n <- nrow(x)
  if (n < 16) {
    par(mfrow = c(n, 1), mai = c(0, .5, .5, 0))
  } else if (n < 64) {
    par(mfrow = c(n / 2, 2), mai = c(0, .5, .5, 0))
  } else {
    par(mfrow = c(n / 4, 4), mai = c(0, .5, .5, 0))
  }
  for (i in 1:nrow(x)) {
    barplot(x[i, ], main = paste0("Cohort ", i), col = color, border = color)
  }
}

shinyServer(function(input, output) {
  # Example state global variable.
  es <- list()

  # Example buttons states.
  ebs <- rep(0, 3)
  
  Params <- reactive({
    list(k = as.numeric(input$size),
         h = as.numeric(input$hashes),
         m = as.numeric(input$instances),
         p = as.numeric(input$p),
         q = as.numeric(input$q),
         f = as.numeric(input$f))
  })
  
  PopParams <- reactive({
    infile <- input$datafile
    filedata <- read.csv(infile$datapath, stringsAsFactors = F)
    filedata <- filedata[,1]
    nsamp <- as.numeric(input$nsamp)
    nct <- as.numeric(input$nstrs)
    if(nsamp==0){
      each_size <- ( max(filedata) - min(filedata) ) / nct
      ctz <- seq(min(filedata),max(filedata),each_size)
      x <- min(filedata):max(filedata)
      plot(x=x, rep(1, length(x) ), type='n', main = 'Recommanded Categorization')
      abline(v=ctz,col='blue')
    }else{
      sampdata <- sample(filedata, nsamp)
      ctz <- A_categorize(sampdata, nct, resamp_cnt=100, plot=F)
    } 
    
    samp <- NULL
    for(x in filedata){
      j <- max(1, which(ctz <= x))
      samp <- c(samp, paste0('V_',j))
    }
    #print(length(samp))
    
    list(length(ctz),#nct,
      1, #as.numeric(input$nonzero),
      'Exponential',#input$decay,
      10,#as.numeric(input$expo),
      1.05,#as.numeric(input$background),
      nsamp,
      filedata,
      ctz,
      samp
      )
  })

  DecodingParams <- reactive({
    list(as.numeric(input$alpha),
         input$correction)
  })

  Sample <- reactive({
    #N <- input$N
    N <- length(PopParams()[[9]])
    params <- Params()
    pop_params <- PopParams()
    decoding_params <- DecodingParams()
    prop_missing <- input$missing
    fit <- GenerateSamples(N, params, pop_params,
                    alpha = decoding_params[[1]],
                    correction = decoding_params[[2]],
                    prop_missing = prop_missing,
                    samp = PopParams()[[9]])
    fit
  })
  
  # Results summary.
  output$pr <- renderTable({
    Sample()$summary
  },
                           include.rownames = FALSE, include.colnames = FALSE)

  # Results table.
  output$tab <- renderDataTable({
     Sample()$fit
   },
                                options = list(iDisplayLength = 100))

  # Epsilon.
  output$epsilon <- renderTable({
    Sample()$privacy
  },
                                include.rownames = FALSE, include.colnames = FALSE, digits = 4)

  # True distribution.
  output$probs <- renderPlot({
    #samp <- Sample()
    #probs <- samp$probs
    #detected <- match(samp$fit[, 1], samp$strs)
    #detection_frequency <- samp$privacy[7, 2]
    #PlotPopulation(probs, detected, detection_frequency)
    params <- PopParams()
    nct <- params[[1]]
    nsamp <- params[[6]]
    filedata <- params[[7]]
    ctz <- params[[8]]

    if(nsamp==0){
      each_size <- ( max(filedata) - min(filedata) ) / nct
      ctz <- seq(min(filedata),max(filedata),each_size)
      x <- min(filedata):max(filedata)
      plot(x=x, rep(1, length(x) ), type='n', main = 'Recommanded Categorization')
      abline(v=ctz,col='blue')
    }else{
      sampdata <- sample(filedata, nsamp)
      A_categorize(sampdata, nct, resamp_cnt=100)
    }
    
  })

  # True bits patterns.
  output$truth <- renderPlot({
    truth <- Sample()$truth
    Plot(truth[, -1, drop = FALSE], color = "darkblue")
  })

  # Lasso plot.
  output$lasso <- renderPlot({
    fit <- Sample()$lasso
    if (!is.null(fit)) {
      plot(fit)
    }
  })

  output$resid <- renderPlot({
    resid <- Sample()$residual
    params <- Params()
    plot(resid, xlab = "Bloom filter bits", ylab = "Residuals")
    abline(h = c(-1.96, 1.96), lty = 2, col = 2)
    sq <- qnorm(.025 / length(resid))
    abline(h = c(sq, -sq), lty = 2, col = 3, lwd = 2)
    abline(h = c(-3, 3), lty = 2, col = 4, lwd = 2)
    abline(v = params$k * (0:params$m), lty = 2, col = "blue")
    legend("topright", legend = paste0("SD = ", round(sd(resid), 2)), bty = "n")
  })

  # Estimated bits patterns.
  output$ests <- renderPlot({
    ests <- Sample()$ests
    Plot(ests, color = "darkred")
  })

  # Estimated vs truth.
  output$ests_truth <- renderPlot({
    plot(unlist(Sample()$ests), unlist(Sample()$truth[, -1]),
         xlab = "Estimates", ylab = "Truth", pch = 19)
    abline(0, 1, lwd = 4, col = "darkred")
  })
  
  output$ests_truth_hist <- renderPlot({
    ctz <- PopParams()[[8]]
    
    tdat <- Sample()$fit
    tdat$string <- as.numeric(gsub("V_", "", tdat$string))
    tdat <- tdat %>% arrange(string)

    decoded_v <- NULL
    for(idx in 1:nrow(tdat)){
      f_v <- as.double( tdat$string[idx]+1 )
      s_v <- as.double( tdat$string[idx] )
      c_v <- as.numeric( tdat$estimate[idx] )
      if( f_v > length(ctz) ) break;
      #decoded_v <- c(decoded_v, rep( mean( ctz[f_v] , ctz[s_v] ), c_v ) )
      decoded_v <- c(decoded_v, sample(  ctz[f_v]:ctz[s_v], c_v , replace=T) )
    }
   
    err_1 <- abs( floor(( (sum(PopParams()[[7]]) - sum(decoded_v)) / sum(PopParams()[[7]]) ) * 100 ) )
    par(mfrow=c(1,2))
    hist(PopParams()[[7]], probability = T, breaks = 100, main='Original Data')
    lines(density(PopParams()[[7]]),col='blue',lwd=2)
    lines(density(decoded_v),col='red',lwd=2)
    hist(decoded_v, probability = T, breaks = 100, xlim = c(min(PopParams()[[7]]),max(PopParams()[[7]])), 
         main=c(
           'Decoded Data',
           paste('Difference rate:',err_1,'%'))
         )
    lines(density(decoded_v),col='red',lwd=2)
  })  

  output$example <- renderPlot({
    params <- Params()
    strs <- Sample()$strs
    map <- Sample()$map
    samp <- Sample()

    # First run on app start.
    value <- sample(strs, 1)
    res <- Encode(value, map, strs, params, N = PopParams()[[9]])

    if (input$new_user > ebs[1]) {
      res <- Encode(es$value, map, strs, params, N = PopParams()[[9]])
      ebs[1] <<- input$new_user
    } else if (input$new_value > ebs[2]) {
      res <- Encode(value, map, strs, params, cohort = es$cohort, id = es$id,
                    N = PopParams()[[9]])
      ebs[2] <<- input$new_value
    } else if (input$new_report > ebs[3]) {
      res <- Encode(es$value, map, strs, params, B = es$B,
                    BP = es$BP, cohort = es$cohort, id = es$id, N = PopParams()[[9]])
      ebs[3] <<- input$new_report
    }
    es <<- res
    ExamplePlot(res, params$k, c(ebs, input$new_user, input$new_value, input$new_report))
  })

})
