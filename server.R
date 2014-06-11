library(shiny)
library(ggplot2)
library(data.table)
source("SABR.R")

# START code that runs once
optionQuotes <- read.csv("data/smile.csv")
# END 

shinyServer(function(input, output) {
  
  output$distPlot <- renderPlot({
    x <- calculateIV()
    print(ggplot(data=x$data, aes(x=Strike, y=IV, colour=Tag))
          + geom_point(size=4)+geom_line(size=1) + theme_grey(base_size=24))
  })
  output$deltaPlot <- renderPlot({
    x <- calculateIV()
    print(ggplot(data=x$dataHedge, aes(x=Strike, y=delta))
          + geom_point(size=4)+geom_line(size=1) + theme_grey(base_size=24))
  })
  
  output$deltaTable <- renderTable({
    x <-calculateIV()
    x$dataHedge
  })
  
  calculateIV <- reactive({
    forward  <- input$forward
    maturity <- input$maturity
    r <- input$r
    #beta <- input$beta
    minVol <- input$minVol
    maxVol <- input$maxVol
    
    optionQuotesClean <- CleanSmile(optionQuotes, minVol, maxVol )
    strike <- optionQuotesClean$Strike
    
    #Calculate implied vols
    Nquotes <-length(strike)
    iv.market <- c()
    sigmaStart <- 0.60
    for (i in 1:Nquotes){
      iv.market[i] <- .ImpliedVolatility( forward , strike[i], maturity , exp(r*maturity)*optionQuotesClean$Mid[i] , 1 , sigmaStart)
    }
    print(iv.market)
    SABR.parameter <- SABR.calibration(maturity, forward, strike, iv.market)
    IV.model <- SABR.iv(
      maturity, forward, strike, SABR.parameter[1], SABR.parameter[2], SABR.parameter[3])
    IV.Hagan <- SABR.HaganIV(
      maturity, forward, strike, SABR.parameter[2], SABR.parameter[4] ,  SABR.parameter[3], SABR.parameter[1])
    delta.Hagan <- SABR.HaganDelta(
      maturity, forward, strike, SABR.parameter[2], SABR.parameter[4] ,  SABR.parameter[3], SABR.parameter[1])
    #
    list(
      parameter=SABR.parameter ,
      data=rbind(
        #data.frame(Strike=strike, IV=IV.model,  Tag="SABRnu"),
        data.frame(Strike=strike, IV=iv.market, Tag="Market"),
        data.frame(Strike=strike, IV=IV.Hagan, Tag="Hagan")
      ),
      dataHedge=data.frame(Strike=strike, delta=delta.Hagan)
    )
  })
  output$summary <- renderPrint({
    x <- calculateIV()
    print(x$parameter)
  })
})
