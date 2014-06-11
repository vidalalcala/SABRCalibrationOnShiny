library(shiny)
library(ggplot2)
library(data.table)
source("SABR.R")

# START code that runs once
optionQuotes <- read.csv("data/smile.csv")
# END 

shinyServer(function(input, output) {
  
  output$distPlotHagan <- renderPlot({
    x <- calculateIV()
    print(ggplot(data=x$dataPlotHagan, aes(x=Strike, y=IV, colour=Tag))
          + geom_point(size=4)+geom_line(size=1) + theme_grey(base_size=24))
  })
  output$deltaPlotHagan <- renderPlot({
    x <- calculateIV()
    print(ggplot(data=x$dataHedgeHagan, aes(x=Strike, y=delta))
          + geom_point(size=4)+geom_line(size=1) + theme_grey(base_size=24))
  })
  
  output$deltaTableHagan <- renderTable({
    x <-calculateIV()
    x$dataHedgeHagan
  })
  
  output$summaryHagan <- renderPrint({
    x <- calculateIV()
    print(x$parameterHagan)
  })
  
  output$distPlotSABR <- renderPlot({
    x <- calculateIV()
    print(ggplot(data=x$dataPlotSABR, aes(x=Strike, y=IV, colour=Tag))
          + geom_point(size=4)+geom_line(size=1) + theme_grey(base_size=24))
  })
  output$deltaPlotSABR <- renderPlot({
    x <- calculateIV()
    print(ggplot(data=x$dataHedgeSABR, aes(x=Strike, y=delta))
          + geom_point(size=4)+geom_line(size=1) + theme_grey(base_size=24))
  })
  
  output$deltaTableSABR <- renderTable({
    x <-calculateIV()
    x$dataHedgeSABR
  })
  
  output$summarySABR <- renderPrint({
    x <- calculateIV()
    print(x$parameterSABR)
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
    price <- optionQuotesClean$Mid
    
    #Calculate implied vols
    Nquotes <-length(strike)
    iv.market <- c()
    sigmaStart <- 0.60
    for (i in 1:Nquotes){
      iv.market[i] <- .ImpliedVolatility( forward , strike[i], maturity , exp(r*maturity)*optionQuotesClean$Mid[i] , 1 , sigmaStart)
    }
    
    Hagan.parameter <- Hagan.calibration(maturity, forward, strike, iv.market)
    SABR.parameter <- SABR.calibration(maturity, forward, strike, price,r)
    IV.Hagan <- Hagan.IV(
      maturity, forward, strike, Hagan.parameter[2], Hagan.parameter[4] ,  Hagan.parameter[3], Hagan.parameter[1])
    IV.SABR <- SABR.iv(
      maturity, forward, strike, SABR.parameter[1], SABR.parameter[2], SABR.parameter[3])
    delta.Hagan <- Hagan.Delta(
      maturity, forward, strike, Hagan.parameter[2], Hagan.parameter[4] ,  Hagan.parameter[3], Hagan.parameter[1])
    delta.SABR <- exp(-r*maturity)*SABR.Delta(
      maturity, forward, strike, SABR.parameter[1], SABR.parameter[2],  SABR.parameter[3])
    #
    list(
      parameterSABR=SABR.parameter ,
      parameterHagan=Hagan.parameter ,
      dataPlotSABR=rbind(
        data.frame(Strike=strike, IV=iv.market, Tag="Market"),
        data.frame(Strike=strike, IV=IV.SABR, Tag="nuSABR")
      ),
      dataPlotHagan=rbind(
        data.frame(Strike=strike, IV=iv.market, Tag="Market"),
        data.frame(Strike=strike, IV=IV.Hagan, Tag="SABRHagan")
      ),
      dataHedgeSABR=data.frame(Strike=strike, delta=delta.SABR),
      dataHedgeHagan=data.frame(Strike=strike, delta=delta.Hagan)
    )
  })
})
