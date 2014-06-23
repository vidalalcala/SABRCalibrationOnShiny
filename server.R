library(shiny)
library(ggplot2)
library(data.table)
source("SABR.R")

# START code that runs once
optionQuotes <- read.csv("data/smile.csv")
# END 

shinyServer(function(input, output) {
  
  # output definition
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
  output$pricesPlotSABR <- renderPlot({
    x <- calculateIV()
    print(ggplot(data=x$dataPlotSABRPrices, aes(x=Strike, y=Price, colour=Tag))
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
  
  # multifractal output
  
  output$distPlotFractal <- renderPlot({
    x <- calculateIV()
    print(ggplot(data=x$dataPlotFractal, aes(x=Strike, y=IV, colour=Tag))
          + geom_point(size=4)+geom_line(size=1) + theme_grey(base_size=24))
  })
  output$pricesPlotFractal <- renderPlot({
    x <- calculateIV()
    print(ggplot(data=x$dataPlotFractalPrices, aes(x=Strike, y=Price, colour=Tag))
          + geom_point(size=4)+geom_line(size=1) + theme_grey(base_size=24))
  })
  
  output$deltaPlotFractal <- renderPlot({
    x <- calculateIV()
    print(ggplot(data=x$dataHedgeFractal, aes(x=Strike, y=delta))
          + geom_point(size=4)+geom_line(size=1) + theme_grey(base_size=24))
  })
  
  output$deltaTableFractal <- renderTable({
    x <-calculateIV()
    x$dataHedgeFractal
  })
  
  output$summaryFractal <- renderPrint({
    x <- calculateIV()
    print(x$parameterFractal)
  })
  
  # the function that calibrates parameters to data
  calculateIV <- reactive({
    # read inputs
    forward  <- input$forward
    maturity <- input$maturity
    r <- input$r
    minVol <- input$minVol
    maxVol <- input$maxVol
    minSigma <- input$minSigma
    maxSigma <- input$maxSigma
    optionQuotesClean <- CleanSmile(optionQuotes, minVol, maxVol )
    strike <- optionQuotesClean$Strike
    price <- optionQuotesClean$Mid
    parameter0 <- data.frame( nu=input$nu0, alpha=input$alpha0, rho=input$rho0, beta=input$beta0, h=input$h0, phi=input$phi0, sig=input$sig0)
    
    #Calculate market implied vols
    Nquotes <-length(strike)
    iv.market <- c()
    for (i in 1:Nquotes){
      iv.market[i] <- .ImpliedVolatilityNewton( maturity, forward , strike[i], exp(r*maturity)*optionQuotesClean$Mid[i] , 10 , parameter0$alpha, minSigma, maxSigma)
    }

    # calibrate models
    Hagan.parameter <- Hagan.calibration(maturity, forward, strike, iv.market, parameter0)
    SABR.parameter <- SABR.calibration(maturity, forward, strike, iv.market, parameter0, minSigma, maxSigma)
    fractal.parameter <- fractal.calibration(maturity, forward, strike, iv.market, parameter0, minSigma, maxSigma)

    # calculate model implied vol
    IV.Hagan <- Hagan.IV(
      maturity, forward, strike, Hagan.parameter[2], Hagan.parameter[4] ,  Hagan.parameter[3], Hagan.parameter[1])
    IV.SABR <- SABR.iv(
      maturity, forward, strike, SABR.parameter[1], SABR.parameter[2], SABR.parameter[3], minSigma, maxSigma)
    IV.fractal <- fractal.iv(
      maturity, forward, strike, fractal.parameter[1], fractal.parameter[2], fractal.parameter[3], fractal.parameter[4], minSigma, maxSigma)
    
    # calculate model hedge
    delta.Hagan <- Hagan.Delta(
      maturity, forward, strike, Hagan.parameter[2], Hagan.parameter[4] ,  Hagan.parameter[3], Hagan.parameter[1])
    delta.SABR <- exp(-r*maturity)*SABR.Delta(
      maturity, forward, strike, SABR.parameter[1], SABR.parameter[2],  SABR.parameter[3])
    delta.fractal <- exp(-r*maturity)*fractal.Delta(
      maturity, forward, strike, fractal.parameter[1], fractal.parameter[2],  fractal.parameter[3], fractal.parameter[4])
    
    # function output
    list(
      parameterSABR=SABR.parameter ,
      parameterFractal=fractal.parameter ,
      parameterHagan=Hagan.parameter ,
      dataPlotSABR=rbind(
        data.frame(Strike=strike, IV=iv.market, Tag="Market"),
        data.frame(Strike=strike, IV=IV.SABR, Tag="nuSABR")
      ),
      dataPlotSABRPrices=rbind(
        data.frame(Strike=strike, Price=optionQuotesClean$Mid, Tag="Market"),
        data.frame(Strike=strike, Price=exp(-r*maturity)*SABR.W(maturity, forward, strike, SABR.parameter[1], SABR.parameter[2], SABR.parameter[3]), Tag="nuSABR")
      ),
      dataPlotFractal=rbind(
        data.frame(Strike=strike, IV=iv.market, Tag="Market"),
        data.frame(Strike=strike, IV=IV.fractal, Tag="multifractal")
      ),
      dataPlotFractalPrices=rbind(
        data.frame(Strike=strike, Price=optionQuotesClean$Mid, Tag="Market"),
        data.frame(Strike=strike, Price=exp(-r*maturity)*fractal.W(maturity, forward, strike, fractal.parameter[1], fractal.parameter[2], fractal.parameter[3], fractal.parameter[4]), Tag="multifractal")
      ),
      dataPlotHagan=rbind(
        data.frame(Strike=strike, IV=iv.market, Tag="Market"),
        data.frame(Strike=strike, IV=IV.Hagan, Tag="SABRHagan")
      ),
      dataHedgeSABR=data.frame(Strike=strike, delta=delta.SABR),
      dataHedgeFractal=data.frame(Strike=strike, delta=delta.fractal),
      dataHedgeHagan=data.frame(Strike=strike, delta=delta.Hagan)
    )
  })
})
