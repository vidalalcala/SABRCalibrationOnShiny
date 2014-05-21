library(shiny)
library(ggplot2)
source("SABR.R")

shinyServer(function(input, output) {
  
  output$distPlot <- renderPlot({
    x <- calculateW()
    print(ggplot(data=x$data, aes(x=Strike, y=W, colour=Tag))
          + geom_point(size=4)+geom_line(size=1) + theme_grey(base_size=24))
  })
  calculateW <- reactive({
    forward  <- input$forward
    maturity <- input$maturity
    x <- input$marketData
    # convert string(ex:"12,0.346\n15,0.28\n17,0.243...") to data.frame
    x <- t(sapply(unlist(strsplit(x, "\n")), function(x) as.numeric(unlist(strsplit(x, ",")))))
    strike    <- x[,1]
    iv.market <- x[,2]
    SABR.parameter <- SABR.calibration(maturity, forward, strike, iv.market)
    W.model <- SABR.W(
      maturity, forward, strike, SABR.parameter[1], SABR.parameter[2], SABR.parameter[3])
    #
    list(
      parameter=SABR.parameter,
      data=rbind(
        data.frame(Strike=strike, W=W.model,  Tag="SABRnu"),
        data.frame(Strike=strike, W=SABR.Black(maturity,forward,strike,iv.market), Tag="Market")
      )
    )
  })
  output$summary <- renderPrint({
    x <- calculateW()
    print(x$parameter)
  })
})
