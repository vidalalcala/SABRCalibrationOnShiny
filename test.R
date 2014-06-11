library(shiny)
library(ggplot2)
library(data.table)
source("SABR.R")

# parametros
minVol <- 1
maxVol <- 100000
optionQuotes <- read.csv("data/smile.csv")
optionQuotesClean <- CleanSmile(optionQuotes, minVol, maxVol )
strike <- optionQuotesClean$Strike
price <- optionQuotesClean$Mid
plot(strike,price, type="l" )