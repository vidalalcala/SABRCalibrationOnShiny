library(shiny)
library(ggplot2)
library(data.table)
source("SABR.R")

# parametros
strike <- seq(50, 200, 1)
vol <- seq(0.001, 10, 0.001)
f <- 100.0
maturity <- 1.0
r <- 0.0015
lambda <- 0.03;
h <- lambda
alpha <- 0.30
phi <- 75
sig <- 0.30
minSigma <- 0.001
maxSigma <- 100
price <- fractal.W(maturity, f, strike, h, alpha, phi, sig)
iv <- fractal.iv(maturity, f, strike, h, alpha, phi, sig, minSigma, maxSigma)
plot(strike, ifelse(price>0, 1.0 , 0.0),type="l")
plot.new()
plot(strike, price,type="l")
plot.new()
plot(strike, iv, type="l" )
print(length(strike))
