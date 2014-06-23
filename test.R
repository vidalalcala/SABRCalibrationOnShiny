library(shiny)
library(ggplot2)
library(data.table)
source("SABR.R")

# parametros
strike <- seq(1, 200, 1)
vol <- seq(0.001, 10, 0.001)
f <- 100.0
maturity <- 1.00
r <- 0.000
nu <- 1.0000
alpha <- 0.20
rho <- -0.75
minSigma <- 0.001
maxSigma <- 100
price <- SABR.W(maturity, f, strike, nu, alpha, rho)
iv <- SABR.iv(maturity, f, strike, nu, alpha, rho, minSigma, maxSigma)

plot(strike, .Black(maturity, f, strike, iv),type="l")
plot.new()
plot(price,iv,type="l")
plot.new()
plot(vol,.Black(maturity, f, 5, vol),type="l")
plot.new()
plot(strike, iv, type="l" )
print(length(strike))
