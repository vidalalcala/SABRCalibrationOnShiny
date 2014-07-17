library(shiny)
library(ggplot2)
library(data.table)
source("SABR.R")

# parametros
strike <- seq(1, 200, 1)
vol <- seq(0.001, 10, 0.001)
f <- 100.0
maturity <- 1.0
r <- 10.10
nu <- 1.00
alpha <- 0.20
rho <- -0.75
minSigma <- 0.001
maxSigma <- 1.0
price <- SABR.W(maturity, f, strike, nu, alpha, rho)
iv <- SABR.iv(maturity, f, strike, nu, alpha, rho, minSigma, maxSigma)
ivHagan <- Hagan.IV(maturity, f, strike, alpha, 1.0, rho, nu)

#data from http://www.wilmott.com/messageview.cfm?catid=34&threadid=94086
ivAlan <- c(0.4065 , 0.2881, 0.1934, 0.1488, 0.1662, 0.1899, 0.2114)
strikeAlan<- c(50, 75, 100, 125, 150, 175, 200)


plot(strike, .Black(maturity, f, strike, ivHagan),type="l")
plot.new()
plot(vol,.Black(maturity, f, 5, vol),type="l")
plot.new()
plot(strike, iv, type="l" )
lines(strikeAlan, ivAlan, col="blue")
lines(strike, ivHagan, col="red")

print(length(strike))
