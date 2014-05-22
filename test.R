library(testthat)
source("SABR.R")
test_that("SBARnu model test", {
  k <- c(12.0, 15.0, 17.0, 19.5, 20.0, 22.0, 22.5,24.5, 25.0,27.0, 27.5, 29.5, 30.0, 32.0, 32.5, 34.5, 37.0)
  iv <- c(0.346,0.280, 0.243, 0.208, 0.203, 0.192, 0.192, 0.201, 0.205, 0.223, 0.228, 0.247, 0.252, 0.271, 0.275, 0.293, 0.313) 
  f <- 22.724
  tau <- 0.583
  a <- 0.317
  rho <- 0.111
  nu <- 0.10
  W.model <- SABR.W(tau, f, k, a, rho, nu)
  params <- SABR.calibration(tau, f, k, iv)
  W.calibrated <- SABR.W(tau, f, k, params[1], params[2], params[3])
  W.market <- SABR.Black(tau,f,k,iv)
  # check wether initial model can produce market IV or not.
  for(i in length(k)){expect_equal(W.model[i],      W.market[i], tolerance=0.001*V.market[i])}
  # check wether calibrated parameter can produce market IV or not.
  for(i in length(k)){expect_equal(W.calibrated[i], W.market[i], tolerance=0.001*V.market[i])}
})