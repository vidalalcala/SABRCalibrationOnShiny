# CONSTANT
EPS <- 10^(-8)

# sub functions for SABR \nu expansion

.N <- function(x){pnorm(x)}
.N1 <- function(x){dnorm(x)}
.N2 <- function(x){-x*dnorm(x)}
.N3 <- function(x){(x^2-1)*dnorm(x)}
.N4 <- function(x){(-x^3+3*x)*dnorm(x)}

SABR.Black <- function(tau, f, K, a){
  y <- (log(f/K)-0.5*a*a*tau)/a
  f*.N(y/sqrt(tau) + a*sqrt(tau))-K*(.N(y/sqrt(tau)))
} 

SABR.W1 <- function(tau, f, K, a, rho){
  y <- (log(f/K)-0.5*a*a*tau)/a
  -0.5*rho*a*K*tau*(.N2(y/sqrt(tau)))
}

# Expansion with error O(nu^2)
SABR.W <- function(tau, f, K, nu, a, rho){
  SABR.Black(tau, f, K, a) + nu * SABR.W1(tau, f, K, a, rho) 
}

# Parameter calibration function for SABRnu
SABR.calibration <- function(t, f, K, iv)
{
  # objective function for optimization
  # variables are transformed because of satisfing the constraint conditions
  objective <- function(x){sum( (SABR.Black(tau,f,K,iv) - SABR.W(tau, f, K, x[1], x[2], x[3]))^2 )}
  x <- nlm(objective, c(0.0, 0.10, 0.0))
  # return the optimized parameters
  parameter <- x$estimate
  names(parameter) <- c("nu", "alpha", "rho")
  parameter
}

