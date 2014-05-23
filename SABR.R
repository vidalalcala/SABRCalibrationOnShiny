# CONSTANT
EPS <- 10^(-8)

# variable transformation function
.t2  <- function(x){2/(1+exp(x)) -1}
.t2inv <- function(r){ log((1-r)/(1+r)) }

# sub functions for SABR \nu expansion

.N <- function(x){pnorm(x)}
.N1 <- function(x){dnorm(x)}
.N2 <- function(x){-x*dnorm(x)}
.N3 <- function(x){(x^2-1)*dnorm(x)}
.N4 <- function(x){(-x^3+3*x)*dnorm(x)}

SABR.Black <- function(tau, f, K, a){
  y <- (log(f/K)-0.5*a*a*tau)/a
  return(f*.N(y/sqrt(tau) + a*sqrt(tau))-K*(.N(y/sqrt(tau))))
} 

SABR.W1 <- function(tau, f, K, a, rho){
  y <- (log(f/K)-0.5*a*a*tau)/a
  return(-0.5*rho*a*K*tau*(.N2(y/sqrt(tau))))
}

# Expansion with error O(nu^2)
SABR.W <- function(tau, f, K, nu, a, rho){
  return(SABR.Black(tau, f, K, a) + nu * SABR.W1(tau, f, K, a, rho) )
}

# Parameter calibration function for SABRnu
SABR.calibration <- function(tau, f, K, iv)
{
  # objective function for optimization
  # variables are transformed because of satisfing the constraint conditions
  objective <- function(x){sum( ( SABR.Black(tau,f,K,iv) - SABR.W(tau, f, K, exp(x[1]), exp(x[2]), .t2(x[3])) )^2 ) }
  x <- optim(c(log(1.00), log(0.20), .t2inv(0.00)), objective,control = list( "maxit" = 100000000) )
  cat("\n"," correction: ",SABR.W1(tau, f, K, exp(x$par[2]), .t2(x$par[3])))
  cat("\n"," price: ",SABR.W(tau, f, K, exp(x$par[1]), exp(x$par[2]), .t2(x$par[3])))
  # return the optimized parameters
  parameter <- x$par
  parameter <- c(exp(parameter[1]),exp(parameter[2]),.t2(parameter[3]))
  cat("\n code: ", x$convergence)
  names(parameter) <- c("nu", "alpha", "rho")
  cat("\n" ,parameter)
  return(parameter)
}

