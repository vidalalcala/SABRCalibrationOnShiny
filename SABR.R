# CONSTANT
EPS <- 10^(-8)

#implied volatility
d1 <- function(s,f,K,tau){(log(f/K)+0.5*s^2*tau)/(s*sqrt(tau))}
d2 <- function(s,f,K,tau){(log(f/K)-0.5*s^2*tau)/(s*sqrt(tau))}
C <- function(s,f,K,tau) {f*pnorm(d1(s,f,K,tau))-K*pnorm(d2(s,f,K,tau))}
dC <- function(s,f,K,tau) {f*dnorm(d1(s,f,K,tau),0,1)*sqrt(tau)}
.ImpliedVolatility <- function( f , K, tau , CM , N, sigma0){
  # check price
  if (CM > f){
    cat("<ImpliedVolatility> The call price is out of bounds")
  }
  
  sigma <- sigma0;
  for (i in 1:N){
    sigma = sigma - (C(sigma,f,K,tau) - CM )/dC(sigma,f,K,tau);  
  }
  return(sigma)
}


# variable transformation function
.t2  <- function(x){2.0/(1.0+exp(x)) - 1.0 }
.t2inv <- function(r){ log( (1.0 - r)/( 1.0 + r) ) }

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
  return(0.5*rho*a*K*tau*(.N2(y/sqrt(tau))))
}

# Expansion with error O(nu^2)
SABR.W <- function(tau, f, K, nu, a, rho){
  return(SABR.Black(tau, f, K, a) + nu * SABR.W1(tau, f, K, a, rho) )
}

SABR.iv <- function(tau, f, K, nu, a, rho){
  Nquotes <-length(K)
  W <- SABR.Black(tau, f, K, a) + nu * SABR.W1(tau, f, K, a, rho)
  iv <- c()
  for (i in 1:Nquotes){
    iv[i] <- .ImpliedVolatility( f , K[i], tau , W[i] , 5 , 0.50)
  }    
  return(iv)
}


# Parameter calibration function for SABRnu
SABR.calibration <- function(tau, f, K, iv)
{
  # objective function for optimization
  # variables are transformed because of satisfing the constraint conditions
  
  objective <- function(x){return(sum( ( iv - SABR.iv(tau, f, K, exp(x[1]), exp(x[2]), .t2(x[3])) )^2 ) )}
  x <- optim(c(log(0.10), log(0.20), .t2inv(-0.05)), objective,control = list( "maxit" = 10000) )

  # return the optimized parameters
  parameter <- x$par
  
  parameter <- c(exp(parameter[1]),exp(parameter[2]),.t2(parameter[3]))
  cat("\n <optim> code: ", x$convergence)
  names(parameter) <- c("nu", "alpha", "rho")
  return(parameter)
}

