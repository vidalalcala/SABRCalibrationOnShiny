# CONSTANT
EPS <- 10^(-8)

# sub function for SABR BS-IV
.x <- function(z, r){
  return(log((sqrt(1-2*r*z+z^2)+z-r)/(1-r)))
}

.dx <- function(z, r){
  return( ( 1.0 + (z-r)/(2*sqrt(1-2*r*z+z^2)) )/(sqrt(1-2*r*z+z^2)+z-r) )
}

.z <- function(f, K, a, nu){(nu/a)*log(f/K)}

#implied volatility
d1 <- function(s,f,K,tau){(log(f/K)+0.5*s^2*tau)/(s*sqrt(tau))}
d2 <- function(s,f,K,tau){(log(f/K)-0.5*s^2*tau)/(s*sqrt(tau))}
C <- function(s,f,K,tau) {f*pnorm(d1(s,f,K,tau))-K*pnorm(d2(s,f,K,tau))}
dC <- function(s,f,K,tau) {K*dnorm(d2(s,f,K,tau),0,1)*sqrt(tau)}
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
.t2  <- function(x){1.0/(1.0+exp(x)) - 1.0 }
.t2inv <- function(r){ log( - r/( 1.0 + r) ) }

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
SABR.W2 <- function(tau, f, K, a, rho){
  y <- (log(f/K)-0.5*a*a*tau)/a
  result1 <- -(1/3) * a * (tau^(2) * .N2(y/sqrt(tau)))
  result1 <- result1 + (1/3) * tau * y * .N2(y/sqrt(tau))
  result1 <- result1 - (1/6) * tau^(3/2) * .N1(y/sqrt(tau))
  result1 <- -0.5 * a * K * result1
  result2 <- -(1/4) * a * (tau^(2) * .N4(y/sqrt(tau)))
  result2 <- result2 + (1/4) * tau * y * .N4(y/sqrt(tau))
  result2 <- result2 - (1/4) * tau^(3/2) * .N3(y/sqrt(tau))
  result2 <- 0.5 * rho * rho * a * K * result2
  return(result1 + result2)
}

# Black-Scholes IV apporoximation formula by Hagan(2002)
SABR.HaganIV <- function(t, f, K, a, b, r, n)
{
  z <- .z(f, K, a, n)
  x <- .x(z, r)
  numerator   <- 1 + ((1-b)^2/24*a^2/(f*K)^(1-b) + 0.25*r*b*n*a/(f*K)^(0.5*(1-b)) + (2-3*r^2)*n^2/24)*t
  denominator <- x*(f*K)^(0.5*(1-b))*(1 + (1-b)^2/24*(log(f/K))^2 + (1-b)^4/1920*(log(f/K))^4)
  ifelse(abs((f-K)/f) < EPS, a*numerator/f^(1-b), z*a*numerator/denominator)
}

SABR.HaganDelta <- function(t, f, K, a, b, r, n)
{ result <- SABR.Black(t, f + EPS, K, SABR.HaganIV(t, f + EPS , K, a, b, r, n))
  result <- result - SABR.Black(t, f , K, SABR.HaganIV(t, f , K, a, b, r, n))
  result <- result/EPS
  return(result)
}

# Expansion with error O(nu^3)
SABR.W <- function(tau, f, K, nu, a, rho){
  return(SABR.Black(tau, f, K, a) + nu * SABR.W1(tau, f, K, a, rho) + nu*nu*SABR.W2(tau, f, K, a, rho) ) 
}

SABR.iv <- function(tau, f, K, nu, a, rho){
  Nquotes <-length(K)
  W <- SABR.W(tau, f, K, nu, a, rho)
  iv <- c()
  for (i in 1:Nquotes){
    iv[i] <- .ImpliedVolatility( f , K[i], tau , W[i] , 100 , 0.30)
  }    
  return(iv)
}


# Parameter calibration function for SABR
SABR.calibration <- function(tau, f, K, iv)
{
  # objective function for optimization
  # variables are transformed because of satisfing the constraint conditions
  
  #objective <- function(x){return(sum( iv - SABR.iv(tau, f, K, exp(x[1]), exp(x[2]), .t2(x[3])) )^2 ) }
  objective <- function(x){return(sum( ( iv - SABR.HaganIV(tau, f, K, exp(x[2]), exp(x[4]), .t2(x[3]) , exp(x[1])  ) )^2 ) ) }
  
  x <- optim(c(log(1.17), log(0.30), .t2inv(-0.30),  0.00 ), objective, list(maxit = 1000))

  # return the optimized parameters
  parameter <- x$par
  
  
  parameter <- c(exp(parameter[1]),exp(parameter[2]),.t2(parameter[3]), exp(parameter[4]) )
  names(parameter) <- c("nu", "alpha", "rho" , "beta")
  return(parameter)
}

# Parameter calibration function for SABR
nuSABR.calibration <- function(tau, f, K, iv)
{
  # objective function for optimization
  # variables are transformed because of satisfing the constraint conditions
  
  objective <- function(x){return(sum( iv - SABR.iv(tau, f, K, exp(x[1]), exp(x[2]), .t2(x[3])) )^2 ) }
  x <- optim(c(log(1.17), log(0.30), .t2inv(-0.30)), objective, list(maxit = 1000))
  
  # return the optimized parameters
  parameter <- x$par
  
  
  parameter <- c(exp(parameter[1]),exp(parameter[2]),.t2(parameter[3]))
  names(parameter) <- c("nu", "alpha", "rho" )
  return(parameter)
}


# clean smile data
CleanSmile <- function(smile, minVol , maxVol)
{ 
  smile <- smile[(smile$Vol>minVol & smile$Vol<maxVol),]
  # collapse rows with same strike and weight according to volume
  smileTable <- data.table(smile)
  setkey(smileTable,Strike)
  smileTable <- smileTable[, list(Last=sum(Last*Vol)/sum(Vol), Bid=sum(Bid*Vol)/sum(Vol), Ask=sum(Ask*Vol)/sum(Vol), Mid=sum(0.5*(Ask+Bid)*Vol)/sum(Vol), Vol=sum(Vol)), by=Strike]  
}

