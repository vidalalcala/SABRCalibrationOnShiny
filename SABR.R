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

# sub functions for SABR \nu expansion

.N <- function(x){pnorm(x)}
.N1 <- function(x){dnorm(x)}
.N2 <- function(x){-x*dnorm(x)}
.N3 <- function(x){(x^2-1)*dnorm(x)}
.N4 <- function(x){(-x^3+3*x)*dnorm(x)}
.N5 <- function(x){(x^4-6*x^2+3)*dnorm(x)}
.N6 <- function(x){(-x^5 + 10*x^3-15*x)*dnorm(x)}

#implied volatility

.Black <- function(tau, f, K, a){
  y <- (log(f/K)-0.5*a*a*tau)/a
  return(f*.N(y/sqrt(tau) + a*sqrt(tau))-K*(.N(y/sqrt(tau))))
}
.dBlack <- function(tau,f,K,a) {
  y <- (log(f/K)-0.5*a*a*tau)/a
  K*.N1(y/sqrt(tau))*sqrt(tau)
}

.ImpliedVolatilityNewton <- function( tau, f, K, WMarket , N, sigma0, minSigma , maxSigma){
  # check price
  if ( WMarket > f ){
    cat("<ImpliedVolatilityNewton> The call price is out of no hedge bounds, strike =", K,"\n")
    return(maxSigma)
  }
  if ( WMarket < ifelse(f-K>0, f-K, 0) ){
    cat("<ImpliedVolatilityNewton> The call price is out of no hedge bounds, strike =", K,"\n")
    return(minSigma)
  }
  
  sigma <- sigma0;
  for (i in 1:N){
    sigma = sigma - (.Black(tau ,f ,K , sigma) - WMarket )/.dBlack(tau,f,K,sigma);  
  }
  error <- abs(.Black( tau, f, K, sigma ) - WMarket)

  if( is.na(sigma) || (error > 0.001) || (is.nan(error))){
    print("<impliedvol> Newton fails !")
    sigma <- .ImpliedVolatilityBisection(tau, f, K, WMarket, minSigma, maxSigma )
    return(sigma)
  }
  else{
    return(sigma)
  }
}

.ImpliedVolatilityBisection <- function(tau, f, K, WM, minSigma, maxSigma){
    sig.up <- maxSigma
    sig.down <- minSigma
    sig <- 0.5*(sig.down + sig.up)
    if ( (.Black(tau,f,K,sig.down)-WM)*(.Black(tau,f,K,sig.up)-WM)>0){
      cat("<bisection> same sign at end points with strike = ", K, "\n")
      return(0.0)
    }
    count <- 0
    err <- .Black(tau,f,K,sig) - WM 
    
    ## repeat until error is sufficiently small or counter hits 1000
    while(abs(err) > 0.001 && count<1000){
      if((.Black(tau,f,K,sig.down)-WM)*(.Black(tau,f,K,sig)-WM)<0){
        sig.up <- sig
        sig <- (sig.down + sig.up)/2
      }else{
        sig.down <- sig
        sig <- (sig.down + sig.up)/2
      }
      err <- .Black(tau,f,K,sig) - WM
      count <- count + 1
    }
    
    ## return NA if counter hit 1000
    if(count > 1000){
      cat("<impliedvol> Bisection reaches 1000 iterations with strike ", K, "\n")
      return(sig)
    }else{
      return(sig)
    }
  }

SABR.W1 <- function(tau, f, K, a, rho){
  y <- (log(f/K)-0.5*a*a*tau)/a
  return(0.5*rho*a*K*tau*(.N2(y/sqrt(tau))))
}

SABR.W2 <- function(tau, f, K, a, rho){
  y <- (log(f/K)-0.5*a*a*tau)/a
  result1 <- (1/3) * a * (tau^(2) * .N2(y/sqrt(tau)))
  result1 <- result1 + (1/3) * tau * y * .N2(y/sqrt(tau))
  result1 <- result1 - (1/6) * tau^(3/2) * .N1(y/sqrt(tau))
  result1 <- -0.5 * a * K * result1
  result2 <- (1/4) * a * (tau^(2) * .N4(y/sqrt(tau)))
  result2 <- result2 + (1/4) * tau * y * .N4(y/sqrt(tau))
  result2 <- result2 - (1/4) * tau^(3/2) * .N3(y/sqrt(tau))
  result2 <- -0.5 * rho * rho * a * K * result2
  return(result1 + result2)
}

SABR.W3 <- function(tau, f, K, a, rho){
  y <- (log(f/K)-0.5*a*a*tau)/a
  result <- (1/24) * (rho*K*a*tau^2) * .N2(y/sqrt(tau))
  result <- result - (1/16) * rho*K*a*tau^(3/2) * ( y + (8/3)*a*tau ) * .N3(y/sqrt(tau))
  result <- result + (1/12) * rho*K*a*tau*( ( y + a*tau )^2 + (1/4)*tau + (1/5)*rho^2 * tau ) * .N4(y/sqrt(tau))
  result <- result - (1/60) * (rho^3)*K*a*tau^(3/2)* ( y + (5/2)*a*tau )*.N5(y/sqrt(tau))
  result <- result + (1/48) * (rho^3)*K*a*tau * ( ( y + a*tau )^2 + (1/5)*tau )*.N6(y/sqrt(tau))
  return(result)
}

fractal.W1 <- function(tau, f, K, a, phi, sig){
  y <- (log(f/K)-0.5*a*a*tau)/a
  result <- K*a*tau*.N1(y/sqrt(tau))*( 0.5*sqrt(tau)*( (-phi)*log(a/sig)+(0.5*a^2)*(1/sig) )-(a/sig)*(2/3) )
  return(result)
}

fractal.W2 <- function(tau, f, K, a, phi, sig){
  y <- (log(f/K)-0.5*a*a*tau)/a
  result <- (-1/sig)*a^2*K*( (-phi)*log(a/sig)+ (-phi) + ((3/2)*a^2)*(1/sig) )*((tau)^(-3/2))*(-y)*.N1(y/sqrt(tau))*((tau)^3)*(1/6) + 1/sig*a^3*K*4/(3*sig)*((tau)^(-3/2))*(-y)*.N1(y/sqrt(tau))*((2/5)*(tau)^(5/2)) +( (-phi)*a*log(a/sig) + (1/2*a^3)*1/sig)*K*((-phi)*log(a/sig)+(-phi)+(3/2*a^2)*1/sig)*((tau)^(-1/2))*.N1(y/sqrt(tau))*((tau)^3/6) -( (-phi)*a*log(a/sig) + (1/2*a^3)*1/sig )*K*4/(3*sig)*a*((tau)^(-1/2))*.N1(y/sqrt(tau))*((2/5)*(tau)^(5/2)) + 1/(2*sig^2)*a^3*K*((tau)^(-1/2))*.N1(y/sqrt(tau))*(((tau)^(-1))*((1/6)*tau^3) + (y^2)*(tau)/3) +1/(2*sig^2)*a^3*K*(1/sqrt(tau))*.N1(y/sqrt(tau))*((y^2)*(tau)+a*y*((tau)^2)/3-2*(y^2)*(tau)/2+(tau)^2/2)
  return(result)
}



fractal.W <- function(tau, f, K, h, a, phi, sig){
  return(.Black(tau, f, K, a) - h * fractal.W1(tau, f, K, a, phi, sig) - h^2 * fractal.W2(tau, f, K, a, phi, sig) )
}

fractal.Delta <- function(tau, f, K, h, a, phi, sig)
{ result <- fractal.W(tau, f + EPS, K, h, a, phi, sig)
  result <- result - fractal.W(tau, f, K, h, a, phi, sig)
  result <- result/EPS
  return(result)
}

fractal.iv <- function(tau, f, K, h, a, phi, sig, minSigma, maxSigma){
  Nquotes <-length(K)
  W <- fractal.W(tau, f, K, h, a, phi, sig)
  iv <- c()
  for (i in 1:Nquotes){
    iv[i] <- .ImpliedVolatilityNewton(tau, f , K[i] , W[i] , 10 , a , minSigma, maxSigma)
  }    
  return(iv)
}

# Expansion with error O(nu^3)
SABR.W <- function(tau, f, K, nu, a, rho){
  return(.Black(tau, f, K, a)  + nu * SABR.W1(tau, f, K, a, rho)   + nu*nu*SABR.W2(tau, f, K, a, rho)  + nu*nu*nu*SABR.W3(tau, f, K, a, rho)) 
}

SABR.Delta <- function(t, f, K, nu, a, rho)
{ result <- SABR.W(t, f + EPS, K, nu, a, rho)
  result <- result - SABR.W(t, f, K, nu, a, rho)
  result <- result/EPS
  return(result)
}

SABR.iv <- function(tau, f, K, nu, a, rho, minSigma, maxSigma){
  Nquotes <-length(K)
  W <- SABR.W(tau, f, K, nu, a, rho)
  iv <- c()
  for (i in 1:Nquotes){
    iv[i] <- .ImpliedVolatilityNewton(tau, f , K[i] , W[i] , 10 , a , minSigma, maxSigma)
  }    
  return(iv)
}

# variable transformation function
.t2  <- function(x){1.0/(1.0+exp(x)) - 1.0 }
.t2inv <- function(r){ log( - r/( 1.0 + r) ) }

# Parameter calibration function for SABR
SABR.calibration <- function(tau, f, K, iv, parameter0, minSigma, maxSigma)
{
  # objective function for optimization
  # variables are transformed because of satisfing the constraint conditions
  
  objective <- function(x){return( sum( (iv - SABR.iv(tau, f, K, exp(x[1]), exp(x[2]), .t2(x[3]), minSigma, maxSigma) )^2 ) ) }
  x <- optim(c(log(parameter0$nu), log(parameter0$alpha), .t2inv(parameter0$rho)), objective, list(maxit = 100000))
  
  # return the optimized parameters
  parameter <- x$par
  parameter <- c(exp(parameter[1]),exp(parameter[2]),.t2(parameter[3]))
  names(parameter) <- c("nu", "alpha", "rho" )
  return(parameter)
}

# Parameter calibration function for multifractal model
fractal.calibration <- function(tau, f, K, iv , parameter0, minSigma, maxSigma)
{
  # objective function for optimization
  # variables are transformed because of satisfing the constraint conditions
  
  objective <- function(x){return( sum( (iv - fractal.iv(tau, f, K, x[1], exp(x[2]), exp(x[3]), exp(x[4]), minSigma, maxSigma) )^2 ) ) }
  x <- optim(c(parameter0$h, log(parameter0$alpha), log(parameter0$phi), log(parameter0$sig)), objective, list(maxit = 1000000))
  
  # return the optimized parameters
  parameter <- x$par
  parameter <- c(parameter[1], exp(parameter[2]), exp(parameter[3]), exp(parameter[4]))
  names(parameter) <- c("h", "alpha", "phi", "sigma" )
  return(parameter)
}

# Black-Scholes IV apporoximation formula by Hagan(2002)
Hagan.IV <- function(t, f, K, a, b, r, n)
{
  z <- .z(f, K, a, n)
  x <- .x(z, r)
  numerator   <- 1 + ((1-b)^2/24*a^2/(f*K)^(1-b) + 0.25*r*b*n*a/(f*K)^(0.5*(1-b)) + (2-3*r^2)*n^2/24)*t
  denominator <- x*(f*K)^(0.5*(1-b))*(1 + (1-b)^2/24*(log(f/K))^2 + (1-b)^4/1920*(log(f/K))^4)
  ifelse(abs((f-K)/f) < EPS, a*numerator/f^(1-b), z*a*numerator/denominator)
}

Hagan.Delta <- function(t, f, K, a, b, r, n)
{ result <- .Black(t, f + EPS, K, Hagan.IV(t, f + EPS , K, a, b, r, n))
  result <- result - .Black(t, f , K, Hagan.IV(t, f , K, a, b, r, n))
  result <- result/EPS
  return(result)
}

# Parameter calibration function for SABR
Hagan.calibration <- function(tau, f, K, iv, parameter0)
{
  # objective function for optimization
  # variables are transformed because of satisfing the constraint conditions
  
  #objective <- function(x){return(sum( iv - SABR.iv(tau, f, K, exp(x[1]), exp(x[2]), .t2(x[3])) )^2 ) }
  objective <- function(x){return(sum( ( iv - Hagan.IV(tau, f, K, exp(x[2]), exp(x[4]), .t2(x[3]) , exp(x[1])  ) )^2 ) ) }
  x <- optim(c(log(parameter0$nu), log(parameter0$alpha), .t2inv(parameter0$rho),  log(parameter0$beta) ), objective, list(maxit = 1000))
  
  # return the optimized parameters
  parameter <- x$par
  parameter <- c(exp(parameter[1]),exp(parameter[2]),.t2(parameter[3]), exp(parameter[4]) )
  names(parameter) <- c("nu", "alpha", "rho" , "beta")
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

