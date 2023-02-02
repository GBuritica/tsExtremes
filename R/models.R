#' ARCHm
#' @description Function to sample from the model
#' X  = AX + B     with log A ~ N-0.5, B ~ Unif(0,1)
#' The model has index of regular variation alpha = 1
#' and extremal index theta = 0.2792
#'
#' @param n0 (Integer with trajectory length)
#'
#' @return (Vector of length n0 sampled from the SRE model)
#' @export
#'
#' @examples
#' path <- ARCHm(1000)
#' plot.ts(path)
ARCHm <- function(n0){
  x0  <- runif(1,0,1)
  nor <- rnorm(2*n0)
  uni <- runif(2*n0,0,1)
  for(i in 1:(2*n0) ) x0  <- c(x0, ( (exp(nor[i] - 0.5)*(x0[i])) + uni[i] ) )
  return(x0[(n0+1):(2*n0)])
}

#' ARCH2m
#' @description Function to sample from the model
#' X= AX + B     with X = (eta + lambda X) Z^2, Z iid Gaussian.
#' For the default values, this model has index of regular variation alpha = 1
#' and extremal index theta = 0.727
#'
#' @param n0 (Integer with trajectory length)
#' @param eta (Real value)
#' @param lambda (Real value)
#'
#' @return (Vector of length n0 sampled from the SRE model)
#'
#' @export
#'
#' @examples
#' path <- ARCH2m(1000)
#' plot.ts(path)
ARCH2m <- function(n0,eta=2*10^{-5},lambda=.5){
  nor   <- rnorm(2*n0)
  x0    <- eta
  for(i in 1:(2*n0) ) x0  <- c(x0, ( (eta+ (lambda*x0[i]) )*nor[i]^2 ) )
  return(x0[(n0+1):(2*n0)])
}


#' AR model
#' @description Function to sample from the model
#' X = (par0)X + Z.
#' For the default values, this model has index of regular variation alpha identical to Z
#' and extremal index 1-(par0)^alpha
#'
#' @param n0 (Integer with trajectory length)
#' @param par0 (Value in (0,1))
#' @param Z.gen (Function of a trajectory length n to sample from the noise)
#'
#' @return (Vector of length n0 sampled from the SRE model)
#' @export
#'
#' @examples
#' path  <- ARm(1000,0.7)
#' plot.ts(path)
ARm  <- function(n0,par0=0.7, Z.gen = function(n) rt(n, df=1)){
  x0 <- stats::arima.sim(n = n0, list(ar=par0, ma=0), rand.gen= Z.gen )
  return(x0)
}
