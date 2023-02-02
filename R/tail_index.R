#' Compute (tail)-index
#' @description This function computes Hill (tail)-index estimate of a stationary time series.
#'
#' @param path (vector of univariate observations)
#' @param k0   (vector of order statistics to compute Hill estimates).
#'             It uses k0 = n*0.05 otherwise
#' @param plot (T or F if a plot should be shown).
#'
#' @return A data.frame of hill estimates as a function of k.
#' @export
#'
#' @example
#' hillestimate(rainfall$LANVEOC[rainfall$SEASON=='SPRING'], 1:200)
#' hillestimate(rainfall$LANVEOC[rainfall$SEASON=='SPRING'], plot=T)

hillestimate <- function(path, k0=NULL, plot=F){
  n      <- length(path)
  if(is.null(k0)){
    if(plot) k0 <- 1:floor(n*0.05)
    else     k0 <- floor(n*0.05)
  }                                   ## Initialize k
  sorted <- rev( sort(path, partial = n-rev(k0)+1 ) )    ## Computes order statistics
  hill   <- mean( log(sorted[1:max(1,(k0[1]-1))]/sorted[k0[1]]) )

  if(length(k0) > 1){
    for(k in 2:length(k0)){
      sums  <- ( hill[(k-1)] + log( sorted[ k0[(k-1)] ] ) )*(k0[(k-1)]-1)
      hill  <- c(hill,  ( (sums + sum(log( sorted[ (k0[(k-1)]):(k0[k]-1) ] )) )/(k0[k]-1)  - log(sorted[k0[k]] ))  )     ## Computes Hill
    }
  }

  hill <- data.frame("k" = k0, "Hill.estimate" = hill)

  if(plot){
    plot(hill, type = "l", ylab = "Tail-index", main = "Hill Plot" )
  }

  return(hill)
}


