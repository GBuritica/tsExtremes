#' Compute (tail)-index
#' @description This function computes Hill (tail)-index estimate of a stationary time series.
#'
#' @param path (vector of univariate observations)
#' @param k1   (vector of order statistics to compute Hill estimates).
#'             It uses k1 = n*0.05 otherwise
#' @param plot (T or F if a plot should be shown).
#'
#' @return A data.frame of hill estimates as a function of k.
#' @export
#'
#' @examples
#' h <- hillestimator(rainfall$LANVEOC[rainfall$SEASON=='SPRING'], 1:200)
#' h <- hillestimator(rainfall$LANVEOC[rainfall$SEASON=='SPRING'], plot=T)
#
hillestimator <- function(path, k1=NULL, plot=F){
  n      <- length(path)
  if(is.null(k1)){
    if(plot) k1 <- 1:floor(n*0.05)
    else     k1 <- floor(n*0.05)
  }                                   ## Initialize k
  sorted <- rev( sort(path, partial = n-rev(k1)+1 ) )    ## Computes order statistics
  hill   <- mean( log(sorted[1:max(1,(k1[1]-1))]/sorted[k1[1]]) )

  if(length(k1) > 1){
    for(k in 2:length(k1)){
      sums  <- ( hill[(k-1)] + log( sorted[ k1[(k-1)] ] ) )*(k1[(k-1)]-1)
      hill  <- c(hill,  ( (sums + sum(log( sorted[ (k1[(k-1)]):(k1[k]-1) ] )) )/(k1[k]-1)  - log(sorted[k1[k]] ))  )     ## Computes Hill
    }
  }

  hill <- data.frame("k" = k1, "Hill.estimate" = hill)

  if(plot){
    plot(hill, type = "l", ylab = "Tail-index", main = "Hill Plot" )
  }

  return(hill)
}


