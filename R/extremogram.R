




#' Row maxima
#'
#' @description
#' `rmax()` computes the maxima row-wise
#'
#' @param mat Matrix
#'
#' @return A numerical vector with the same number of columns as `mat`
#'
#' @export
#'
#' @examples
#' x <- rbind( c(1,2), c(1,1) )
#' rmax(x)
rmax <- function(mat){
  apply(mat, 1, function(x) max(x) )
}


#' Empirical temporal extremogram
#'
#' @param sample A numerical vector to compute the extremogram
#' @param maxlag An integer with the
#' @param q A quantile level in (0,1).
#' @param plot A logical value to plot the results
#'
#' @return A vector with the extremogram values
#' @export
#'
#' @examples
#' pre    <- tsExtremes::rainfall ## load rainfall data
#' sample <- na.omit(sample)      ## removes columns with missing values
#' sample <- tsExtremes::rmax(pre[pre$SEASON=='FALL' , 5:7]) ## computes vector-wise maxima
#' tsExtremes::extremogram(sample)
#'
#' @details
#' The temporal extremogram \eqn{\chi_t} is defined by
#' \deqn{\chi_t = \lim_{t \to \infty} \mathbb{P}( X_t > x | X_0 > x).}
#' Its empirical version, computes the average number of exceedances of the q-th order statistic.
#' As a baseline, the extremogram takes the value of 1-q at independent lags.

extremogram <- function(sample,maxlag=45,q=.95,plot=T){
  sorted <- sort(sample,decreasing=T)
  quant  <- sorted[floor(n*(1-q))]
  ext    <- function(sample, maxlag=10, q) sapply( 1:maxlag, function(k) mean( sample[ (which(sample > q) + k) ] > q , na.rm = T) )
  extremo <- ext(sample, maxlag, q=quant)
  if(plot){
    plot(x=NULL,y=NULL, xlim = c(1,maxlag),
         ylim = c(0,max(0.6,max(extremo)) ), xlab = "time lag",
         ylab = "Extremogram" , main = '')
  }
  for(k in 1:(maxlag+1)) segments( k,0,k,extremo[k])
  abline(h=1-q,lty=2)
  ext
}

