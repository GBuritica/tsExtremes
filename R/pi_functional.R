#' piCP
#' @description Estimates of the cluster legths and the extrema index based on the alpha cluster Process
#'
#' @param path0 (Matrix with the sample trajectory of a time series)
#' @param alpha0 (Integer with an estimate of the tail index)
#' @param klim0 (Integer with the largest order statistic k to consider, or a vector with values k to consider)
#' @param n0 (Integer with the path length)
#' @param plot (F or T if the plot with estimates as a function of k must be shown)
#'
#' @return data frame with estimates of the the cluster lengths pi_j with values of k and block lengths b used for inference
#' An estimate of the asymptotic variance is also provided using cluster process inference.
#' @export
#'
#' @examples
#' path  <- ARCHm(20000)
#' alpha <- 1/alphaestimator(path,k1=1200)$xi ## The real value should be one
#' pi    <- piCP(path,alpha,plot=T)
piCP    <- function(path0,alpha0,klim0=100,n0=length(path0),plot=F){
  path0 <- abs(path0)
  if(length(klim0)==1){
    b            <- unique(floor(sqrt(n0/1:klim0)))   ## identifies unique block lengths
    klim0        <- floor(n0/b^2)                     ## Determines the k sequence
  }else{
    b          <- floor(sqrt(n0/(klim0)))
  }

  estimate <- NULL
  variance <- NULL
  for(k in 1:length(klim0)){
    b0 <- b[k]
    ## computes paths
    sumaalpha <- sapply( 1:floor(n0/b0),
                         function(l) sum(path0[((l-1)*b0+1):(l*b0)]^alpha0) )

    maxaalpha <- sapply( 1:floor(n0/b0), function(l)
      sort(path0[((l-1)*b0+1):(l*b0)]^alpha0, decreasing = T )[1:5])

    ordered       <- order(sumaalpha, decreasing=T) ## returns de order (1), (2), (3)...
    ind           <- ordered[1:klim0[k]]           ## I chose the ones less than k
    estimate      <- rbind(estimate,  c( mean( maxaalpha[1,ind]/sumaalpha[ind]),
                                         sapply(1:4, function(k)
                                           mean((maxaalpha[k,ind]-maxaalpha[(k+1),ind])/sumaalpha[ind] ) ) ) )  ## among these I compute the mean
    colnames(estimate)  <- c('Theta', 'Pi1', 'Pi2', 'Pi3', 'Pi4')
    estimate        <- as.data.frame(estimate)
    variance      <- rbind(variance, c( mean( maxaalpha[1,ind]^2/sumaalpha[ind]^2),
                                        sapply(1:4, function(k)
                                          mean( (maxaalpha[k,ind]-maxaalpha[(k+1),ind])^2/sumaalpha[ind]^2 )  )))
    colnames(variance)  <- c('Theta', 'Pi1', 'Pi2', 'Pi3', 'Pi4')
    variance          <- as.data.frame(variance)
  }
  variance <- variance-(estimate)^2
  estimate <- list('estimate' = estimate, 'variance' = variance, 'k' = klim0, 'b' = b)
  if(plot){
    par(mfrow = c(2,1))
    ks           <- estimate$k
    plot(ks, estimate$estimate[,1],type = "l", xlim = c(0,100) , ylim = c(0,1),
         xlab = "k", ylab = "Extremal index")
    lines(ks, estimate$estimate[,1] + qnorm(0.975)*sqrt(abs(estimate$variance[,1])/ks) , lty = 3, col = 'black' )
    lines(ks, estimate$estimate[,1] - qnorm(0.975)*sqrt(abs(estimate$variance[,1])/ks) , lty = 3, col = 'black' )
    points(ks,estimate$estimate[,1], pch = 16, cex = 0.5)
    abline(h=estimate$estimate[8,1], lty = 2, col = 'blue')

    plot(ks, estimate$estimate[,2],type = "l", xlim = c(0,100) , ylim = c(0,1),
         xlab = "k", ylab = "Pi1")
    lines(ks, estimate$estimate[,2] + qnorm(0.975)*sqrt(abs(estimate$variance[,2])/ks) , lty = 3, col = 'black' )
    lines(ks, estimate$estimate[,2] - qnorm(0.975)*sqrt(abs(estimate$variance[,2])/ks) , lty = 3, col = 'black' )
    points(ks,estimate$estimate[,2], pch = 16, cex = 0.5)
    abline(h=estimate$estimate[8,2], lty = 2, col = 'blue')
    par(mfrow = c(1,1))
  }
  return( estimate )
}

