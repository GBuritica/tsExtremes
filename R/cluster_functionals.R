#######################################################################
#######################################################################
#######################################################################
#######################################################################
### Extremal Index Cluser process function
#' eiCP
#' @description Estimates of the extremal index based on the alpha cluster Process
#'
#' @param path0 (Matrix with the sample trajectory of a time series)
#' @param alpha0 (Integer with an estimate of the tail index)
#' @param klim0 (Integer with the largest order statistic k to consider, or a vector with values k to consider)
#' @param n0 (Integer with the path length)
#' @param plot (F or T if the plot with estimates as a function of k must be shown)
#'
#' @return data frame with estimates of the extremal index with values of k and block lengths b used for inference
#' An estimate of the asymptotic variance is also provided using cluster process inference.
#' @export
#'
#' @examples
#' path  <- ARCHm(20000)
#' alpha <- 1/alphaestimator(path,k1=1200)$xi ## The real value should be one
#' ei    <- eiCP(path,alpha,plot=T)
#' pi    <- piCP(path,alpha,plot=T)

eiCP   <- function(path0,alpha0,klim0=100,n0=length(path0),plot=F){

  path0 <- abs(path0)
  if(length(klim0)==1){
    b            <- unique(floor(sqrt(n0/1:klim0)))   ## identifies unique block lengths
    klim0        <- floor(n0/b^2)                     ## Determines the k sequence
  }else{
    b          <- floor(sqrt(n0/(klim0)))
  }

  eireturn   <- vector(mode="double", length=length(klim0) )
  eivariance <- vector(mode="double", length=length(klim0) )
  for(k in 1:length(klim0)){
    b0 <- b[k]
    ## computes paths
    sumaalpha <- sapply( 1:floor(n0/b0),
                         function(l) sum(path0[((l-1)*b0+1):(l*b0)]^alpha0) )

    maxaalpha <- sapply( 1:floor(n0/b0),
                         function(l) max(path0[((l-1)*b0+1):(l*b0)]^alpha0) )

    ordered   <- order(sumaalpha, decreasing=T) ## returns de order (1), (2), (3)...
    ind       <- ordered[1:klim0[k]]                   ## I chose the ones less than k
    estimate  <- mean( maxaalpha[ind]/sumaalpha[ind] ) ## among these I compute the mean
    variance  <- mean( maxaalpha[ind]^2/sumaalpha[ind]^2)
    ### moves b forward
    eireturn[k]   <- estimate
    eivariance[k] <- variance
  }
  eivariance <- eivariance-(eireturn)^2
  estimate <- data.frame('b'=b,'k'=klim0, 'extremal_index'=eireturn,'ei_variance'=eivariance)
  if(plot){
    ks           <- estimate$k
    plot(ks, estimate$extremal_index,type = "l", xlim = c(2,100) , ylim = c(0,1),
         xlab = "k", ylab = "Extremal index")
    lines(ks, estimate$extremal_index + qnorm(0.975)*sqrt(abs(estimate$ei_variance)/ks) , lty = 3, col = 'black' )
    lines(ks, estimate$extremal_index - qnorm(0.975)*sqrt(abs(estimate$ei_variance)/ks) , lty = 3, col = 'black' )
    points(ks,estimate$extremal_index, pch = 16, cex = 0.5)
    abline(h=estimate$extremal_index[8], lty = 2, col = 'blue')
  }

  return(estimate )
}


#' functionalCP
#' @description Estimates of the cluster statistics and based on the alpha cluster Process
#'
#' @param path0 (Matrix with the sample trajectory of a time series)
#' @param alpha0 (Integer with an estimate of the tail index)
#' @param klim0 (Integer with the largest order statistic k to consider, or a vector with values k to consider)
#' @param n0 (Integer with the path length)
#' @param plot (F or T if the plot with estimates as a function of k must be shown)
#' @param cluster_func (A function to evaluate on extremal blocks. As a default we provide the functional for computing the extremal index)
#'
#' @return data frame with estimates of the the cluster functional estimates values of k and block lengths b used for inference
#' An estimate of the asymptotic variance is also provided using cluster process inference.
#' @export
#'
#' @examples
#' path  <- ARCHm(20000)
#' alpha <- 1/alphaestimator(path,k1=1200)$xi ## The real value should be one
#' two options are available for computing the extremal index
#' ei  <- functionalCP(path,1,100,plot=T)
#' ei2 <- eiCP(path,1,100,plot=T)


functionalCP <- function(path0,alpha0,klim0=100,n0=length(path0),p0='alpha',
                            cluster_func = function(block,p0,alpha0) max(abs(block)^alpha0) , plot=F){

  if(length(klim0)==1){
    b            <- unique(floor(sqrt(n0/1:klim0)))   ## identifies unique block lengths
    klim0        <- floor(n0/b^2)                     ## Determines the k sequence

  }else{
    b          <- floor(sqrt(n0/(klim0)))
  }
  estimate <- NULL
  variance <- NULL

  if(p0=='infty'){
    b0 <- b[k]
    ## computes paths
    maxalpha <- sapply( 1:floor(n0/b0),
                        function(l) max(abs(path0[((l-1)*b0+1):(l*b0)])^alpha0) )

    ordered       <- order(maxalpha, decreasing=T) ## returns de order (1), (2), (3)...
    ind           <- ordered[1:klim0[k]]           ## I chose the ones less than k
    estimate      <- mean(sapply(ind, function(l) cluster_func( path0[((l-1)*b0+1):(l*b0)]/maxalpha[l],p0,alpha0))) ## among these I compute the mean
    variance      <- mean(sapply(ind, function(l) cluster_func( path0[((l-1)*b0+1):(l*b0)]/maxalpha[l],p0,alpha0)))
  }else{
    if(p0=='alpha') p0 = alpha0
    for(k in 1:length(klim0)){
      b0 <- b[k]
      ## computes paths
      suma_p <- sapply( 1:floor(n0/b0),
                           function(l) sum(abs(path0[((l-1)*b0+1):(l*b0)])^p0)^(1/p0))

      ordered       <- order(suma_p, decreasing=T) ## returns de order (1), (2), (3)...
      ind           <- ordered[1:klim0[k]]           ## I chose the ones less than k
      estimate      <- rbind(estimate,mean(sapply(ind, function(l) cluster_func( path0[((l-1)*b0+1):(l*b0)]/suma_p[l],p0,alpha0)))) ## among these I compute the mean
      variance      <- rbind(variance,mean(sapply(ind, function(l) cluster_func( path0[((l-1)*b0+1):(l*b0)]/suma_p[l],p0,alpha0)^2)))

    }
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

  }
  return( estimate )
}
