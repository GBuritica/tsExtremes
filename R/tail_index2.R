#######################################################################
#' alphaestimator2
#' @description
#'  This function computes the (tail)-index.
#'  It returns the mean estimate from the Pickands, Hill and DeHaan estimators
#'  at 96th empirical quantile.
#'  It uses the fExtreme code.
#'
#' @param path0 (Vector of univariate observations)
#'
#' @return An integer with the tail - index estimate
#' @export
#'
#' @examples
#' sample   <- abs( arima.sim(n = 8000, list(ar=0.5, ma=0), rand.gen=function(n) rt(n,df=4) ) )
#' alphaestimator2(sample)
alphaestimator2 <- function(path0){
  alpha  <- fExtremes::shaparmPlot(-path0, plottype = "upper", doplot=FALSE) ## Upper is blue
  alpha1 <- 1/alpha$Upper[4,2:4]
  alpha  <- sum(alpha1)/3
  return(alpha)
}
#######################################################################
#' alphaestimator
#' @description
#' This function computes the (tail)-index as a function of k
#' It returns the unbiased Hill estimator from de Haan et al. with tuning parameter rho = 2
#'
#'
#' @param sample (Vector of nonnegative univariate entries)
#' @param k1 (Integer indicating the number of high order statistics to consider for inference)
#' @param plot (T or F indicate if the plot of Hill estimates as a function of k must be shown)
#' @param R0 (Integer with the number of bootstrap replicates to consider for computing confidence intervals)
#' @param hill (T or F indicate if the classical Hill-plot is also plotted)
#' @param ylim0 (Vector with lower and upper bounds of the y-axis for the plot )
#'
#' @return A data.frame with the xi estimate and the confidence intervals lower and upper bounds.
#' @export
#'
#' @examples
#' sample   <- abs( arima.sim(n = 8000, list(ar=0.5, ma=0), rand.gen=function(n) rt(n,df=4) ) )
#' alphaestimator(sample, plot=TRUE , R0 = 100,  hill=TRUE,   k1 = 1000 )
#' alphaestimator(sample, plot=TRUE , R0 = 100,  hill=FALSE,  k1 = 1000 )
#' abline(h=0.25,col = "red")
alphaestimator  <- function(sample,k1=floor(n^(0.7)),plot=FALSE,R0=100,hill=FALSE,ylim0=NULL){
  ### n + transforms to log
  n          <- length(sample)
  lsorted    <- log(sort(abs(sample)))
  #################### Estimating rho parameter
  krhomax   <- floor( min( (sum(!lsorted==-Inf)-1),n, 2*n/log(log(n)) ) ) ## limits for the k(rho)
  rhoes     <- sapply( 2:min(krhomax, max(floor(n^0.7),k1)  ) , function(l) rho_Estimate2( l, lsorted , n  ) )
  rhohat    <- median(rhoes, na.rm = TRUE)
  ################### Estimating gamma
  es        <- gammaes2(lsorted, n,rhohat,k1)
  al        <- (es$hill-es$biais)                   ## Unbiased estimator
  al2       <- es$hill                              ## Hill estimator
  ################## Defyining Bootstrap statistic
  stathill2        <- function(data){
     lpath         <- log(sort(data))
     ind <- c(0.6,0.7,0.8,0.9)
     j   <- 1
     es  <- NULL
     while(k1 > n^ind[j]){
       krhoes     <- sapply( 2:min(krhomax, floor(n^0.7) ) , function(i) rho_Estimate2( i, lpath , n  ) )
       krhohat    <- median(krhoes, na.rm = TRUE)
       esp        <- gammaes2(lpath,n,krhohat, 1:floor(n^0.7) )

       krhoes     <- sapply( 2:min(krhomax, floor(n^ind[j])  ) , function(i) rho_Estimate2( i, lpath , n  ) )
       krhohat    <- median(krhoes, na.rm = TRUE)
       if(j==1)      esp    <- gammaes2(lpath,n,krhohat, 1:floor(n^ind[j]) )
       if(j > 1)     esp    <- rbind(esp,gammaes2(lpath,n,krhohat, (floor(n^ind[j-1])+1):floor(n^ind[j]) ))
       j <- j+1
     }
       krhoes      <- sapply( 2:min(krhomax, k1  ) , function(i) rho_Estimate2( i, lpath , n  ) )
       krhohat    <- median(krhoes, na.rm = TRUE)
       esp        <- rbind(esp,gammaes2(lpath,n,krhohat,  (floor(n^ind[j-1])+1):k1 ))
      return(  c( (esp$hill-esp$biais) , esp$hill ) )
  }           ## Takes longer to run. Also bootstrapping the estimate krho.
  stathill         <- function(data){
    lpath      <- log(sort(data))
    esp        <- gammaes2(lpath,n, rhohat,  1:k1 )
    return(  c( (esp$hill-esp$biais) , esp$hill ) )
  }
  ################## Bootstrap replicates
  if(R0 > 10){
    ## Computing the paths Hill and Unbiased Hill
    ## For different estimates of kro parameter
    ind <- c(0.6,0.7,0.8,0.9)
    j   <- 1
    es  <- NULL
    while(k1 > n^ind[j]){
        rhoesp     <- sapply( 2:min(krhomax, floor(n^ind[j])  ) , function(l) rho_Estimate2( l, lsorted , n  ) )
        rhohatp    <- median(rhoes, na.rm = TRUE)
        if(j==1)      es    <- gammaes2(lsorted, n,rhohatp, 1:floor(n^ind[j]) )
        if(j > 1)     es    <- rbind(es,gammaes2(lsorted, n,rhohatp, (floor(n^ind[j-1])+1):floor(n^ind[j]) ))
        j <- j+1
    }
    es   <- rbind(es,gammaes2(lsorted, n,rhohat, (floor(n^ind[j-1])+1):k1 ))
    ##################
    ################## Bootstrap
    b          <-  boot::tsboot(sample,statistic=stathill, R=R0, sim = "geom", l = 200 )
    IC1        <-  sapply(1:length(es$hill) , function(l)  boot::boot.ci(b,  type = "perc", index = l )$percent[4:5] )
    IC2        <-  sapply(1:length(es$hill) , function(l)  boot::boot.ci(b,  type = "perc", index = (length(es$hill) + l)  )$percent[4:5] )
    ################## Plots
    if(plot == TRUE){
      if(length(ylim0)==0) ylim0 <- c(min(IC1, IC2),max((IC1),IC2))
      if(hill){
          plot.ts((es$hill-es$biais), ylim = ylim0 , main = "Hill plot" ,xlab = "k", ylab = " ",col = "darkblue")
          lines( es$hill, col = "grey")
          for(i in 1:2) lines(IC2[i,],lty=3, col = "grey")
          for(i in 1:2) lines(IC1[i,],lty=2, col = "darkblue")

          #polygon( x=c(1:length(es$hill), rev(1:length(es$hill) ) ), y=c(IC2[2,],rev(IC2[1,]) ), col = "grey" , lty=1, density=20, angle = 40)
          #polygon( x=c(1:length(es$hill), rev(1:length(es$hill) ) ), y=c(IC1[2,],rev(IC1[1,]) ), col = "skyblue" , lty=1, density=20, angle = 90)
          #sapply(1:2, function(k) lines( IC2[k,], lty = 1, col = "grey"))

          #lines( (es$hill), col = "grey", lty = 2)
          #lines( (es$hill-es$biais), col = "darkblue", lty=1)
          legend("topright", legend = c("Unbiased Hill", "Hill"), col =c("darkblue","grey"), lty=1 ,cex=0.7)
      }
      else{
        plot.ts(NA, xlim = c(10,k1), ylim = ylim0 , main = "Hill plot" ,xlab = "k", ylab = " ",col = "darkblue")
        #polygon( x=c(1:length(es$hill), rev(1:length(es$hill) ) ), y=c(IC1[2,],rev(IC1[1,]) ), col = "grey" , lty=1, density=20, angle = 90)
        for(k in 1:2 ) lines( IC1[k,], lty = 2, col = "darkblue")
        lines((es$hill-es$biais), col = "darkblue" )
        #abline(h=al); print(al)
      }
    }
      return(data.frame("lbound"=IC1[1,k1] , "xi" = al, "ubound"=IC1[2,k1]))
  }

  return( data.frame("xi"=al))
}

#######################################################################
## Auxiliar fct. for unbiased Hill estimator in De Haan - Mercadier - Zhou
##
Mka   <- function(k,a,lsortedpath,n) return( mean( ( lsortedpath[(n-k+1):n] - lsortedpath[(n-k)])^a ))
Rka   <- function(k,v,lsortedpath,n){
  mk2 <- Mka(k,2,lsortedpath,n)
  mk1 <- Mka(k,1,lsortedpath,n)
  na <- length(v)
  res <- 1:na
  for( i in 1:na){
    mka    <- Mka(k,v[i],lsortedpath,n)
    res[i] <- ( mka - ( gamma( (v[i]+1) )*((mk1)^(v[i]) ) ) )/( mk2 - ( 2*((mk1)^2) ) )
  }
  names(res) <- v
  return(res)
}
Ska   <- function(k,a,lsortedpath,n){
  na <- length(a); v <- NULL
  for( i in 1:na) v <- c(v, c( (2*a[i]) ,(a[i]+1) ) )
  rk <- Rka(k,v,lsortedpath,n)
  sk <- sapply(1:na ,
               function(i) ( a[i]*( (a[i]+1)^2 )*( (gamma(a[i]))^2 )*rk[(2*(i-1) +1)] )/( 4*gamma( (2*a[i]) )*( rk[(2*i)]^2 ) )   )
  names(sk) <- a
  return(sk)
}
sarho <- function(p,a){
  u <- sapply(a, function(al) (p^2)*( 1- ((1-p)^(2*al)) - ( 2*al*p*( (1-p)^( (2*al)-1 ) ) ) ) )
  l <- sapply(a, function(al)  ( 1 - ( (1-p)^(al+1) ) - ( (al+1)*p*( (1-p)^(al) ) ) )^2 )
  res <- u/l; names(res) <- a
  return(res)
}  ## ok decroit
rho_Estimate  <- function(k,lsortedpath,n){
  s2  <- Ska(k,2,lsortedpath,n)
  rho <- NA
  if(s2 > 2/3 && s2 < 3/4){
    rho <- ( ( 2*( (3*s2) - 2 ) ) + sqrt( (3*s2) - 2 ) )/( 3 - ( 4*s2 ) )
  }
  return( -1*rho )
}
#### Following the results for the simulation study of de Haan we take alpha = 2 Then,
rho_Estimate2 <- function(k,lsortedpath,n){
  mk4 <- Mka(k,a=4,lsortedpath,n)
  mk1 <- Mka(k,a=1,lsortedpath,n)
  mk2 <- Mka(k,a=2,lsortedpath,n)
  mk3 <- Mka(k,a=3,lsortedpath,n)
  if( (mk3 - ( 6*( (mk1)^3 ) ))^2  > 0 ){
    s2  <- ( (0.75)*( mk4 - ( (24)*( (mk1)^4 ) ) )*( mk2 - ( 2*( (mk1)^2 ) ) ) )
    s2 <- s2/( mk3 - ( 6*( (mk1)^3 ) ))^2
    rho <- NA
    if(s2 > 2/3 && s2 < 3/4) rho <- ( -4 + (6*s2) + sqrt( ((3*s2)-2) ) )/( (4*s2)-3 )
    return( rho )
  }
  else return(NA)
}
gammaes2      <- function(lsorted,n,rho,k0=1:(n^0.9)){
  if(is.na(rho)) rho <- 1         ## returns Hill estimator
  hill  <-    sapply( k0, function(k) Mka(k,1,lsorted,n))
  biais <-    sapply( 1:length(k0), function(k) ( ( Mka(k0[k],2,lsorted,n) - 2*( (hill[k])^2 ))*(1-rho)*(1/( 2*hill[k]*rho )) ))
  return( data.frame("hill" = hill, "biais" = biais)  )
}
#plot(seq(-2,-0.01,0.01),sapply(seq(-2,-0.01,0.01), function(p) sarho(p,1.5) )) # function sarho
#######################################################################
#######################################################################

