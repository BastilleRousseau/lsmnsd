
#' Calculate NSD 
#'
#' This function calculates NSD from a time-series of spatial locations 
#' @param x,y Spatial locations
#' @keywords nsd
#' @export
#' @examples
#' data(Christian)
#' ts.plot(NSD_fct(Christian$x, Christian$y))

NSD_fct<-function(x, y) {
	(x-x[1])^2+(y-y[1])^2
}

#' Calculate euclidean distance 
#'
#' This function calculates euclidean distance between spatial locations 
#' @param x,y Spatial locations
#' @keywords distance
#' @export
#' @examples
#' data(Christian)
#' hist(Dist_fct(Christian$x, Christian$y))
Dist_fct<-function(x, y) {
	sqrt((x[2:length(x)]-x[1:(length(x)-1)])^2+(y[2:length(y)]-y[1:(length(y)-1)])^2)
}

#' Mode  
#'
#' This function finds the mode of a vector 
#' @param x A vector
#' @keywords mode
#' @export
#' @examples
#' v<-c(1,2,2,3,4,4,4,4,5,6)
#' Mode(v)
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#' Range standardisation (0,1)  
#'
#' This function standardises a vector between 0 and 1
#' @param x A vector
#' @keywords range standardisation
#' @export
#' @examples
#' v<-c(1,2,2,3,4,4,5,6)
#' range01(v)
range01 <- function(x){(x-min(x))/(max(x)-min(x))}


#' WAIC calculation  
#'
#' This function calculates WAIC from clustNSD output. Based on script initially written by Andrew Gelman
#' @param jagsfit Output of a clustNSD call
#' @keywords waic
#' @export
#' @examples
#' data(Christian_rjags)
#' waic(Christian_rjags)
waic <- function (jagsfit){
  log_lik <- jagsfit$BUGSoutput$sims.list$loglik
  lppd <- sum (log (colMeans(exp(log_lik))))
  p_waic_1 <- 2*sum (log(colMeans(exp(log_lik))) - colMeans(log_lik))
  p_waic_2 <- sum (colVars(log_lik))
  waic_2 <- -2*lppd + 2*p_waic_2
  return(waic_2)
}

#' Posterior variances 
#'
#' This function calculates posterior variances from simulation. Based on script initially written by Andrew Gelman
#' @param a A matrix
#' @keywords waic
#' @export
#' @examples
#' mat<-matrix(rnorm(20), nrow=10, ncol=2)
#. colVars(mat)
colVars <- function (a){
  diff <- a - matrix (colMeans(a), nrow(a), ncol(a), byrow=TRUE)
  vars <- colMeans (diff^2)*nrow(a)/(nrow(a)-1)
  return (vars)
}
