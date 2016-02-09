#' Rerun clustNSD  
#'
#' This function uses the autojags function to rerun the function clustNSD to facilitate convergence
#' @param rjags Output from clustNSD
#' @param n.update the max number of updates, default =3 
#' @param inc.fact Factor by which the number of iterations will be increased, default=1 (no increase)
#' @keywords rerun
#' @export
#' @examples
#' # DO NOT RUN - TAKES A LONG TIME
#' #data(Christian_rjags)
#' #rerun(Christian_rjags, n.update=2, inc.fact=1)
rerun<-function(rjags, n.update=3, inc.fact=1) {
	n.thin<-rjags$BUGSoutput$n.thin
	n.iter<-rjags$BUGSoutput$n.keep*inc.fact
	out<-autojags(rjags, n.iter=n.iter, n.thin=n.thin, n.update=n.update)
}


#' Diagnostics plots 
#'
#' This function uses the traceplot to display a plot of iterations vs. sampled values for each variable in the chain with a separate plot per variable.
#' @param out Output from clustNSD
#' @keywords t.plot traceplot
#' @export
#' @examples
#' data(Christian_rjags)
#' t.plot(Christian_rjags)
t.plot<-function(out){
		traceplot(out, mfrow=c(4,4), varname=c("mu", "sigma", "deviance", "trans.mat"))
		}

	

