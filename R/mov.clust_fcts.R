#' Convert to mov.clust class
#'
#' This function simplifies a rjags output from clustNSD to a mov.clust object
#' @param out Output of clustNSD function
#' @keywords mov.clust simple.clust
#' @export
#' @return A mov.clust object
#' @examples
#' data(Christian_rjags)
#' summary(simple.clust(Christian_rjags))
#' plot(simple.clust(Christian_rjags))
#' data(Zelfa_rjags)
#' summary(simple.clust(Zelfa_rjags))
#' plot(simple.clust(Zelfa_rjags))

simple.clust<-function(out) {
	if(class(out)!="rjags"){
			    print("ERROR: object must be of class rjags")
			} else { 
				n.clust<-length(out$BUGSoutput$mean$mu)-1
				Results<-out$BUGSoutput$summary
				idx <- apply(out$BUGSoutput$sims.list$idx, 2, Mode)
				criteria<-c(out$BUGSoutput$pD,out$BUGSoutput$DIC,out$BUGSoutput$WAIC); names(criteria)<-c("pd", "DIC", "WAIC")	
				param<-c(out$BUGSoutput$n.chains, out$BUGSoutput$n.iter, out$BUGSoutput$n.burnin, out$BUGSoutput$n.thin, out$BUGSoutput$n.sims)	
				NSD<-out$BUGSoutput$data		
				Convergence<-out$BUGSoutput$Convergence
				out<-list(n.clust, Results, idx, criteria, Convergence, param, NSD); names(out)<-c("n.clust", "Results", "idx", "Criteria", "Convergence", "Param_jags", "NSD")
				class(out)<-"mov.clust"
				return(out)
			}
		}


#' Summarizing mov.clust object
#'
#' summary method for class "mov.clust"
#' @param object an object of class "mov.clust"
#' @keywords mov.clust simple.clust
#' @export
#' @examples
#' data(Christian_rjags)
#' summary(simple.clust(Christian_rjags))
summary.mov.clust<-function(object) {
cat("Inference for Bugs model fit using jags", "\n", "n.chains= ", object$Param_jags[1], ", n.iter= ", object$Param_jags[2], ", n.burnin= ", object$Param_jags[3], "\n", "iterations saved = ", object$Param_jags[5], "\n", "\n")
print(object$Results[c("deviance", "mu[1]","mu[2]","mu[3]","sigma[1]","sigma[2]", "sigma[3]", "trans.mat[1,1]","trans.mat[2,2]","trans.mat[3,3]"), c(1:3, 7:9)], digits=3)
cat("\n")
print(object$Criteria, digits=4)
if(object$Convergence==0) {print("Warning: Problem of convergence")}
}     
		
		
		
#' Print mov.clust object
#'
#' printing method for class "mov.clust"
#' @param object an object of class "mov.clust"
#' @keywords mov.clust simple.clust
#' @export
#' @examples
#' data(Zelfa_rjags)
#' simple.clust(Zelfa_rjags)
print.mov.clust <- function(object) {
print(object$Results, digits=3)
print(object$Criteria, digits=3)
}



#' Plot mov.clust object
#'
#' plotting method for class "mov.clust"
#' @param object an object of class "mov.clust"
#' @keywords mov.clust simple.clust
#' @export
#' @examples
#' data(Zelfa_rjags)
#' plot(simple.clust(Zelfa_rjags))
plot.mov.clust<-function(object) {
	if(class(object)=="rjags") {object<-simple.clust(object)}	
	par(mfrow=c(2,1))
	mat<-cbind(object$NSD, object$idx[object$NSD[,2]])
	min1<-min(mat[mat[,3]==1,1])
	max1<-max(mat[mat[,3]==1,1])
	min2<-min(mat[mat[,3]==2,1])
	max2<-max(mat[mat[,3]==2,1])
	prop<-table(object$idx)/length(object$idx)
	n.clust<-object$n.clust
	
	dens<-density(object$NSD[,1])
	plot(dens$x, range01(dens$y), main="Density", type="l", ylab="Density", xlab="Scaled NSD", lwd=4, xlim=c(0,1))
	abline(v=min1, lwd=2, col="Red")
	abline(v=max1, lwd=3, col="Red")
	if(n.clust==2) abline(v=min2, lwd=2, col="Blue")
	if(n.clust==2) abline(v=max2, lwd=2, col="Blue")
		
	plot(mat[,2], mat[,1],type="l", main="Time-series", xlab="Time", ylab="NSD")
	points(mat[mat[,3]==1,2], mat[mat[,3]==1,1], col="Red") 
	abline(h=min1, lwd=2, col="Red")
	abline(h=max1, lwd=3, col="Red")
	
	points(mat[mat[,3]==2,2], mat[mat[,3]==2,1], col="Blue")
	points(mat[mat[,3]==3,2], mat[mat[,3]==3,1], col="Black", pch=17, ps=2)
	abline(h=min2, lwd=2, col="Blue")
	abline(h=max2, lwd=2, col="Blue")	
	legend("topleft", c("Cluster1", "Cluster2", "Uniform"), pch=c(1,1,17), lty=1, lwd=c(1,1,0), col=c("Red", "Blue", "Black"), bty="n", cex=1.25)
		
}

