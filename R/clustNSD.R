#' Latent-state model of NSD
#'
#' Perform latent-state model to characterize movement patterns based on NSD
#' @param data a maxtrix with x,y, and time columns
#' @param WAIC save log-likelihood of every iteration to allow calculations of WAIC, default=FALSE
#' @param n.iter number of total iterations per chain (including burn in; default: 5000)
#' @param n.chains number of Markov chains (default: 3)
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the beginning. Default is n.iter/2, that is, discarding the first half of the simulations. If n.burnin is 0, jags() will run 100 iterations for adaption.
#' @param n.thin thinning rate. Must be a positive integer. Set n.thin > 1 to save memory and computation time if n.iter is large. Default is max(1, floor(n.chains * (n.iter-n.burnin) / 1000)) which will only thin if there are at least 2000 simulations.
#' @param simplify Convert output to mov.clust object. Default=FALSE. See simple.clust for details
#' @param sigma1.min Lower limit of uniform prior for SD of first normal distribution (Default=0.001)
#' @param sigma1.max Upper limit of uniform prior for SD of first normal distribution (Default=0.1)
#' @param sigma2.min Lower limit of uniform prior for SD of second normal distribution (Default=0.001)
#' @param sigma2.max Upper limit of uniform prior for SD of second normal distribution (Default=0.1)
#' @param mu1.min Lower limit of uniform prior for mean of first normal distribution (Default=0.001)
#' @param mu1.max Upper limit of uniform prior for mean of first normal distribution (Default=0.5)
#' @param mu2.min Lower limit of uniform prior for difference between mean of first and second normal distribution (Default=0)
#' @param mu2.max Upper limit of uniform prior for difference between mean of first and second normal distribution (Default=1)
#' @keywords mov.clust clustNSD
#' @return A rjags or mov.clust object
#' @export
#' @examples
#' data(Christian)
#' nsd1<-NSD_fct(Christian$x, Christian$y)
#' Christian_rjags<-clustNSD(cbind(range01(nsd1), Christian$Time), n.iter=10000, WAIC=T, simplify=F)
#' summary(simple.clust(Christian_rjags))
#' data(Zelfa)
#' nsd2<-NSD_fct(Zelfa$x, Zelfa$y)
#' Zelfa_rjags<-clustNSD(cbind(range01(nsd2), Zelfa$Time), n.iter=10000, WAIC=F, simplify=F)
#' summary(simple.clust(Zelfa_rjags))



clustNSD<-function(data, WAIC=FALSE,  n.iter=5000, n.chains=3, n.burnin=floor(n.iter/2), n.thin=max(1, floor((n.iter - n.burnin) / 1000)), simplify=FALSE,sigma1.max=0.1,sigma2.max=0.1,sigma1.min=0.001,sigma2.min=0.001,mu1.max=0.5,mu2.max=1,mu1.min=0,mu2.min=0,  ...){
				if(require("R2jags")){
			    print("R2jags is loaded correctly")
			} else {
			    print("trying to install R2jags")
			    install.packages("R2jags")
			    if(require(R2jags)){
			        print("R2jags installed and loaded")
			    } else {
			        stop("could not install R2jags")
			    }
			}
data[,1]<-range01(data[,1])				

	 if(WAIC==F){out<-jags(list(NSD=data[,1],day=data[,2], ndays=max(data[,2]), npts=nrow(data), sigma1.max=sigma1.max,sigma2.max=sigma2.max,sigma1.min=sigma1.min,sigma2.min=sigma2.min,mu1.max=mu1.max,mu2.max=mu2.max,mu1.min=mu1.min,mu2.min=mu2.min), NULL, c('mu', 'sigma',  "idx","trans.mat"), textConnection(M1), n.chains, n.iter)
				close(out$model.file); out$BUGSoutput$WAIC<-NA; out$BUGSoutput$data<-data
				out$BUGSoutput$Convergence<-ifelse(max(out$BUGSoutput$summary[c("deviance", "mu[1]", "mu[2]", "sigma[1]", "sigma[2]","trans.mat[1,2]","trans.mat[1,3]","trans.mat[2,1]","trans.mat[2,2]","trans.mat[2,3]","trans.mat[3,1]","trans.mat[3,2]","trans.mat[3,3]"), 8])<1.1, 1,0) 
				
				}
				
	if(WAIC==T){out<-jags(list(NSD=data[,1],day=data[,2], ndays=max(data[,2]), npts=nrow(data), sigma1.max=sigma1.max,sigma2.max=sigma2.max,sigma1.min=sigma1.min,sigma2.min=sigma2.min,mu1.max=mu1.max,mu2.max=mu2.max,mu1.min=mu1.min,mu2.min=mu2.min), NULL, c('mu', 'sigma',  "idx", "loglik","trans.mat"), textConnection(M1), n.chains, n.iter)
				close(out$model.file); 	out$BUGSoutput$WAIC<-waic(out); out$BUGSoutput$data<-data
				out$BUGSoutput$Convergence<-ifelse(max(out$BUGSoutput$summary[c("deviance", "mu[1]", "mu[2]", "sigma[1]", "sigma[2]", "trans.mat[1,1]","trans.mat[1,2]","trans.mat[1,3]","trans.mat[2,1]","trans.mat[2,2]","trans.mat[2,3]","trans.mat[3,1]","trans.mat[3,2]","trans.mat[3,3]"), 8])<1.1, 1,0) 
				}			


if (out$BUGSoutput$Convergence==0){print("PROBLEM WITH CONVERGENCE")}
if(simplify) {out<-simple.clust(out)}
	return(out)
}




# Bi-modal
M1 <- 'model  {
    mu[1] ~ dunif(mu1.min,mu1.max) 
    mu[2] <- min(mu[1] + eps[1], 1)
    mu[3] <- 0.5
    eps[1]~dunif(mu2.min,mu2.max)
    sigma[1] ~ dunif(sigma1.min,sigma1.max)  
    sigma[2] ~ dunif(sigma2.min,sigma2.max) 
    sigma[3]<-10    
    tau[1] <- pow(sigma[1],-2)
    tau[2] <- pow(sigma[2],-2)
    tau[3] <- pow(sigma[3],-2) 
    alpha1[1] <- 1
  	alpha1[2] <- 1
  	alpha1[3] <- 1
 	p.state[1:3] ~ ddirch(alpha1[])
    idx[1] ~ dcat(p.state[])	
    alpha[1]<-1
	alpha[2]<-1
	alpha[3]<-1
	trans.mat[1,1:3]~ddirch(alpha[])
	trans.mat[2,1:3]~ddirch(alpha[])
	trans.mat[3,1:3]~ddirch(alpha[])
	pi<-3.141593
		
	for (t in 2:ndays){
		idx[t]~dcat(trans.mat[idx[t-1],])
       }
    for (i in 1:npts) {
 	    loglik[i] <- -0.5*log(2*pi*sigma[idx[day[i]]]) - (NSD[i]-mu[idx[day[i]]])/(2*sigma[idx[day[i]]]) # Necessary for the WAIC
 	    NSD[i]~dnorm(mu[idx[day[i]]],tau[idx[day[i]]])
 	    }
}'
