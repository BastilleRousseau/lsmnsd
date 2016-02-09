#' Classify movement strategies 
#'
#' Extract information to classify movement strategies
#' @param object an object of class "mov.clust" or "rjags"
#' @param grph Produce a graph of classification (Default = F) 
#' @keywords classify clust.classify mov.clust
#' @export
#' @return A clust.classify object
#' @examples
#' data(Christian_rjags)
#' summary(classify(Christian_rjags))	
#' data(Zelfa_rjags)
#' summary(classify(Zelfa_rjags))
classify<-function(object, grph=F,...){
	if(class(object)=="rjags") {object<-simple.clust(object)}	
	mat<-cbind(	1:nrow(object$NSD), object$NSD, object$idx[object$NSD[,2]])
	mean<-tapply(mat[,2], mat[,3],mean)
	sd<-tapply(mat[,2], mat[,3],sd)
	predict<-rbind(object$Results[c("mu[1]","sigma[1]","mu[2]","sigma[2]"), 1], c(mean[1], sd[1], mean[2], sd[2]))
	rownames(predict)<-c("MCMC", "Observed")
	TS<-time.spent(object)
	brk<-breaks(object)
	overlap<-integrate(min.f1f2, -Inf, Inf, mu1=predict[1,1], mu2=predict[1,3], sd1=predict[1,2], sd2=predict[1,4])[[1]]
	switch<-switch.matrix(object)
	Strategy<- ifelse(switch[[1]][2,2]<0.90 & switch[[1]][3,3]<0.90, "Resident",
				ifelse(switch[[1]][1,1]>0.95 & switch[[1]][2,2]>0.95 & switch[[1]][3,3]>0.85, "Migration",
				ifelse(switch[[1]][1,1]>0.95 & switch[[1]][2,2]>0.90 & switch[[1]][3,3]<0.85, "Nomadism- recommend visual inspection", "Inspection required"))) 
	Strategy<-ifelse(as.numeric(brk[2])==0 & Strategy=="Migration", "Dispersal", Strategy)				
	out<-list(predict, TS, brk, overlap, switch, Strategy)
	names(out)<-c("Parameters", "Time spent", "Transitions", "Weitzman's overlap", "Switching probabilities", "Movement strategy")
	class(out)<-"clust.classify"
	if (grph==T) {plot(object)}
	return(out)
	}
	

#' Summarizing clust.classify object
#'
#' summary method for class "clust.classify"
#' @param object an object of class "clust.classify"
#' @keywords classify clust.classify
#' @export
#' @examples
#' data(Christian_rjags)
#' summary(classify(Christian_rjags))	
summary.clust.classify <- function(object) {
cat("Summary information for clustering applied to NSD time-series", "\n", "\n", "Estimates of mean and standard deviation for each cluster", "\n")
print(object[[1]], digits=3)
cat("\n", "Overlap between the two normal distributions:")
print(object[[4]], digits=3)
cat("\n", "Total time spent in each cluster (3 is associated to the uniform distribution):")
print(object[[2]][[1]], digits=0)
cat("\n", "Average time spent in each cluster (3 is associated to the uniform distribution):", "\n")
print(object[[2]][[2]], digits=3)
cat("\n", "Number of transition from first to second cluster: ", object[[3]][[1]],  "\n", "Number of transition from second to first cluster: ", object[[3]][[2]], "\n","Number of departure from first cluster: ", object[[3]][[3]], "\n")
cat("\n", "Matrix of switching probability:", "\n")
print(object[[5]][[1]], digits=3)
cat("\n", "Movement strategy:", "\n")
print(object[[6]][1], digits=3)
}	
	

#' Switching probabilties matrix
#'
#' Extract matrix of switching probabilities
#' @param object an object of class "mov.clust" or "rjags"
#' @keywords classify clust.classify
#' @export
#' @examples
#' data(Zelfa_rjags)
#' switch.matrix(Zelfa_rjags)	
switch.matrix<-function(object){
tt<-object$Results[c("trans.mat[1,1]","trans.mat[1,2]","trans.mat[1,3]","trans.mat[2,1]","trans.mat[2,2]","trans.mat[2,3]","trans.mat[3,1]","trans.mat[3,2]","trans.mat[3,3]"), 1:2]	
mat<-matrix(tt[,1], 3,3, byrow=T)
rownames(mat)<-c("From N1", "From N2", "From Uni")
colnames(mat)<-c("To N1", "To N2", "To Uni") 
mat2<-matrix(tt[,2], 3,3, byrow=T)
rownames(mat2)<-c("From N1", "From N2", "From Uni")
colnames(mat2)<-c("To N1", "To N2", "To Uni") 
out<-list(mat, mat2); names(out)<-c("Mean", "SD")
return(out)
}



# For overlap
min.f1f2 <- function(x, mu1, mu2, sd1, sd2) {
    f1 <- dnorm(x, mean=mu1, sd=sd1)
    f2 <- dnorm(x, mean=mu2, sd=sd2)
    pmin(f1, f2)
}



#' Calculating transitions among clusters
#'
#' Calculating pattern of transitions among clusters
#' @param object an object of class "mov.clust" or "rjags"
#' @keywords classify clust.classify
#' @export
#' @examples
#' data(Zelfa_rjags)
#' breaks(Zelfa_rjags)		
breaks<-function(object) {
	if(class(object)=="rjags") {object<-simple.clust(object)}	
	if (object$n.clust==1) {print("WARNING: Function is designed for two clusters")}
	AB<-sum(diff(object$idx[object$idx!=3])==1)
	BA<-sum(diff(object$idx[object$idx!=3])==-1)
	idx2<-ifelse(object$idx==3, 2, object$idx)
	AC<-sum(diff(idx2)==1)
	bb<-c(AB, BA, AC)
	names(bb)<-c("Nb_1st_2nd", "Nb_2nd_1st", "Nb_departure_1st")
	return(bb)
}
	

#' Calculating time-spent in cluster
#'
#' Calculating pattern of time spent in each cluster
#' @param object an object of class "mov.clust" or "rjags"
#' @keywords classify clust.classify
#' @export
#' @examples
#' data(Zelfa_rjags)
#' time.spent(Zelfa_rjags)
time.spent<-function(object){	
	if(class(object)=="rjags") {object<-simple.clust(object)}	
	if (object$n.clust==1) {print("WARNING: Function is designed for two clusters")}
	TS_total<-table(object$idx)
	last<-tail(object$idx,1)
	ff1<-cumsum(ifelse(object$idx==1, 1, 0))
	ff2<-cumsum(ifelse(object$idx==2, 1, 0))
	ff3<-cumsum(ifelse(object$idx==3, 1, 0))
	tt1<-unique(ff1[duplicated(ff1)])
	TS1<-c(tt1[1], tt1[-1]-tt1[-length(tt1)])
	if (last==1) {TS1<-c(TS1, TS_total[1]-sum(TS1))}
	tt2<-unique(ff2[duplicated(ff2)])
	TS2<-tt2[-1]-tt2[-length(tt2)]
	if (last==2) {TS2<-c(TS2, TS_total[2]-sum(TS2))}
	tt3<-unique(ff3[duplicated(ff3)])
	TS3<-tt3[-1]-tt3[-length(tt3)]
	if (last==3) {TS3<-c(TS3, TS_total[3]-sum(TS3))}
	Average<-c(mean(TS1), mean(TS2), mean(TS3))
	names(Average)<-1:3
	out<-list(TS_total, Average, TS1, TS2, TS3)
	names(out)<-c("Total_time", "Average", "Time_in_1", "Time_in_2","Time_in_3")
	return(out) 
}



