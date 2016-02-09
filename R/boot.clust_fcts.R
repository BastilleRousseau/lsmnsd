#' Iterative latent-state model of NSD
#'
#' Bootstrap latent-state model by changing starting locations 
#' @param x, y, time,  information on locations and time (continuous)
#' @param interval Amount of variation in time to be tested
#' @param by Time step between each iteration (number of iterations is a function of the interval and time step)
#' @param n.iter number of total iterations per chain (including burn in; default: 5000)
#' @keywords clustNSD boot.clust bootNSD
#' @return A boot.clust object
#' @export
#' @examples
#' data(Christian)
#' #May results in low convergence given low number of iterations
#' boot1<-bootNSD(Christian$x[1:300], Christian$y[1:300], Christian$Time[1:300], interval=20, by=5, n.iter=10000) 
#' summary(boot1)
#' plot(boot1)
bootNSD<-function(x,y, day, interval=30, by=5, n.iter=20000, ...) {
	out<-list()
	out2<-data.frame(day=day, x=x, y=y)
	start<-seq(1, interval, by)
	end<-length(x)-rev(start)+1
	out2[,seq(4,length(start)*2)]<-NA
	for (i in 1:length(start)) {
		id<-seq(start[i], end[i], 1)
		ttt<-cbind(range01(NSD_fct(x[id], y[id])), day[id]-min(day[id])+1)  
		out[[i]]<-clustNSD(ttt, n.clust=2, n.iter=n.iter, WAIC=F, simplify=T)
		out2[id,2*i+2]<-out[[i]]$NSD[,1]
		out2[id,2*i+3]<-out[[i]]$idx[ttt[,2]] ### Problematic with missing data		
}
out2$prct1<-rowSums(out2[,seq(5,ncol(out2), 2)]==1, na.rm=T)/rowSums(out2[,seq(5,ncol(out2), 2)]>0, na.rm=T)
out2$prct2<-rowSums(out2[,seq(5,ncol(out2), 2)]==2, na.rm=T)/rowSums(out2[,seq(5,ncol(out2), 2)]>0, na.rm=T)
out2$prct3<-rowSums(out2[,seq(5,ncol(out2), 2)]==3, na.rm=T)/rowSums(out2[,seq(5,ncol(out2), 2)]>0, na.rm=T)
out3<-list(out2, out)
class(out3)<-"boot.clust"
return(out3) 
}



#' Plotting boot.clust object
#'
#' plot method for class "boot.clust"
#' @param object an object of class "boot.clust"
#' @keywords bootNSD boot.clust
#' @export
#' @examples
#' #plot(boot1)
plot.boot.clust<-function(object) {
	n.iter<-(ncol(object[[1]])-3)/2
	plot(object[[1]]$day, seq(0,1, length.out=nrow(object[[1]])), type="n", xlab="Time", ylab="NSD", main="Bootstrap with different starting date")
	lines(object[[1]]$day, apply(object[[1]][,c("prct1", "prct2", "prct3")], 1, max, na.rm=T), col="Green", lwd=2, lty=2)	
	for (i in 1:n.iter) {
		lines(object[[1]]$day, object[[1]][,2*i+2], col="grey")
	}
	}

#' Summarizing boot.clust object
#'
#' summary method for class "boot.clust"
#' @param object an object of class "boot.clust"
#' @keywords bootNSD boot.clust
#' @export
#' @examples
#' summary(boot1)
summary.boot.clust<-function(object) {
cat("Summary information of bootstrapping for clustering applied to NSD time-series", "\n", "\n","%Days with 100% agreement")
print(sum(rowSums(object[[1]][,c("prct1", "prct2", "prct3")]==1, na.rm=T))/nrow(object[[1]]) , digits=3)
cat("%Days with 95% agreement")
print(sum(rowSums(object[[1]][,c("prct1", "prct2", "prct3")]>=0.95, na.rm=T))/nrow(object[[1]]) , digits=3)
cat("%Days with 80% agreement")
print(sum(rowSums(object[[1]][,c("prct1", "prct2", "prct3")]>=0.80, na.rm=T))/nrow(object[[1]]) , digits=3)
cat("\n", "Average agreement")
print(mean(apply(object[[1]][,c("prct1", "prct2", "prct3")], 1, max, na.rm=T), na.rm=T), digits=3)
cat("\n", "% of iterations converged")
sum(unlist(lapply(object[[2]], function(x) x$Convergence)))/length(object[[2]])
}
