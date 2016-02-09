## lsmnsd  ##

**lsmnsd** is a R package that analyze animal movement using latent-state model and net squared displacement.

This is the development area for the package `lsmnsd`, which provides a series of functions to analyze and classify animal movement strategies based on analysis of Net Squared Displacement (NSD) and latent-state model. 

*References: Bastille-Rousseau, G., Potts, J., Yackulic, C., Frair, J., Ellington, E.H., Blake, S. (In review) Characterizing movement strategies of Galapagos giant tortoises using a Bayesian mixture distribution model and net squared displacement. Movement Ecology.* 

## Installation of the development version  ##

You need to use the package `devtools`from Hadley Wickham. 
    
    library(devtools)
    install_github("BastilleRousseau/lsmnsd")

You also need a version of [JAGS](http://mcmc-jags.sourceforge.net/) and the package `R2jags` installed before using the package.  

## Getting started ##

The package main functions are `clustNSD`, `simple.clust`, and `classify`. The package comes with movement data from two giant tortoises : Christian and Zelfa. For a list of documented functions see the Reference manual. 

Alternatively, here is a quick example to get you going: 

    library(lsmnsd)
    data(Christian) 
    nsd1<-NSD_fct(Christian$x, Christian$y)
    Christian_out<-clustNSD(cbind(range01(nsd1), Christian$Time), n.iter=10000, WAIC=F, simplify=F)
    summary(simple.clust(Christian_out))
    Christian_class<-classify(Christian_out)
    summary(Christian_class)

