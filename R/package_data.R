#' lsmnsd: A package for classifying animal movement strategies based on a latent-space model and net squared displacement.
#'
#' The lsmnsd package provides four categories of important functions:
#' clustNSD, simple.clust, classify, bootNSD.
#' @docType package
#' @name lsmnsd
NULL


#' Christian - daily spatial locations 
#'
#' A dataset containing the x, y locations (UTM), and time of Christian, a male giant tortoise inhabiting Alcedo volcano on Isabela Island, Galapagos. 
#'
#' @format A data frame with 1536 rows and 3 variables (x, y, time)
"Christian"


#' Zelfa - daily spatial locations 
#'
#' A dataset containing the x, y locations (UTM), and time of Zelfa, a female giant tortoise inhabiting Espanola Island, Galapagos. 
#'
#' @format A data frame with 1537 rows and 3 variables (x, y, time)
"Zelfa"


#' Latent-state model applied to Christian  
#'
#' The output of clustNSD function applied to Christian. 
#'
#' @format A "rjags" object 
"Christian_rjags"

#' Latent-state model applied to Zelfa  
#'
#' The output of clustNSD function applied to Zelfa. 
#'
#' @format A "rjags" object 
"Zelfa_rjags"