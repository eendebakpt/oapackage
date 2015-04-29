# Testing functions
#foo <- function(x){x*2}

#cfoo=function(a,b){.C('cfoo',as.double(a),as.double(b),c=as.double(0))$c}
#cppfoo=function(a,b){.C('cppfoo',as.double(a),as.double(b),c=as.double(0))$c}

#### Wrap Doptimize function

#' Wrapper function for OApackage Doptimize function.
#'
#' This function generates
#' a single 2-level design with the specified number of runs and factors. The
#' design is optimized for the specified criterium. 
#'
#' The criterium that is optimized is:
#'	alpha1*D + alpha2*Ds + alpha3 * D1
#'
#' Here D is the D-efficiency of the design and Ds and D1 other efficiency measures.
#' For more details see http://pietereendebak.nl/oapackage/Doptim.html
#' 
#' @param N Number of runs
#' @param k Number of factors
#' @param nrestarts Number of restarts to generated an optimal design
#' @param alpha1 Parameter
#' @param alpha2 Parameter
#' @param alpha3 Parameter
#' @param verbose Integer that determines the amount of debug output
#' @param method Integer, default: 0
#' @return A matrixs containing the generated design
Doptimize=function(N, k, nrestarts, alpha1=1, alpha2=0, alpha3=0, verbose=1, niter=20000,  method=0, maxtime=500, nabort=6000) {

nn <- N*k
#print('call')
tmp <- .C('DoptimizeR', as.integer(N), as.integer(k), as.integer(nrestarts), as.integer(niter), as.double(alpha1), as.double(alpha2), as.double(alpha3), as.integer(verbose), as.integer(method), as.double(maxtime), as.integer(nabort), result=double(nn) ) 
#print(tmp)
#print('done')
p = tmp[['result']]

A <- array(p, dim=c(N,k) )

 }


# To compile documentation from the source code:
#
# setwd(...)
# devtools::document() 
