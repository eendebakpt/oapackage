# Testing functions
#foo <- function(x){x*2}

#cfoo=function(a,b){.C('cfoo',as.double(a),as.double(b),c=as.double(0))$c}
#cppfoo=function(a,b){.C('cppfoo',as.double(a),as.double(b),c=as.double(0))$c}

#### Wrap Doptimize function

#' Wrapper function for OApackage Doptimize function
Doptimize=function(N, k, nrestarts, alpha1, alpha2, alpha3, verbose=1, niter=20000,  method=0, maxtime=500, nabort=6000) {

nn <- N*k
#print('call')
tmp <- .C('DoptimizeR', as.integer(N), as.integer(k), as.integer(nrestarts), as.integer(niter), as.double(alpha1), as.double(alpha2), as.double(alpha3), as.integer(verbose), as.integer(method), as.double(maxtime), as.integer(nabort), result=double(nn) ) 
#print(tmp)
#print('done')
p = tmp[['result']]

A <- array(p, dim=c(N,k) )

 }


