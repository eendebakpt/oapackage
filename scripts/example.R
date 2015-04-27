## Example code for calling Doptimize function from R
library(devtools)
load_all('oapackage')

# set variables

N = 40
k=7
N = 8
k=3
nrestarts=16
niter=10000
alpha1=1
alpha2=1
alpha3=0
verbose=1
method=0
maxtime=100 # max running time in seconds

p = Doptimize(N, k, nrestarts, alpha1, alpha2, alpha3, verbose, niter, method, maxtime ) 

print('resulting design:')
print(p)


