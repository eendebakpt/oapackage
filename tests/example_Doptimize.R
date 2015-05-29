## Example code for calling Doptimize function from R

library(oapackage)

# set variables

N = 32   # number of runs
k=7      # number of factors
nrestarts=16  # number of restarts 
alpha1=1      # parameters of the optimization function
alpha2=1
alpha3=0
verbose=1    # specified the output level of the function

# run the function
p = Doptimize(N, k, nrestarts, alpha1, alpha2, alpha3, verbose ) 

print('resulting design:')
print(p)


