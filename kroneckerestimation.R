#####################################################################
### estimating parameters!
#####################################################################
# see Find_scale_params() below
# I found that functions from dfoptim worked better than from optim
# this is for kronecker(C, R), not the other way around!
# starting parameters goes first, then name of function to optimize
# distC = distance matrix for C
# distR = distnace matrix for R
# z is response vector
# lower =  lower bounds
# upper = upper bounds
# nelder mead method

library(dfoptim)
Fast_kronecker_quadratic=function(C, R, Uc, Sc, Uc_t, Ur, Sr, Ur_t, z, nugget, marginal_variance){
  p = nrow(C)
  k = nrow(R)
  Zmat = matrix(z, nrow = k, ncol = p)
  # could be faster using Dave's solve()...
  Myvec = as.vector(Ur_t %*% Zmat %*% Uc) 
  
  eigenvals = marginal_variance * as.vector((outer(Sr, Sc))) + nugget #ok this works for sure
  
  
  InvEigenMat = diag((1 / eigenvals), nrow = p*k, ncol = p*k)
  
  ssqKronecker = (t(Myvec) %*% InvEigenMat %*% Myvec) 
  LogDet = sum(log(eigenvals)) 
  
  return(list(quadratic = ssqKronecker, LogDet = LogDet))
}


LogLik_MVN <- function(C, R, Uc, Sc, Uc_t, Ur, Sr, Ur_t, z, nugget, marginal_variance) {
  p = nrow(C) 
  k = nrow(R) 
  n = p*k #total length of z vec
  
  # calculates the quadratic term
  kronResult = Fast_kronecker_quadratic(C, R, Uc, Sc, Uc_t, Ur, Sr, Ur_t, z, nugget, marginal_variance) 
  ssqKronecker = kronResult$quadratic 
  
  LogDet = kronResult$LogDet
  
  # return log likelihood
  loglike = - .5*LogDet -n/2*log(2 *pi) - .5*ssqKronecker

  return(loglike)
}

# distC and distR are distance matrices
# z is response
# parameters are ... parameters!
# this function assumes you want kronecker(C, R), NOT kronecker(R, C)
Estimate_params = function(parameters, distC, distR, z){
  # parameters we want to estimate
  scaleC = parameters[1]
  scaleR = parameters[2]
  nugget = parameters[3]
  marginal_variance = parameters[4]
  
  #calculate eigendecomp of C and R 
  C = exp(-(distC / scaleC)^2) 
  R = exp(-(distR / scaleR)^2) 
  
  # do svd (eigen) decomp of both cov mats
  svdC = svd(C)
  Uc = svdC$u
  Sc = svdC$d
  Uc_t = t(svdC$v)
  
  R = exp(-(distR / scaleR)^2) 
  svdR = svd(R)
  Ur = svdR$u
  Sr = svdR$d
  Ur_t = t(svdR$v)
  
  MylogLik = LogLik_MVN(C, R, Uc, Sc, Uc_t, Ur, Sr, Ur_t, z, nugget, marginal_variance)
  # print(MylogLik)
  return(-MylogLik) #return negative log likelihood
}

# a function to make the Gaussian covariance matrix
covpow <- function(locs,pow=2,scale=5){
  # browser()
  d1 <- dist(locs)
  n <- dim(locs)[1]
  mat <- matrix(0,n,n)
  mat[lower.tri(mat)] <- d1
  mat <- mat+t(mat)
  cc <- exp(-(mat/scale)^pow)
  return(cc)
}

make_distmat = function(locs){
  d1 <- dist(locs)
  n <- dim(locs)[1]
  mat <- matrix(0,n,n)
  mat[lower.tri(mat)] <- d1
  mat <- mat+t(mat)
  return(mat)
}

############
# Example 1
############
# try estimating using our deterministic fucntion f
# this does not work so well...
n1 = 55; n2 = 50
t1=seq(0,1,length=n1)
x1=seq(0,1,length=n2)
x=expand.grid(t1,x1)

# a deterministic function to be emulated
# adding noise seems to help estimation
f = (x[,2]+1)*cos(pi*x[,1]) + rnorm(nrow(x), sd= 0.0001)

# C = x
# R = t
Tdist = make_distmat(cbind(t1,0))
Xdist = make_distmat(cbind(0, x1))
persp(t1,x1,matrix(f,nrow=n1),theta = 130-90, phi = 10,xlab='t',ylab='x',zlab='f',zlim=c(-3,3)) -> res
points(trans3d(x[,1], x[,2], f, pmat = res), col = 'black', pch = 16,cex=.7)

testmkb = nmkb(c(0.5, .5, .001, 1), Estimate_params, distC = Xdist, distR = Tdist, z=f,
               lower = c(0.1, 0.1, 0.00000001, .5), upper = c(50, 10, 1, 5))
# HMMMMMM.....not sure these make much sense....
testmkb

#################
# Example 2
#################
# Make up some data that comes from the true generating process 
n1 = 44; n2 = 40
t1=seq(0,1,length=n1)
x1=seq(0,1,length=n2)
x=expand.grid(t1,x1)
Tdist = make_distmat(cbind(t1,0))
Xdist = make_distmat(cbind(0, x1))

# NOTE
# C = Covx
# R = Covt
Covt = covpow(cbind(t1,0),scale=.5)
k=nrow(Covt)
Covx = covpow(cbind(0,x1),scale=.8)
p=nrow(Covx)
Cnug = .0003
margvar = 1.5
KronCovmat = margvar * kronecker(Covx, Covt) +  diag(Cnug, nrow = p*k, ncol = p*k) 
mySim = rmultnorm(1, rep(0, p*k), KronCovmat)
mySim = as.vector(mySim)
#persp(t1,x1,matrix(mySim,nrow=n1),theta = 130-90, phi = 10,xlab='t',ylab='x',zlab='f',zlim=c(-3,3)) -> res
#points(trans3d(x[,1], x[,2], mySim, pmat = res), col = 'black', pch = 16,cex=.7)
#range(mySim)


# ok, this seems to work pretty well. Better to use more data!
testmkb = nmkb(c(0.5, .5, .001, 1), Estimate_params, distC = Xdist, distR = Tdist, z=mySim,
               lower = c(0.1, 0.1, 0.000001, .5), upper = c(10, 10, 1, 5))
testmkb
#dmvnorm(x=z, mean = rep(0, nrow(SigGGT)), sigma = SigGGT, log = T)
#######################################################################################################
