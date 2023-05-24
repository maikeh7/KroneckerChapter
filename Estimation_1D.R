#####################################################################
### Parameter estimation
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