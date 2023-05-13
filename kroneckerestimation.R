Fast_kronecker_quadratic=function(C, R, Uc, Sc, Uc_t, Ur, Sr, Ur_t, z, nugget, marginal_variance){
  p=nrow(C)
  k=nrow(R)
  Zmat = matrix(z, nrow = k, ncol = p)
  Myvec = as.vector(Ur_t %*% Zmat %*% Uc) 
  
  # Scmin = min(Sc)
  # 
  # Srmin = min(Sr)
  # print(c(Scmin, Srmin))
  
  eigenvals = marginal_variance * as.vector((outer(Sr, Sc))) + nugget #ok this works for sure
  
  
  InvEigenMat = diag((1 / eigenvals), nrow = p*k, ncol = p*k)
  
  ssqKronecker = (t(Myvec) %*% InvEigenMat %*% Myvec) 
  LogDet = sum(log(eigenvals)) 
  
  return(list(quadratic = ssqKronecker, LogDet = LogDet))
}


LogLik_MVN <- function(C, R, Uc, Sc, Uc_t, Ur, Sr, Ur_t, z, nugget, marginal_variance, Temporal_nugget) {
  p=nrow(C) #temporal
  k=nrow(R) #spatial
  n = p*k #total length of z vec
  
  #calculates the quadratic term
  kronResult = Fast_kronecker_quadratic(C, R, Uc, Sc, Uc_t, Ur, Sr, Ur_t, z, nugget, marginal_variance) 
  ssqKronecker = kronResult$quadratic 
  
  LogDet = kronResult$LogDet
  #return log likelihood
  loglike = - .5*LogDet -n/2*log(2 *pi) - .5*ssqKronecker
  
  #SigSquared = 1000000
  # loglike = (-n/2) *log(2 *pi) - .5*LogDet - .5*ssqKronecker - 
  #    .5*SigSquared*(.285^2 - marginal_variance - nugget - Temporal_nugget)^2 -
  #    .5*SigSquared*(0.0807932 - nugget - Temporal_nugget)^2
  
  
  return(loglike)
}

Find_scale_params = function(parameters, distC, distR, z){
  P = 365
  #P= 55
  ell = parameters[1]
  #ell = 1.043448
  #scaleR = parameters[2]
  scaleR = 120
  nugget = parameters[2]
  marginal_variance = parameters[3]
  Temporal_nugget = parameters[4]
  
  
  #calculate eigendecomp of C (temporal ) and R (spatial) correlation matrices
  C = exp(-2*ell^2 * (sin(pi * distC/P))^2 ) + diag(Temporal_nugget, nrow = nrow(distC), ncol = ncol(distC))
  #Cmin = matrix(min(C), nrow = nrow(C), ncol = ncol(C))
  #C2 = C - Cmin
  #C3 = C2 / max(C2)
  #C=C3 + diag(Temporal_nugget, nrow = nrow(distC), ncol = ncol(distC))
  
  
  svdC = svd(C)
  Uc = svdC$u
  Sc = svdC$d
  Uc_t = t(svdC$v)
  
  R = exp(-(distR / scaleR)^2) 
  svdR = svd(R)
  Ur = svdR$u
  Sr = svdR$d
  Ur_t = t(svdR$v)
  
  MylogLik = LogLik_MVN(C, R, Uc, Sc, Uc_t, Ur, Sr, Ur_t, z, nugget, marginal_variance, Temporal_nugget)
  # print(MylogLik)
  return(-MylogLik) #return negative log likelihood
}

.26^2
nugget=.1
marginal_variance=.01
dmvnorm(x=z, mean = rep(0, nrow(SigGGT)), sigma = SigGGT, log = T)
#######################################################################################################
dim(distC)
length(z)
range(distC)
range(distR)
distC = sigmaT
distR = SigmaNN
z = monthAves$MeanResid
z = monthAves$MeanTMAXwhite

#try dfoptim --nelder mead method (this works better than optim)
#this converges and doesn't give an error message
library(dfoptim)
# ell = parameters[1]
# scaleR = parameters[2]
# nugget = parameters[3]
# marginal_variance = parameters[4]
# Temporal_nugget = parameters[5]

#0.672604743 47.258957106  0.007379727  0.999999993  0.013868505
testmkb = nmkb(c(900, .0001, .0001, .0001), Find_scale_params,distC = distC, distR=distR, z=z,
               lower = c(10, 0.00001, 0.00001, .00001), upper = c(1000, .1, .1, .1))