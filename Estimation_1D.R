#####################################################################
### Parameter estimation
#####################################################################
# I found that functions from dfoptim worked better than from optim
# this is for kronecker(C, R), not the other way around!
# starting parameters goes first, then name of function to optimize
# distC = distance matrix for C
# distR = distnace matrix for R
# z is response vector
# lower =  lower bounds
# upper = upper bounds
# nelder mead method
######################################################################
library(dfoptim)


#' Fast_kronecker_quadratic: efficiently carries out the quadratic equation for MVN likelihood as well as the log determinant
#'
#' @param C covariance matrix C
#' @param R covariance matrix R
#' @param Uc svd(C)$u
#' @param Sc svd(C)$d
#' @param Uc_t t(svd(C)$u)
#' @param Ur svd(R)$u
#' @param Sr svd(R)$d
#' @param Ur_t t(svd(R)$u)
#' @param z vector of data
#' @param nugget nugget value
#' @param marginal_variance marginal variance value
#'
#' @return a list containing the quadratic calculation and log determinant
#' @export 
#'
#' @examples
Fast_kronecker_quadratic=function(C, R, Uc, Sc, Uc_t, Ur, Sr, Ur_t, z, nugget, marginal_variance){
  p = nrow(C)
  k = nrow(R)
  Zmat = matrix(z, nrow = k, ncol = p)
  
  Myvec = as.vector(Ur_t %*% Zmat %*% Uc) 
  
  eigenvals = marginal_variance * as.vector((outer(Sr, Sc))) + nugget #ok this works for sure
  inv_eigs = 1/eigenvals
  
  ssqKronecker = (Myvec)^2 %*% inv_eigs
  
  # slower:
  # InvEigenMat = diag((1 / eigenvals), nrow = p*k, ncol = p*k)
  # ssqKronecker = (t(Myvec) %*% InvEigenMat %*% Myvec) 
  
  LogDet = sum(log(eigenvals)) 
  
  return(list(quadratic = ssqKronecker, LogDet = LogDet))
}


#' LogLik_MVN: this function calls Fast_kronecker_quadratic() and calculates the log likelihood
#'
#' @param C covariance matrix C
#' @param R covariance matrix R
#' @param Uc svd(C)$u
#' @param Sc svd(C)$d
#' @param Uc_t t(svd(C)$u)
#' @param Ur svd(R)$u
#' @param Sr svd(R)$d
#' @param Ur_t t(svd(R)$u)
#' @param z vector of data
#' @param nugget nugget value
#' @param marginal_variance marginal variance value
#'
#' @return log likelihood
#' @export
#'
#' @examples
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

#' Estimate_params: estimate parameters for 1D estimation of marginal variance, scales, and nugget given data z. This
#' function assumes the udnerlying distribution is MVN(0, kronecker(C, R))
#'
#' @param parameters vector of parameters to be estimated
#' @param distC distance matrix for covariance matrix C
#' @param distR distance matrix for covariance matrix R
#' @param z vector of data
#'
#' @return negative log likelihood value
#' @export
#'
#' @examples
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
