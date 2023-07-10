###########################################################################################
# function that carries out optimization
# f -- computer model output (vector)
# x1 -- space coordinates of f x1
# q --  number of bases (>0)
# This function returns a list where each element of the list are the phi_j's and sig2_j's
# corresponding to q_j
# nu is set to 1e-6
###########################################################################################
ML_pcEMU = function(fn, x, q=2){
  stuff = pre_process(fn, x, q)
  K = stuff$K
  paramList = list(NULL)
  for (i in 1:q){
    myParams = nmkb(par = c(1, 1), fn = Estimate_params_pcEm, lower=.001, upper=100,
                    myObject = stuff, iterNum = i)
    pvec = myParams$par; names(pvec)=c("phi","sigma^2")
    paramList[[i]] = pvec
  }
  return(paramList)
}



#####################################################################################
# main function for optimization
# myObject comes from pre_process(), 
# iterNum is current basis element (j = 1..q), 
#####################################################################################
Estimate_params_pcEm = function(parameters, myObject, iterNum){
  
  phi = parameters[1]
  sig2 = parameters[2]
  
  nu = myObject$nu 
  K = myObject$K
  fmat0 = myObject$fmat0
  fsvd = myObject$fsvd
  s = myObject$s
  nc = myObject$nc
  x1= myObject$x1
  Xdist = myObject$xdist
  L2 = myObject$L2
  
  kj = K[, iterNum]
  
  nug = nu / sum(kj^2)
  
  covmat_j  = sig2*exp(-phi*Xdist^2) + diag(rep(nug, nc))
  
  chCov_j = chol(covmat_j)
  
  LogDet = 2*sum(log(diag(chCov_j))) 
  
  what_j = fsvd$v[ ,iterNum] * sqrt(nc)
  
  vec = backsolve(chCov_j, what_j)
  
  ssq  = sum(vec^2)
  
  loglike = - .5*LogDet - .5*ssq - 0.5 * L2
  
  return(-loglike) #return negative log likelihood
}

######################################################
# pre process f (computer output), 
# make K and compute second part of likelihood
######################################################
pre_process = function(fn, x, q, nu=1e-6){
  # make f into matrix of appropriate dims
  nc = length(x)
  fmat = matrix(fn, ncol=nc)
  xdist = as.matrix(dist(x))
  
  # subtract means of rows
  meanf = apply(fmat, 1, mean)
  fmat0 = fmat - meanf
  
  # do svd
  fsvd = svd(fmat0)
  
  # make K
  K = fsvd$u[,1:q] %*% diag(fsvd$d[1:q]) / sqrt(nc)
  s = nrow(K)
  
  # necessary to calculate second chunk of likelihood
  Inc = diag(rep(1, nc))
  Kbig = cbind(kronecker(Inc, K[, 1]), kronecker(Inc, K[, 2]))
  
  fvec = as.vector(fmat0)
  
  # this is stupid but wtf can't quickly see another way to get this 
  # fcalc2 = fvec %*% (diag(1, nrow = length(f), ncol = length(f)) - 
  #                  Kbig %*% solve(t(Kbig) %*% Kbig,t(Kbig)))  %*% fvec
  # v1 = Kbig %*% solve(t(Kbig) %*% Kbig,t(Kbig)%*%fvec) 
  # v2 = fvec - as.vector(v1)
  # fcalc = sum(fvec*v2)
  
  # here's a simpler (coding-wise) calculation
  aa = lm(fvec ~ Kbig - 1)
  fcalc = sum(aa$residuals^2)
  
  L2 = ((s - q)*nc / 2) * log(nu) - 1/2 * (1 / nu) * fcalc
  
  return(list(L2= L2, K = K, fmat0=fmat0, fsvd=fsvd, s = s, nu = nu,
              nc = nc, x=x, xdist=xdist, q=q))
}
