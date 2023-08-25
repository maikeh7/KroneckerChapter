#' ML_pcEMU: function that estimates phi_j's and sig2_j's using maximum likelihood
#'
#' @param fn computer model output
#' @param x inputs to computer model that are to be varied
#' @param q number of bases desired
#'
#' @return a list where each element of the list are the phi_j's and sig2_j's corresponding to q_j (note nu is set to 1e-6)
#' @export
#'
#' @examples
#' s = 11         # size of output space
#' nc = 10        # number of computer model runs

#' t1 = seq(0,1,length = s)
#' x1 = seq(0,1,length = nc)
#' x = expand.grid(t1, x1)
#'
#' # function to emulate
#' f = (x[, 2]+1)*cos(pi*x[, 1]) + .03*(exp(x[, 2]))
#' 
#' params = ML_pcEMU(fn = f, x = x1, q = 2)
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
# main function for optimization -- called by ML_pcEMU
# parameters: parameters to be estimated
# myObject: object returned from pre_process(), 
# iterNum: current basis element (1..q) 
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
  
  covmat_j  = sig2*exp(-(phi*Xdist)^2) + diag(rep(nug, nc))
  
  chCov_j = chol(covmat_j)
  
  LogDet = 2*sum(log(diag(chCov_j))) 
  
  what_j = fsvd$v[ ,iterNum] * sqrt(nc)
  
  vec = backsolve(chCov_j, what_j)
  
  ssq  = sum(vec^2)
  
  loglike = - .5*LogDet - .5*ssq - 0.5 * L2
  
  return(-loglike) #return negative log likelihood
}

#' pre_process: pre process computer model output to facilitate parameter estimation. This function is called by ML_pcEMU()
#'
#' @param fn computer model output
#' @param x inputs to computer model that are to be varied
#' @param q number of bases desired
#' @param nu nugget, default is 1e-6
#'
#' @return a list containing components necessary for estimating parameters via Estimate_params_pcEm()
#' @export
#'
#' @examples
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
  
  klist = list()
  for (i in 1:q){
    klist[[i]] = kronecker(Inc, K[, i])
  }
  Kbigtest = do.call("cbind", klist)
  
  fvec = as.vector(fmat0)
  
  # SLOW VERSION
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


#' get_wstar_distr_preds: estimate phi and sig2 for w*_j
#'
#' @param f f =  function output at train input locations
#' @param xnew vector of new input locations
#' @param nc number of computer runs
#' @param q number of bases
#' @param x1 original x train input locations
#' @param t1 original t train input locations
#' @param params object returned from you get from ML_pcEMU
#'
#' @return a list containing means and covariance matrices for w*_j's and predictions at xnew and 
#' the f matrix as well as the mean of f
#' @export
#'
#' @examples
get_wstar_distr_preds = function(f, xnew, nc, q, x1,t1, params){
  
  nc_new = length(xnew)
  
  fmat = matrix(f, ncol=nc)
  
  xdist_aug = as.matrix(dist( c(x1, xnew)) )
  
  # subtract means of rows
  meanf = apply(fmat, 1, mean)
  fmat0 = fmat - meanf
  
  # do svd
  fsvd = svd(fmat0)
  
  # make K
  K = fsvd$u[,1:q] %*% diag(fsvd$d[1:q]) / sqrt(nc)
  
  
  s = nrow(K)
  wlist = list()
  
  preds = rep(0, length(t1)*nc_new)
  for (i in 1:q){
    phi = params[[i]][1]
    sig2 = params[[i]][2]
    
    nu = 1e-6
    
    kj = K[, i]
    
    Rcov = sig2*exp(-(phi*xdist_aug)^2)
    nug = nu / sum(kj^2)
    
    # get the parts of V
    sigma11 = nug * diag(1, nrow = nc, ncol = nc) +  Rcov[1:nc, 1:nc] 
    sigma12 = Rcov[1:nc,(nc+1):(nc+nc_new) ]
    sigma22 = Rcov[(nc+1):(nc+nc_new), (nc+1):(nc+nc_new) ]
    sigma21 = t(sigma12)
    
    sig11_inv = solve(sigma11)
    
    what_j = fsvd$v[ ,i] * sqrt(nc) 
    
    wmean = sigma21 %*% sig11_inv %*% what_j
    
    wcov = sigma22 - sigma21 %*%  sig11_inv %*% sigma12
    sublist = list(wmean = wmean, wcov = wcov)
    wlist[[i]] = sublist
    
  }
  
  return(list(wlist = wlist, K = K, meanf = meanf, fmat0 = fmat0))
  
}
