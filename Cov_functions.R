covpow <- function(locs,pow=2,scale=5){
  d1 <- dist(locs)
  n <- dim(locs)[1]
  mat <- matrix(0,n,n)
  mat[lower.tri(mat)] <- d1
  mat <- mat+t(mat)
  cc <- exp(-(mat/scale)^pow)
  return(cc)
}

eudis = function(x, y) { sqrt( sum( (x-y)^2 ) ) } 
calc_general_distmat = function(locs1, locs2){
  seq1 = 1:length(locs1)
  seq2 = 1:length(locs2)
  outer(seq1, seq2, FUN=Vectorize(function(x, y) eudis(locs1[x], locs2[y])))
}

# function to generate MVN realizations (using svd decomposition, not cholesky!)
rmultnorm <- function(n,mu,sigma){
  # returns n rows of multinormal mu, sigma random vectors
  # ie: returns a n x length(mu) matrix
  p <- length(mu)
  z <- matrix(rnorm(n * p),nrow=n)
  # ch <- chol(sigma,pivot=T)
  svdsig <- svd(sigma)
  ch <- sqrt(diag(svdsig$d)) %*% t(svdsig$u)  
  #piv <- attr(ch,"pivot")
  #zz <- (z%*%ch)
  zz <- z %*% ch
  #zzz <- 0*zz
  #zzz[,piv] <- zz
  zz + matrix(mu,nrow=n,ncol=p,byrow=T)
}

# generate efficient realizations from 0 mean GP
# w/ covariance matrix kronecker(R, C)
# may not work if matrix is ill conditioned (use rmultnorm_svd())
rmultnorm_Cholesky = function(R, C){
  LR = chol(R, pivot=T)
  LC = chol(C, pivot =T)
  pivR = attr(LR, "pivot")
  pivC = attr(LC, "pivot")
  ordR = order(pivR)
  ordC = order(pivC)
  LR = LR[,ordR]
  LC = LC[,ordC]
  p = dim(R)[1]
  k = dim(C)[1]
  
  z = matrix(rnorm(p*k), nrow=k, ncol = p)
  
  #mean matrix has same dimensions as z
  mu = matrix(rep(0, p*k), nrow = k, ncol = p)
  mysim = mu + t(LC) %*% z %*% LR
  
  #h2 = mu + %*% z %*% t(LsigT) 
  mysim = as.vector(mysim)
  return(mysim)
}

# generate efficient realizations from 0 mean GP
# w/ covariance matrix kronecker(R, C)
# using svd rather than cholesky
rmultnorm_SVD = function(R, C){
  p = dim(R)[1]
  k = dim(C)[1]
  svdR = svd(R)
  svdC = svd(C)
  R_U = svdR$u
  C_U = svdC$u
  R_d = diag(sqrt(svdR$d))
  C_d = diag(sqrt(svdC$d))
  z = matrix(rnorm(p*k), nrow= k, ncol = p)
  mu = matrix(rep(0, p*k), nrow = k, ncol = p)
  h3 = (C_U %*% C_d) %*% z %*% t(R_U %*% R_d) + mu
  mysim = as.vector(h3)
  return(mysim)
}

make_distmat = function(locs){
  d1 <- dist(locs)
  n <- dim(locs)[1]
  mat <- matrix(0,n,n)
  mat[lower.tri(mat)] <- d1
  mat <- mat+t(mat)
  return(mat)
}