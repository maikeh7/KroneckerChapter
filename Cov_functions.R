covpow <- function(locs,pow=2,scale=5){
  d1 <- dist(locs)
  n <- dim(locs)[1]
  mat <- matrix(0,n,n)
  mat[lower.tri(mat)] <- d1
  mat <- mat+t(mat)
  cc <- exp(-(mat/scale)^pow)
  return(cc)
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


make_distmat = function(locs){
  d1 <- dist(locs)
  n <- dim(locs)[1]
  mat <- matrix(0,n,n)
  mat[lower.tri(mat)] <- d1
  mat <- mat+t(mat)
  return(mat)
}