
#' eudis: calculate euclidean distance
#'
#' @param x x coordinate
#' @param y y coordinate
#'
#' @return distance between x and y coordinates
#' @export
#'
#' @examples
#' eudis(x = 1, y = 2)
eudis = function(x, y) { sqrt( sum( (x-y)^2 ) ) } 


#' get_distmat: construct a distance matrix (distance is euclidean distance)
#'
#' @param l1 vector of coordinates
#' @param l2 vector of coordinates
#'
#' @return a matrix of dimension length(l1) x length(l2)
#' @export
#'
#' @examples
#' get_distmat(l1 = 1:5, l2 = 6:15)
get_distmat = function(l1, l2){
  seq1 = 1:length(l1)
  seq2 = 1:length(l2)
  outer(seq1, seq2, FUN = Vectorize(function(x, y) eudis(l1[x], l2[y])))
}


#' rmultnorm generates realizations from a multivariate normal distribution
#'
#' @param n number of realizations to generate
#' @param mu mean vector
#' @param sigma a covariance matrix of dimension (length(mu) x length(mu))
#'
#' @return an n x length(mu) matrix of realizations
#' @export
#'
#' @examples
#' n = 1
#' mu = rep(0, 10)
#' mydist = get_distmat(l1 = 1:10, l2 = 1:10)
#' sigma = exp( -(mydist / 2)^2 )
#' mysim = rmultnorm(n = 1, mu = mu, sigma = sigma)
rmultnorm <- function(n,mu,sigma){
  p <- length(mu)
  z <- matrix(rnorm(n * p),nrow=n)
  svdsig <- svd(sigma)
  ch <- sqrt(diag(svdsig$d)) %*% t(svdsig$u)  
  zz <- z %*% ch
  zz + matrix(mu,nrow=n,ncol=p,byrow=T)
}


#' rmultnorm_SVD efficiently generates realizations from a MVN normal distribution where the covariance 
#' matrix is represented as kronecker(cov1, cov2)
#'
#' @param cov1 an nxn covariance matrix
#' @param cov2 a pxp covariance matrix
#'
#' @return an nxp vetor of realizations
#' @export
#'
#' @examples
#' t1 = seq(0,1,length = 11)
#' x1 = seq(0,1,length = 10)
#' 
#' Tdist = get_distmat(t1, t1)
#' Xdist = get_distmat(x1, x1)
#' 
#' Rt = exp( -(Tdist / .5)^2 )
#' Rx = exp( -(Xdist / 1.2)^2 )
#' sim = rmultnorm_SVD(cov1 = Rx, cov2 = Rt)
rmultnorm_SVD = function(cov1, cov2){
  R = cov1
  C = cov2
  p = dim(R)[1]
  k = dim(C)[1]
  svdR = svd(R)
  svdC = svd(C)
  R_U = svdR$u
  C_U = svdC$u
  
  z = rnorm(p*k)
  Rd = sqrt(svdR$d)
  Cd = sqrt(svdC$d)
  eigs = kronecker(Rd, Cd)
  ztilde = eigs*z
  mysim = as.vector(C_U %*% matrix(ztilde, nrow= k, ncol = p) %*% t(R_U))
  # slower alternative:
   # p = nrow(R) 
   # k = nrow(C) 
   # R_d = diag(sqrt(svdR$d))
   # C_d = diag(sqrt(svdC$d))
   # z = matrix(rnorm(p*k), nrow= k, ncol = p)
   # mu = matrix(rep(0, p*k), nrow = k, ncol = p)
   # h3 = (C_U %*% C_d) %*% z %*% t(R_U %*% R_d) + mu
   # mysim = as.vector(h3)
  
  return(mysim)
}
#cov1 = R_xstar_xstar; cov2 = Rt

