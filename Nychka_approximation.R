source("Cov_functions.R")
source("kroneckerestimation.R")
nt = 11; nx = 10       # good for plots
#n1 = 11*4; n2 = 10*4    # good for timing
t1 = seq(0,1,length=nt)
x1 = seq(0,1,length=nx)
x = expand.grid(t1,x1)
Covt = covpow(cbind(t1,0),scale=1)
Covx = covpow(cbind(0,x1),scale=1)
Cov_x_t = kronecker(Covx, Covt)

pow = 2
scale = 1
nstar = 25
xstar = seq(0,1, length = nstar)

# grid for plotting 
plot_grid = expand.grid(t1, xstar)

xstar_grid = c(x1, xstar)
Cov_xstar_xstar = covpow(cbind(0,xstar_grid), scale=1)

# this is Nychka's sigma
Sigma22 = kronecker(Cov_xstar_xstar, Covt)

train_idx = 1:length(x1)
train_idx_kron = 1:(length(train_idx) * length(t1))
test_idx_kron = (train_idx_kron[length(train_idx_kron)] + 1): nrow(Sigma22)

kmat = matrix(0, nrow = length(train_idx_kron), ncol = nrow(Sigma22))
dim(kmat)
# b/c the training samples are stacked first, kmat will be identiy matrix for
# the first several row/cols and then just zeros
kmat[train_idx_kron, train_idx_kron] = diag(1, nrow= length(train_idx_kron))

# nugget-sigma^2 in Nychka's paper
sigmasq = 0.001

#checked--this is kronecker(covx, covt)
Sigma11 = kmat %*% Sigma22 %*% t(kmat) + diag(sigmasq,
                    nrow = nrow(Covx)*nrow(Covt), ncol = nrow(Covx)*nrow(Covt))
dim(Sigma11)
Sigma12 = Sigma22 %*% t(kmat)
dim(Sigma12)
Sigma21 = kmat %*% Sigma22
dim(Sigma21)

# make a 'real' y
real_y = rmultnorm(1, mu=rep(0, nrow(Sigma11)), Sigma11)

# Nychka's zhat --  conditional mean
zhat = Sigma12 %*% solve(Sigma11) %*% as.vector(real_y)

# make the actual conditional covariance matrix
True_cond_cov =  Sigma22 - Sigma12 %*% solve(Sigma11) %*% Sigma21

# True realization from the conditional cov matrix
true_realization = as.vector(rmultnorm(1, zhat, True_cond_cov))
persp(t1, xstar, matrix(true_realization[test_idx_kron], nrow = nt),theta = 130-90, phi = 10,
      xlab='t',ylab='x',zlab='f',zlim=c(-5,5)) -> res
points(trans3d(plot_grid[,1], plot_grid[,2], true_realization[test_idx_kron], pmat = res),
       col = 'black', pch = 16,cex=.7)

# for comparison to Nychka
Lots_of_realizations = rmultnorm(200, zhat, True_cond_cov)
dim(Lots_of_realizations)
est_mean_true =apply(Lots_of_realizations, 2, mean)
est_sds_true = apply(Lots_of_realizations, 2, sd)

#############################
# vectorized but still slow #
#############################
epsilon = t(rmultnorm(1000, rep(0, nrow(Sigma11)), diag(sigmasq, nrow = nrow(Sigma11))))
u = rmultnorm(200, mu = rep(0, nrow(Sigma22)), Sigma22)
y_star = kmat %*% t(u) + epsilon
ustar = t(u) - Sigma12 %*% solve(Sigma11) %*% y_star 
conditional_samps = as.vector(zhat) + ustar
test = t(apply(ustar, 2, function(x) x + zhat))
est_mean = apply(test, 2, mean)
est_sds = apply(test, 2, sd)

# set up
R = Cov_xstar_xstar
C = Covt
p = dim(R)[1]
k = dim(C)[1]
svdR = svd(R)
svdC = svd(C)
R_U = svdR$u
C_U = svdC$u
R_d = diag(sqrt(svdR$d))
C_d = diag(sqrt(svdC$d))

pre_process_svd = function(z, p, k, C_U, C_d, R_U, R_d){
  zmat = matrix(z, nrow = k, ncol= p)
  res = as.vector((C_U %*% C_d) %*% zmat %*% t(R_U %*% R_d) )
}
allZs = matrix(rnorm(p*k*200), nrow=p*k)

umat = apply(allZs,2, function(x) testfun(z=x, p=p, k=k, C_U = C_U, C_d = C_d, R_U = R_U, R_d = R_d))



persp(t1, xstar, matrix(test[test_idx_kron], nrow = nt),theta = 130-90, phi = 10,
      xlab='t',ylab='x',zlab='f',zlim=c(-2,2)) -> res
points(trans3d(plot_grid[,1], plot_grid[,2], test[test_idx_kron], pmat = res),
       col = 'black', pch = 16,cex=.7)

persp(t1, xstar, matrix(comp[test_idx_kron], nrow = nt),theta = 130-90, phi = 10,
      xlab='t',ylab='x',zlab='f',zlim=c(-2,2)) -> res
points(trans3d(plot_grid[,1], plot_grid[,2], comp[test_idx_kron], pmat = res),
       col = 'black', pch = 16,cex=.7)


rootT = svd_T$u %*% sqrt(diag(svd_T$d))
rootXX = svd_xx$u %*% sqrt(diag(svd_xx$d))
dim(rootXX)

test= rootT %*% Zmat %*% rootXX
range(test)
dim(Zmat)
hist(test[test_idx_kron])

svd_T = svd(Covt)
svd_xx = svd(Cov_xstar_xstar)

test2=svd_T$u %*% diag(svd_T$d) %*% t(svd_T$v)
image(test2)
image(Chol_t)
# try this tomorrow
test = rmultnorm_Cholesky(R= Cov_xstar_xstar, C = Covt)
test = rmultnorm_SVD(Cov_xstar_xstar, Covt)
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
  #spatial = rows, temporal = columns
  z = matrix(rnorm(p*k), nrow=k, ncol = p)
  
  #mean matrix has same dimensions as z
  mu = matrix(rep(0, p*k), nrow = k, ncol = p)
  mysim = mu + t(LC) %*% z %*% LR
  
  #h2 = mu + %*% z %*% t(LsigT) 
  mysim = as.vector(mysim)
  return(mysim)
}

start = Sys.time()
for (i in 1:200){
  sim = rmultnorm_SVD(R, C)
}
stop = Sys.time()
stop - start
start = Sys.time()
sim = rmultnorm(200, mu = rep(0, nrow(Sigma22)), Sigma22)
stop = Sys.time()
stop - start

dim(zmat)
pre_process_svd = function(z, p, k, C_U, C_d, R_U, R_d){
  zmat = matrix(z, nrow = k, ncol= p)
  res = as.vector((C_U %*% C_d) %*% zmat %*% t(R_U %*% R_d) )
}
testfun(z=z, p=p, k=k, C_U = C_U, C_d = C_d, R_U = R_U, R_d = R_d)
z = rnorm(p*k)
R = Cov_xstar_xstar
C = Covt
allZs = matrix(rnorm(p*k*200), nrow=p*k)
stop = Sys.time()
testing = apply(allZs,2, function(x) testfun(z=x, p=p, k=k, C_U = C_U, C_d = C_d, R_U = R_U, R_d = R_d))
dim(testing)
stop = Sys.time()
stop - start


Myapplyfunction <- function(){
    test=apply(allZs,2, function(x) testfun(z=x, p=p, k=k, C_U = C_U, C_d = C_d, R_U = R_U, R_d = R_d))
}

forloopfunction = function(){
  umat = matrix(NA, ncol = p*k, nrow = 200)
  for (i in 1:200){
    umat[i, ] = rmultnorm_SVD(R, C)
  }
}

regularFunc = function(){
  test=t(rmultnorm(1000, rep(0, nrow(Sigma11)), diag(sigmasq, nrow = nrow(Sigma11))))
}

res <- microbenchmark(NULL, f(), times=100L)

## Print results:
print(res)
library(microbenchmark)
results = microbenchmark(forloopfunction(), Myapplyfunction(), regularFunc()) 
results






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

rmultnorm <- function(n, mu, sigma){
  # returns n rows of multinormal mu,sigma random vectors
  # ie: returns a n x length(mu) matrix
  p <- length(mu)
  z <- matrix(rnorm(n * p), nrow=n)
  ch <- chol(sigma, pivot=T)
  piv <- attr(ch, "pivot")
  zz <- (z%*%ch)
  zzz <- 0*zz
  zzz[, piv] <- zz
  zzz + matrix(mu, nrow=n, ncol=p, byrow=T)
}


dvals = sqrt(outer(svd_xx$d, svd_T$d))
dim(svd_xx$u)

