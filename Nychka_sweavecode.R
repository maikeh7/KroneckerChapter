source("Cov_functions.R")
source("Estimation_1D.R")
source("Estimation_Basis.R")
# make the 2-d x support
s = 11; nc = 10     

t1 = seq(0,1,length = s)
x1 = seq(0,1,length = nc)
x = expand.grid(t1, x1)

# a deterministic function to be emulated
f = (x[,2]+1)*cos(pi*x[,1]) + .03*(exp(x[, 2]))

persp(t1, x1, 
      matrix(f, nrow = s),
      theta = 130-90, 
      phi = 10,
      xlab = 't', ylab ='x', zlab = 'f',
      zlim = c(-2.2, 2.4)) -> res

points(trans3d(x[, 1], x[, 2], 
               f,
               pmat = res),
       col = 'black',
       pch = 16,
       cex = 0.7)

Tdist = get_distmat(t1, t1)
Xdist = get_distmat(x1, x1)

# this function assumes the kronecker product: kronecker(C, R), so order matters
library(dfoptim)
params = nmkb(c(0.5, .5, .0001, .5), 
              Estimate_params, distC = Xdist, distR = Tdist, z=f,
              lower = c(0.1, 0.1, 1e-6, .1), upper = c(100, 100, 100, 100))

print(params$par[1])
print(params$par[2])
print(params$par[3]) # nu = nugget
print(params$par[4]) # sig2 = margvar

scaleX = params$par[1]
scalet = params$par[2]

Rt = exp( -(Tdist / scalet)^2 )
Rx = exp( -(Xdist / scaleX)^2 )
#Covt = make_covmat_Gauss(distmat = Tdist, scale = scalet, pow = 2) 
#Covx = make_covmat_Gauss(distmat = Xdist, scale = scaleX, pow = 2)

# prediction grid for x_star
#xnew = seq(0,1, length = nstar)

xnew = c(.13, .72, .97)
# in Nychka's setup, combine train and test points
xstar_grid = c(x1, xnew)

# This will be used for getting Nychka's 'Sigma'
# (Sigma = kronecker(Cov_(x*,x*), Cov_t))
dist_xstar_xstar = get_distmat(xstar_grid, xstar_grid)

#Cov_xstar_xstar = make_covmat_Gauss(distmat = dist_xstar_xstar, scale = scaleX, pow = 2)
R_xstar_xstar = exp( -(dist_xstar_xstar / scaleX)^2 )

# this is Nychka's Sigma
length(xstar_grid) * s

# build cross covariance matrix
dist_xstar_x = get_distmat(xstar_grid, x1)

#Cov_xstar_x = make_covmat_Gauss(distmat = distmat_xstar_x, scale = scaleX, pow = 2) 

R_xstar_x = exp( -(dist_xstar_x / scaleX)^2 )


svdx = svd(Rx)
svdt = svd(Rt)
svdG = svd(R_xstar_x)

Ux = svdx$u
Sx = svdx$d
Ux_t = t(svdx$v)

Ut = svdt$u
St = svdt$d
Ut_t = t(svdt$v)


#p = nc
#k = s
#p = nrow(C) 
#k = nrow(R)  
nu = params$par[3] #sigma_sq
sig2 = params$par[4] #margvar

# inverse of kronecker of singular values of C and R
# equivalent to kronecker(Sc, Sr)! 
lambda=diag(1 / (sig2*kronecker(Sx, St) + nu))

# reshape the response to be k * p
zreshape = matrix(f, nrow = s, ncol = nc)

res1 = lambda %*% as.vector(Ut_t %*% zreshape %*% Ux)

resmat = matrix(res1, nrow = s, ncol = nc)

cond_mean = sig2 * as.vector(Rt %*% Ut %*% resmat %*% Ux_t %*% t(R_xstar_x) )

# equivalent to cond_mean--Dave has this in the chapter
#test = margvar * as.vector(Ur %*% diag(Sr, nrow=nrow(Ur)) %*% resmat %*% Uc_t %*% t(G) )
#head(cond_mean)
#test- cond_mean

# grid for plotting 
# note we only plot test locations here
plot_grid = expand.grid(t1, xnew)
plot_grid_train = expand.grid(t1, x1)

# this is for differentiating train/test indices
train_idx = 1:length(x1)
train_idx_kron = 1:(length(train_idx) * length(t1))
# the end was previously nrow(Sigma22), but we don't actually ever want to make Sigma22
test_idx_kron = (train_idx_kron[length(train_idx_kron)] + 1): (length(xstar_grid) * s)


persp(t1, x1,
      matrix(cond_mean[train_idx_kron], nrow = s),
      theta = 130-90,
      phi = 10,
      xlab = 't', ylab = 'x', zlab = 'f', 
      zlim = c(-2.2, 2.4)) -> res

points(trans3d(plot_grid_train[,1],
               plot_grid_train[,2], 
               cond_mean[train_idx_kron],
               pmat = res),
       col = 'black',
       pch = 16,
       cex = 0.7)

points(trans3d(plot_grid[, 1],
               plot_grid[, 2], 
               cond_mean[test_idx_kron],
               pmat = res),
       col = 'red',
       pch = 16,
       cex = 0.7)
legend("topright", legend = c("train", "test"), col = c("black", "red"), pch=19)

# can be a little faster if we do this
# set up--this is based on rmultnorm_SVD in Cov_functions.R

svdR_xstar_xstar = svd(R_xstar_xstar)
Rstst_U = svdR_xstar_xstar$u
Rststd = svdR_xstar_xstar$d

svdt = svd(Rt)
Ut = svdt$u
Ut_t = t(svdt$v)
St = svdt$d

svdx = svd(Rx)
Sx = svdx$d
Ux = svdx$u
Ux_t = t(svdx$v)

#svdG = svd(R_xstar_x)

R_xstar_x # this is Sigma12
nc_prime = length(Rststd) # 13

nreals = 50
pre_process_svd = function(z, s, Ut, St, Rstst_U, Rststd){
 # z = rnorm(s*nc_prime)
  nc_prime = length(Rststd)
  # need to take square roots!!
  eigs = kronecker(sqrt(Rststd), sqrt(St))
  ztilde = eigs*z
  #zmat = matrix(z, nrow = s, ncol= nc)
  #res = as.vector((Rt_U %*% Rt_d) %*% zmat %*% t(Rstst_U %*% Rstst_d) )
  mysim = as.vector(Ut %*% matrix(ztilde, nrow = s, ncol = nc_prime) %*% t(Rstst_U))
}

# generate a matrix of N(0,1)'s
allZs = matrix(rnorm(s*nc_prime*nreals), nrow=s*nc_prime)

# here we draw from Sigma22, but note that we do not
# actually construct Sigma22!
kmat = matrix(0, nrow = length(train_idx_kron), ncol = s*nc_prime)
kmat[train_idx_kron, train_idx_kron] = diag(1, nrow= length(train_idx_kron))

epsilon = t(rmultnorm(nreals, rep(0, s*nc), diag(nu, nrow = s*nc)))

# use the pre_process_svd() function to get u and collect into a matrix
umat = apply(allZs,2, function(x) pre_process_svd(z=x, s=s, 
                                                  Ut = Ut, St = St, Rstst_U = Rstst_U, Rststd = Rststd))
# get y_tilde
y_tilde = kmat %*% umat + epsilon
# get w_star
w_star = get_cond_mean(svdt, svdx, R_xstar_x, y_tilde, nu, sig2)
ustar = umat - w_star
conditional_samps = as.vector(cond_mean) + ustar


persp(t1, x1,
      matrix(cond_mean[train_idx_kron], nrow = s),
      theta = 130-90,
      phi = 10,
      xlab = 't', ylab = 'x', zlab = 'f', 
      zlim = c(-2.2, 2.4)) -> res

for (i in 1:nreals){
  points(trans3d(plot_grid_train[,1],
                 plot_grid_train[,2], 
                 conditional_samps[train_idx_kron, i],
                 pmat = res),
         col = 'black',
         pch = 16,
         cex = 0.7)
  
  points(trans3d(plot_grid[, 1],
                 plot_grid[, 2], 
                 conditional_samps[test_idx_kron, i],
                 pmat = res),
         col = 'red',
         pch = 16,
         cex = 0.7)
}


legend("topright", legend = c("train", "test"), col = c("black", "red"), pch=19)



get_cond_mean = function(svdt, svdx, R_xstar_x, y_tilde, nu, sig2){

  Ut = svdt$u
  Ut_t = t(svdt$v)
  St = svdt$d
  
  Sx = svdx$d
  Ux = svdx$u
  Ux_t = t(svdx$v)
  
  # carry out conditional mean calculation for each 
  # column of y_star, where each col of y_star is a single y* as in Table .01
  # if we do it this way, we can vectorize calculations as much as possible
  # except for this one for loop....
  lambda = diag(1 / (sig2*kronecker(Sx, St) + nu))
  
  cond_mean_mat = matrix(nrow = nrow(R_xstar_x)*length(St), ncol = ncol(y_tilde))
  for (i in 1:ncol(y_tilde)){
    y_star_samp = y_tilde[, i]
  # reshape the response to be k * p
  zreshape = matrix(y_star_samp, nrow = s, ncol = nc)
  
  res1 = lambda %*% as.vector(Ut_t %*% zreshape %*% Ux)
  
  resmat = matrix(res1, nrow = s, ncol = nc)
  
  cond_mean_samp = sig2 * as.vector(Rt %*% Ut %*% resmat %*% Ux_t %*% t(R_xstar_x) )
  cond_mean_mat[,i] = cond_mean_samp
  }
  return(cond_mean_mat)
}

########################################################################################################
s = 11         # size of output space
nc = 10        # number of computer model runs

t1 = seq(0,1,length = s)
x1 = seq(0,1,length = nc)
x = expand.grid(t1, x1)

# function to emulate
f = (x[, 2]+1)*cos(pi*x[, 1]) + .03*(exp(x[, 2]))


persp(t1, x1, 
        matrix(f, nrow = s),
        theta = 130-90, 
        phi = 10,
        xlab = 't', ylab ='x', zlab='f',
        zlim = c(-2.2, 2.4)) -> res

points(trans3d(x[, 1], x[, 2], 
               f,
               pmat = res),
       col = 'black',
       pch = 16,
       cex = 0.7)

params = ML_pcEMU(fn = f, x = x1, q = 2)

s = 11  # size of output space

q = 2
xnew = c(.13, .72, .97)
nc_new = length(xnew)
x_grid = expand.grid(t1, x1)
f = (x_grid[, 2]+1)*cos(pi*x_grid[, 1]) + .03*(exp(x_grid[, 2]))
# should be 33 w new points

# f =  function output at train points
# xnew = vector of new locations
# nc = number of computer runs
# q = number of bases
# x1 = original x train locations
# params = object you get from the basis estimation function

w_stuff = get_wstar_distr_preds(f=f,xnew = xnew, nc=nc,q=2, x1=x1, t1=t1, params=params)

w1 = as.vector( w_stuff$wlist[[1]]$wmean)
w2 = as.vector( w_stuff$wlist[[2]]$wmean)
w1cov = w_stuff$wlist[[1]]$wcov
w2cov = w_stuff$wlist[[2]]$wcov

K = w_stuff$K
meanf = w_stuff$meanf

w_all = rbind(w1, w2)

basis_preds = K %*% w_all
basis_preds = as.vector(basis_preds + meanf)



persp(t1, x1,
      matrix(f, nrow = s),
      theta = 130-90,
      phi = 10,
      xlab = 't', ylab = 'x', zlab = 'f', 
      zlim = c(-2.4, 2.4)) -> res

points(trans3d(x[, 1], x[, 2], 
               f,
               pmat = res),
       col = 'black',
       pch = 16,
       cex = 0.7)

points(trans3d(plot_grid[, 1],
                 plot_grid[, 2], 
                 basis_preds,
                 pmat = res),
         col = 'red',
         pch = 16,
         cex = 0.7)

fmat0 + meanf
# this is what we should be getting!!
cond_mean[test_idx_kron]


w1 = as.vector( w_stuff$wlist[[1]]$wmean)
w2 = as.vector( w_stuff$wlist[[2]]$wmean)
w1cov = w_stuff$wlist[[1]]$wcov
w2cov = w_stuff$wlist[[2]]$wcov

K = w_stuff$K
meanf = w_stuff$meanf

w_all = rbind(w1, w2)

basis_preds = K %*% w_all
basis_preds = as.vector(basis_preds + meanf)

w1reals = rmultnorm(200, mu = w1, sigma = w1cov)
w2reals = rmultnorm(200, mu = w2, sigma = w2cov)

reals_matrix = matrix(nrow = 200, ncol = length(basis_preds))
for (i in 1:200){
  w_all = rbind(w1reals[i,], w2reals[i, ])
  bpreds = as.vector(K %*% w_all + meanf)
  reals_matrix[i, ] = bpreds
}


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


#######################################################################################
# Old stuff
#######################################################################################
nreals = 200

# can be a little faster if we do this
# set up--this is based on rmultnorm_SVD in Cov_functions.R
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

# a function to generate a matrix of N(0,1) variates and 
# compute L %*% Z %*% R as in (2)
pre_process_svd = function(z, p, k, C_U, C_d, R_U, R_d){
  zmat = matrix(z, nrow = k, ncol= p)
  res = as.vector((C_U %*% C_d) %*% zmat %*% t(R_U %*% R_d) )
}

# generate a vector of N(0,1) variates
allZs = matrix(rnorm(p*k*200), nrow=p*k)

# here we draw from Sigma11
# Sigma11 can be constructed as kronecker(Covx, Covt)
epsilon = t(rmultnorm(nreals, rep(0, nrow(Sigma11)), diag(sigma_sq, nrow = nrow(Sigma11))))

# generate u as a matrix
umat = apply(allZs, 2, function(x) pre_process_svd(z = x, p = p, k = k, 
                                                   C_U = C_U, C_d = C_d, R_U = R_U, R_d = R_d))

y_star = kmat %*% umat + epsilon

ustar = umat - Sigma12 %*% solve(Sigma11) %*% y_star 

conditional_samps = as.vector(zhat) + ustar

###############################
# APPENDIX
###############################
# this is Nychka's Sigma
Sigma22 = margvar*kronecker(R_xstar_xstar, Rt)

# construct K
kmat = matrix(0, nrow = length(train_idx_kron), ncol = nrow(Sigma22))
dim(kmat)
# b/c the training samples are stacked first, K will be identiy matrix for
# the first (train_idx_kron, train_idx_kron) row/cols and then just zeros
kmat[train_idx_kron, train_idx_kron] = diag(1, nrow= length(train_idx_kron))


# Note: this is kronecker(Covx, Covt)
Sigma11 = kmat %*% Sigma22 %*% t(kmat) + diag(sigma_sq,
                                              nrow = nrow(Rx)*nrow(Rt),
                                              ncol = nrow(Rx)*nrow(Rt))

Sigma12 = Sigma22 %*% t(kmat)

Sigma21 = kmat %*% Sigma22

# Nychka's zhat = conditional mean
zhat = cond_mean

# number of realizations
nreals = 200

# make the actual conditional covariance matrix
True_cond_cov =  Sigma22 - Sigma12 %*% solve(Sigma11) %*% Sigma21

# True realization from the conditional cov matrix
true_realization = as.vector(rmultnorm(1, zhat, True_cond_cov))


persp(t1, x1, 
      matrix(f, nrow = s),
      theta = 130-90, 
      phi = 10,
      xlab = 't', ylab ='x', zlab = 'f',
      zlim = c(-2.2, 2.4)) -> res

points(trans3d(x[, 1], x[, 2], 
               f,
               pmat = res),
       col = 'black',
       pch = 16,
       cex = 0.7)


points(trans3d(plot_grid[,1], plot_grid[,2],
               true_realization[test_idx_kron],
               pmat = res),
       col = 'red',
       pch = 16,
       cex = 0.7)

legend("topright", legend = c("f", "realization"), col = c("black", "red"), pch=19)

# for comparison to Nychka
Lots_of_realizations = rmultnorm(nreals, zhat, True_cond_cov)
est_mean_true = apply(Lots_of_realizations, 2, mean)
est_sds_true = apply(Lots_of_realizations, 2, sd)


epsilon = t(rmultnorm(nreals, rep(0, nrow(Sigma11)), diag(sigma_sq, nrow = nrow(Sigma11))))
u = rmultnorm(200, mu = rep(0, nrow(Sigma22)), Sigma22)
y_star = kmat %*% t(u) + epsilon
ustar = t(u) - Sigma12 %*% solve(Sigma11) %*% y_star 
conditional_samps = as.vector(zhat) + ustar
#test = t(apply(ustar, 2, function(x) x + zhat))
#est_mean = apply(test, 2, mean)
#est_sds = apply(test, 2, sd)


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