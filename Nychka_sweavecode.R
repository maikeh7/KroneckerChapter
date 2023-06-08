source("Cov_functions.R")
source("Estimation_1D.R")
# make the 2-d x support
nt = 11; nx = 10     

t1 = seq(0,1,length=nt)
x1 = seq(0,1,length=nx)
x = expand.grid(t1,x1)

# a deterministic function to be emulated
f = (x[,2]+1)*cos(pi*x[,1]) + rnorm(nrow(x), sd = 1e-6)

Tdist = make_distmat(cbind(t1,0))
Xdist = make_distmat(cbind(0, x1))

# this function assumes the kronecker product: kronecker(C, R), so order matters
params = nmkb(c(0.5, .5, .0001, .5), 
              Estimate_params, distC = Xdist, distR = Tdist, z=f,
              lower = c(0.1, 0.1, 1e-6, .1), upper = c(100, 100, 100, 100))

print(params$par[1])
print(params$par[2])
print(params$par[3])
print(params$par[4])

scaleX = params$par[1]
scalet = params$par[2]
Covt = covpow(cbind(t1, 0), scale = scalet)
Covx = covpow(cbind(0, x1), scale = scaleX)

nt = 11; nstar = 25

# prediction grid for x^*
xstar = seq(0,1, length = nstar)

# in Nychka's setup, combine train and test points
xstar_grid = c(x1, xstar)

# This will be used for getting Nychka's Sigma 
# (Sigma = kronecker(Cov_(x*,x*), Cov_t))
Cov_xstar_xstar = covpow(cbind(0, xstar_grid), scale = scaleX)

# make distance matrices (used for obtaining conditional mean)
distmat_xstar_x = abs(outer(xstar_grid, x1, "-"))

# build cross covariance matrix
Cov_xstar_x = exp(-(distmat_xstar_x/scaleX)^2) 

# eigen decomposition of Covx = C and
# Covt = R and Cov_xstar_x = G
R = Covt
C = Covx
G = Cov_xstar_x

svdC = svd(C)
svdR = svd(R)
svdG = svd(G)

Uc = svdC$u
Sc = svdC$d
Uc_t = t(svdC$v)

Ur = svdR$u
Sr = svdR$d
Ur_t = t(svdR$v)

p = nrow(C) 
k = nrow(R) 
Cnug = params$par[3]
margvar = params$par[4]

# inverse of kronecker of singular values of C and R
# equivalent to kronecker(Sc, Sr)! 
InvSSmat = diag(1/(margvar*as.vector(outer(Sr, Sc)) + Cnug))

# reshape the response to be k * p
zreshape = matrix(f, nrow = k, ncol = p)

res1 = InvSSmat %*% as.vector(solve(Ur) %*% zreshape %*% Uc)

resmat = matrix(res1, nrow = k, ncol = p)

cond_mean = margvar * as.vector(R %*% Ur %*% resmat %*% Uc_t %*% t(G) )

# grid for plotting 
# note we only plot test locations here
plot_grid = expand.grid(t1, xstar)

persp(t1, xstar,
      matrix(cond_mean[test_idx_kron], nrow = nt),
      theta = 130-90,
      phi = 10,
      xlab = 't', ylab = 'x', zlab = 'f', 
      zlim = c(-3,3)) -> res

points(trans3d(plot_grid[,1],
               plot_grid[,2], 
               cond_mean[test_idx_kron],
               pmat = res),
       col = 'black',
       pch = 16,
       cex = 0.7)


# this is Nychka's Sigma
Sigma22 = kronecker(Cov_xstar_xstar, Covt)

train_idx = 1:length(x1)
train_idx_kron = 1:(length(train_idx) * length(t1))
test_idx_kron = (train_idx_kron[length(train_idx_kron)] + 1): nrow(Sigma22)

# construct K
kmat = matrix(0, nrow = length(train_idx_kron), ncol = nrow(Sigma22))
dim(kmat)
# b/c the training samples are stacked first, K will be identiy matrix for
# the first (train_idx_kron, train_idx_kron) row/cols and then just zeros
kmat[train_idx_kron, train_idx_kron] = diag(1, nrow= length(train_idx_kron))

# Note: this is kronecker(Covx, Covt)
Sigma11 = kmat %*% Sigma22 %*% t(kmat) + diag(Cnug,
                                              nrow = nrow(Covx)*nrow(Covt),
                                              ncol = nrow(Covx)*nrow(Covt))

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

persp(t1, xstar, matrix(true_realization[test_idx_kron], nrow = nt),theta = 130-90, phi = 10,
      xlab='t',ylab='x',zlab='f',zlim=c(-5,5)) -> res
points(trans3d(plot_grid[,1], plot_grid[,2], true_realization[test_idx_kron], pmat = res),
       col = 'black', pch = 16,cex=.7)

# for comparison to Nychka
Lots_of_realizations = rmultnorm(200, zhat, True_cond_cov)
est_mean_true =apply(Lots_of_realizations, 2, mean)
est_sds_true = apply(Lots_of_realizations, 2, sd)

# vectorized, still kinda slow
epsilon = t(rmultnorm(nreals, rep(0, nrow(Sigma11)), diag(sigmasq, nrow = nrow(Sigma11))))
u = rmultnorm(nreals, mu = rep(0, nrow(Sigma22)), Sigma22)
y_star = kmat %*% t(u) + epsilon
ustar = t(u) - Sigma12 %*% solve(Sigma11) %*% y_star 
conditional_samps = as.vector(zhat) + ustar
test = t(apply(ustar, 2, function(x) x + zhat))
est_mean = apply(test, 2, mean)
est_sds = apply(test, 2, sd)

# can be a little faster if we do this
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
umat = apply(allZs,2, function(x) pre_process_svd(z=x, p=p, k=k, C_U = C_U, C_d = C_d, R_U = R_U, R_d = R_d))
y_star = kmat %*% umat + epsilon
ustar = umat - Sigma12 %*% solve(Sigma11) %*% y_star 
# this is much faster
conditional_samps = as.vector(zhat) + ustar

