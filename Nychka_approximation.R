source("Cov_functions.R")
n1 = 11; n2 = 10       # good for plots
#n1 = 11*4; n2 = 10*4    # good for timing
t1 = seq(0,1,length=n1)
x1 = seq(0,1,length=n2)
x = expand.grid(t1,x1)
Covt = covpow(cbind(t1,0),scale=1)
Covx = covpow(cbind(0,x1),scale=1)
Cov_x_t = kronecker(Covx, Covt)

pow = 2
scale = 1
nstar = 20
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
# b/c the training samples are stacked first, kmat will be identiy matrix for
# the first several row/cols and then just zeros
kmat[train_idx_kron, train_idx_kron] = diag(1, nrow= length(train_idx_kron))
# get the other components

sigmasq = 0.001

#checked--this is kronecker(covx, covt)
Sigma11 = kmat %*% Sigma22 %*% t(kmat) + diag(sigmasq,
                    nrow = nrow(Covx)*nrow(Covt), ncol = nrow(Covx)*nrow(Covt))
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
persp(t1, xstar, matrix(true_realization[test_idx_kron], nrow = n1),theta = 130-90, phi = 10,
      xlab='t',ylab='x',zlab='f',zlim=c(-3,3)) -> res
points(trans3d(plot_grid[,1], plot_grid[,2], true_realization[test_idx_kron], pmat = res),
       col = 'black', pch = 16,cex=.7)


Lots_of_realizations = rmultnorm(200, zhat, True_cond_cov)
dim(Lots_of_realizations)
est_mean_true =apply(Lots_of_realizations, 2, mean)
est_sds_true = apply(Lots_of_realizations, 2, sd)

##############
# vectorized #
##############
epsilon = t(rmultnorm(200, rep(0, nrow(Sigma11)), diag(sigmasq, nrow = nrow(Sigma11))))
u = rmultnorm(200, mu = rep(0, nrow(Sigma22)), Sigma22)
y_star = kmat %*% t(u) + epsilon
ustar = t(u) - Sigma12 %*% solve(Sigma11) %*% y_star 
conditional_samps = zhat + ustar
test = t(apply(ustar, 2, function(x) x + zhat))
est_mean = apply(test, 2, mean)
est_sds = apply(test, 2, sd)
