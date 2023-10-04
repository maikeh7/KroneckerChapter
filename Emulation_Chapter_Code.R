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
params = nmkb(c(0.5, .5, .0001, .5), 
                Estimate_params, distC = Xdist, distR = Tdist, z=f,
                lower = c(0.1, 0.1, 1e-6, .1), upper = c(100, 100, 100, 100))
# scale for Rx
print(params$par[1])

# scale for Rt
print(params$par[2])

# nugget nu
print(params$par[3])

#  marginal variance: sigma2
print(params$par[4])
scaleX = params$par[1]
scalet = params$par[2]

Rt = exp( -(Tdist / scalet)^2 )
Rx = exp( -(Xdist / scaleX)^2 )

s = 11

# prediction grid of new locations, xnew
xnew = c(.135, .72, .97)

# make the augmented x grid
xstar_grid = c(x1, xnew)

# build covariance matrix, R_x*,x*
dist_xstar_xstar = get_distmat(xstar_grid, xstar_grid)
R_xstar_xstar = exp( -(dist_xstar_xstar / scaleX)^2 )

# build cross covariance matrix, R_x*,x
dist_xstar_x = get_distmat(xstar_grid, x1)

R_xstar_x = exp( -(dist_xstar_x / scaleX)^2 ) 

svdx = svd(Rx)
svdt = svd(Rt)

Ux = svdx$u
Sx = svdx$d
Ux_t = t(svdx$v)

Ut = svdt$u
St = svdt$d
Ut_t = t(svdt$v)

  
nu = params$par[3] 
sig2 = params$par[4]

# calculate kronecker product of singular values Sx and St
lambda = diag(1 / (sig2*kronecker(Sx, St) + nu))

# reshape the response to be s * nc
f_reshape = matrix(f, nrow = s, ncol = nc)

res1 = lambda %*% as.vector(Ut_t %*% f_reshape %*% Ux)

resmat = matrix(res1, nrow = s, ncol = nc)

cond_mean = sig2 * as.vector(Rt %*% Ut %*% resmat %*% Ux_t %*% t(R_xstar_x) )

# grid for plotting 
# train locations in black; test locations in red
  
  plot_grid_train = expand.grid(t1, x1)
plot_grid = expand.grid(t1, xnew)

# this is for differentiating train/test indices
train_idx = 1:length(x1)
train_idx_kron = 1:(length(train_idx) * length(t1))
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
legend("bottomleft", legend = c("train", "test"), col = c("black", "red"), pch = 19)

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

pre_process_svd = function(z, s, Ut, St, Rstst_U, Rststd){
    nc_prime = length(Rststd)
    # need to take square roots!!
    eigs = kronecker(sqrt(Rststd), sqrt(St))
    ztilde = eigs*z
    mysim = as.vector(Ut %*% matrix(ztilde, nrow = s, ncol = nc_prime) %*% t(Rstst_U))
}

  

nreals = 100
nc_prime = length(Rststd)

# generate a matrix of N(0,1)'s
allZs = matrix(rnorm(s*nc_prime*nreals), nrow=s*nc_prime)

  
kmat = matrix(0, nrow = length(train_idx_kron), ncol = s*nc_prime)
kmat[train_idx_kron, train_idx_kron] = diag(1, nrow= length(train_idx_kron))

  
get_cond_mean = function(svdt, svdx, R_xstar_x, y_tilde, nu, sig2){
    
    Ut = svdt$u
    Ut_t = t(svdt$v)
    St = svdt$d
    
    Sx = svdx$d
    Ux = svdx$u
    Ux_t = t(svdx$v)
    
    # carry out conditional mean calculation for each 
    # column of y_tilde, where each col of y_tilde is a single y_tilde as in Table .01
    # if we do it this way, we can vectorize calculation steps as much as possible
    # except for this one for loop....
    lambda = diag(1 / (sig2*kronecker(Sx, St) + nu))
    
    cond_mean_mat = matrix(nrow = nrow(R_xstar_x)*length(St), ncol = ncol(y_tilde))
    for (i in 1:ncol(y_tilde)){
      y_star_samp = y_tilde[, i]
      # reshape to be s * nc
      zreshape = matrix(y_star_samp, nrow = s, ncol = nc)
      
      res1 = lambda %*% as.vector(Ut_t %*% zreshape %*% Ux)
      
      resmat = matrix(res1, nrow = s, ncol = nc)
      
      cond_mean_samp = sig2 * as.vector(Rt %*% Ut %*% resmat %*% Ux_t %*% t(R_xstar_x) )
      cond_mean_mat[,i] = cond_mean_samp
    }
    return(cond_mean_mat)
  }

# generate epsilon ~ MVN(0, nu)
epsilon = t(rmultnorm(nreals, rep(0, s*nc), diag(nu, nrow = s*nc)))

# use the pre_process_svd() function to get u and collect into a matrix
umat = apply(allZs,2, function(x) pre_process_svd(z=x, s=s, 
                                                  Ut = Ut, St = St, Rstst_U = Rstst_U, Rststd = Rststd))
# get y_tilde
y_tilde = kmat %*% umat + epsilon
# get w_star
w_star = get_cond_mean(svdt, svdx, R_xstar_x, y_tilde, nu, sig2)
ustar = umat - w_star
# we calculated cond_mean previously
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

# carry out optimization
  params = ML_pcEMU(fn = f, x = x1, q = 2)
print(params)

  
  # f =  function output at train points
  # xnew = vector of new locations
  # nc = number of computer runs
  # q = number of bases
  # x1 = original x train locations
  # params = object you get from the basis estimation function
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

  w_stuff = get_wstar_distr_preds(f=f,xnew = xnew, nc=nc,q=2, x1=x1, t1=t1, params=params)
# get predictive mean of w_1 and w_2
w1 = as.vector( w_stuff$wlist[[1]]$wmean)
w2 = as.vector( w_stuff$wlist[[2]]$wmean)

# grab the K matrix
K = w_stuff$K

# mean of f (f is demeaned prior to parameter estimation)
meanf = w_stuff$meanf

# make w's into a q x nc matrix 
w_all = rbind(w1, w2)

# get predictions, add mean back 
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

# get covariance matrices for w's
w1cov = w_stuff$wlist[[1]]$wcov
w2cov = w_stuff$wlist[[2]]$wcov

w1reals = rmultnorm(200, mu = w1, sigma = w1cov)
w2reals = rmultnorm(200, mu = w2, sigma = w2cov)

reals_matrix = matrix(nrow = 200, ncol = length(basis_preds))
for (i in 1:200){
  w_all = rbind(w1reals[i,], w2reals[i, ])
  bpreds = as.vector(K %*% w_all + meanf)
  reals_matrix[i, ] = bpreds
}

