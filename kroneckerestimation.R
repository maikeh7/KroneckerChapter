source("Cov_functions.R")
source("Estimation_1D.R")
library(dfoptim)
############
# Example 1
############
# Estimating params using our deterministic fucntion f
nt = 11; nx = 10
t1=seq(0,1,length=nt)
x1=seq(0,1,length=nx)
x=expand.grid(t1,x1)


# a deterministic function to be emulated
f = (x[,2]+1)*cos(pi*x[,1]) + rnorm(nrow(x), sd= 0.0001)

# C = x
# R = t
Tdist = make_distmat(cbind(t1,0))
Xdist = make_distmat(cbind(0, x1))
persp(t1,x1,matrix(f,nrow=nt),theta = 130-90, phi = 10,xlab='t',ylab='x',zlab='f',zlim=c(-3,3)) -> res
points(trans3d(x[,1], x[,2], f, pmat = res), col = 'black', pch = 16,cex=.7)

testmkb = nmkb(c(0.5, .5, .0001, .5), Estimate_params, distC = Xdist, distR = Tdist, z=f,
               lower = c(0.1, 0.1, 1e-8, .1), upper = c(100, 100, 100, 100))
testmkb
# HMMMMMM.....not sure these make much sense....but conditional mean below looks ok.
params= testmkb$par
# scaleC           ScaleR      Nugget      Marginal var
# 100.00000000   0.56668151   0.00000001  99.99999787
scaleX = params[1]
scalet = params[2]
Cnug = params[3]
margvar = params[4]
pow=2

Covt = covpow(cbind(t1,0),scale= scalet)
Covx = covpow(cbind(0,x1),scale= scaleX)
CovKr = kronecker(Covx,Covt) # so we can check results later

nt = 11; nstar = 25

# prediction grid for x
xstar = seq(0,1, length = nstar)
xstar_grid = expand.grid(t1, xstar)
# make distance matrices
distmat_xstar_x = abs(outer(xstar, x1, "-"))
distmat_x_xstar = abs(outer(x1, xstar, "-"))
distmat_xstar_xstar = abs(outer(xstar, xstar, "-"))

# build cross covariance matrices
Cov_xstar_x = exp(-(distmat_xstar_x/scaleX)^pow) 
Cov_x_xstar = exp(-(distmat_x_xstar/scaleX)^pow)
Cov_xstar_xstar = exp(- (distmat_xstar_xstar/scaleX)^pow)

# svd/eigen decomps of covx = C and Covt = R
R = Covt
C = Covx
G = Cov_xstar_x

#svd C and R
svdC = svd(C)
svdR = svd(R)
svdG = svd(G)

Uc = svdC$u
Sc = svdC$d
Uc_t = t(svdC$v)

Ur = svdR$u
Sr = svdR$d
Ur_t = t(svdR$v)

# This is the very slow way of getting post conditional mean
p = nrow(C) 
k = nrow(R) 

KronCovmat = margvar * kronecker(Covx, Covt) + diag(Cnug, nrow = p*k, ncol = p*k) 
#mySim = rmultnorm(1, rep(0, p*k), KronCovmat)

# this is the slow way
postmeantruth = as.vector(margvar*kronecker(Cov_xstar_x, Covt) %*% solve(KronCovmat) %*% f)

# inverse of kronecker of singular values of C and R
# equivalent to kronecker(Sc, Sr)! but faster obviously
InvSSmat = diag(1/(margvar*as.vector(outer(Sr, Sc)) + Cnug))

# reshape z to be kxp
zreshape = matrix(f, nrow=k, ncol=p)

#res1 = InvSSmat %*% as.vector(solve(Uc) %*% zreshape %*% Ur)
res1 = InvSSmat %*% as.vector(solve(Ur) %*% zreshape %*% Uc)

#resmat = matrix(res1, nrow = p, ncol = k)
resmat = matrix(res1, nrow=k, ncol=p)

#result = margvar * as.vector(C %*% Uc %*% resmat %*% Ur_t %*% t(G) ) 
result = margvar * as.vector(R %*% Ur %*% resmat %*% Uc_t %*% t(G) )

# inspect resulting predicted mean surface
persp(t1,xstar,matrix(result, nrow = nt),theta = 130-90, phi = 10,xlab='t',ylab='x',zlab='f',zlim=c(-3,3)) -> res
points(trans3d(xstar_grid[,1], xstar_grid[,2], result, pmat = res), col = 'black', pch = 16,cex=.7)

length(result)
# check if correct
# these are NOT 100% identical. Not sure why. Rounding error maybe? 
all.equal(postmeantruth, result)
postmeantruth

#################
# Example 2
#################
# actually I don't think we need this...
# Make up some data that comes from the true generating process 

#######################################################################################################
