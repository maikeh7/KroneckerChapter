library(raster)
# Dave's function to make a periodic distance matrix
# will be periodic over daysInYear
make_temporal_distmat = function(daysInYear){
  halfyear = round(daysInYear/2)
  t1=1:daysInYear
  t2=1:daysInYear
  t1c <- t1%%daysInYear
  t2c <- t2%%daysInYear
  dist1 <- outer(t1c,t2c,'-')
  #
  t1cs <- (t1+halfyear)%%daysInYear
  t2cs <- (t2+halfyear)%%daysInYear
  dist2 <- outer(t1cs,t2cs,'-')
  #
  dist3 <- ifelse(abs(dist1)<abs(dist2),dist1, dist2)
  dist3 = abs(dist3)
  return(dist3)
}

##########################################################
# -------  Dave's rmultnorm function
#Generates variates from multivariate normal distribution
#with mean vector mu and covariance matrix sigma
#n is the number of desired realizations
#########################################################
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


# function to generate MVN realizations (using svd decomposition, not cholesky)
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


svdinv <- function(a,tol=1e-6){
  sva <- svd(a)
  maxsv <- max(sva$d)
  dd <- ifelse(sva$d<tol*maxsv,0,1/sva$d)
  t(sva$u %*% (t(sva$v)*dd))
}

# Make covariance matrix of pred locs (Grid) to training data (N)
# Gaussian cov function
MakeSigGridN = function(sampleLocs, finegrid, pow = 2, scale=.4) {
  distStation2fineGrid = matrix(NA, nrow = nrow(sampleLocs), ncol = nrow(finegrid))
  for (i in 1:nrow(sampleLocs)){
    x1 = as.numeric(sampleLocs[i, 1])
    y1= as.numeric(sampleLocs[i, 2])
    for (j in 1:nrow(finegrid)){
      x2 = finegrid[j, 1]
      y2= finegrid[j, 2]
      distStation2fineGrid[i,j] = sqrt(((x1-x2)^2 + (y1-y2)^2))
    }
  }
  
  SigmaGrid_N=t(distStation2fineGrid)
  
  SigmaGrid_N = exp(-(SigmaGrid_N/scale)^pow)
  return(SigmaGrid_N)
}

#SigGN = MakeSigGridN(sampleLocs = locs, finegrid = dat)
#SigGN = Sigtest
#dim(SigGN)
# Gaussian covariance function
covpow <- function(locs,pow=2,scale=5){
  # browser()
  # create distance matrix from observation data matrix
  d1 <- dist(locs)
  # get an n value from matrix dimensions
  n <- dim(locs)[1]
  # create nXn empty matrix
  mat <- matrix(0,n,n)
  # indexing the distance matrix to empty nXn matrix as a lower triangular?
  
  #Maike: yeah, the dist() function creates a dist object, not a matrix, and we need a matrix, so this is a way to do that
  mat[lower.tri(mat)] <- d1
  
  # make mat the sum of mat and mat transpose so that it is symmetric which is a covariance property
  # Maike: yes!
  mat <- mat+t(mat)
  # take the exponential of negative (mat/scale)^pow which is a guassian covariance matrix?
  # Gaussian covariance function results in a valid covariance matrix --> C() = exp(-(d_{i,j} / scale)^2),
  
  cc <- exp(-(mat/scale)^pow)
  # makes sense
  return(cc)
}

#make temporal covariance matrix
# lengthscale parameter
sc = 0.6
# make prediction grid
xfine = seq(0, 1, by = .05)

yfine =  seq(0, 1, by = .05)

dat = expand.grid(xfine, yfine)
nrow(dat)

locs = dat

# training data (take some samples)
numsamps = 150
sampindx = sample(1:nrow(dat), numsamps)
locs = dat[sampindx, ]
head(locs)
locs = cbind(locs$Var1, locs$Var2)

# make SigmaGridN (GN) using training data and prediction grid
SigGN = MakeSigGridN(sampleLocs = locs, finegrid = dat, pow = 2, scale = sc)

# marginal variance
margvar = .5
# nugget for spatial dimension
nugget =.001
# nugget for the temporal dimension
temporal_nugget = 0
# smoothness parameter for periodic cov matrix
ell = .3

# this is the temporal distance matrix (it will be periodic)
numdays=5
mydistC = make_temporal_distmat(daysInYear = numdays)

# this is the temporal covariance matrix. It is periodic 
# ell controls smoothness
SigmaT = exp(-2*ell^2 * (sin(pi * mydistC/numdays))^2 ) +
  diag(temporal_nugget, nrow = nrow(mydistC), ncol = ncol(mydistC))
dim(SigmaT)
z <- rmultnorm(1,rep(0,nrow(SigmaT)), SigmaT)
plot(as.vector(z), type = "l")

#SigmaGG = make_SigmaGG(gridsize = .05, scale = sc) 
SigmaN = covpow(locs, pow = 2, scale = sc) 

p = dim(SigmaT)[1]
k = dim(SigmaN)[1]

KronCovmat = margvar * kronecker(SigmaT, SigmaN) + diag(nugget, nrow = p*k, ncol = p*k) 
mySim = rmultnorm(1, rep(0, p*k), KronCovmat)
dim(mySim)
day1 = mySim[1:numsamps]
plot(locs[,1], locs[,2])
spdf = data.frame(x = locs[,1], y = locs[,2], z = as.vector(mySim))
coordinates(spdf) = ~x+y
spplot(spdf, "z")

# this is brute force post mean (you don't want to actually calculate this!)
postmeantruth = as.vector(margvar*kronecker(SigmaT, SigGN) %*% solve(KronCovmat) %*% t(mySim))
for (i in 1:numdays){
  mystart = start
  end = mystart + (dim(SigGN)[1]-1)
  day1 = postmeantruth[mystart:end]
  dmat = matrix(day1, nrow = 21, ncol = 21)
  rmat = raster(dmat)
  plot(rmat, main = paste("cond mean, day ", i))
  points(locs[,1], locs[,2], pch=19)
  #print(mean(day1))
  start = end + 1
  #hist(day1)
  
  #dev.copy(png, paste("TEST", i, ".png", sep = ""))
  #ev.off()
}
#####################################################################################################
# marginal variance
margvar = 2
# nugget for spatial dimension
nugget =.01
# nugget for the temporal dimension
temporal_nugget = 0.0001
# smoothness parameter for periodic cov matrix
ell = .3
# this is the temporal distance matrix (it will be periodic)
numdays=5
# scale 
sc = 0.6
# make prediction grid
xfine = seq(0, 1, by = .05)

yfine =  seq(0, 1, by = .05)

dat = expand.grid(xfine, yfine)
nrow(dat)

locs = dat

# training data (take some samples)
numsamps = 150
sampindx = sample(1:nrow(dat), numsamps)
locs = dat[sampindx, ]
head(locs)
locs = cbind(locs$Var1, locs$Var2)

SigmaN = covpow(locs, pow = 2, scale = sc) 
# make SigmaGridN (GN) using training data and prediction grid
SigGN = MakeSigGridN(sampleLocs = locs, finegrid = dat, pow = 2, scale = sc)

mydistC = make_temporal_distmat(daysInYear = numdays)

# this is the temporal covariance matrix. It is periodic 
# ell controls smoothness
SigmaT = exp(-2*ell^2 * (sin(pi * mydistC/numdays))^2 ) +
  diag(temporal_nugget, nrow = nrow(mydistC), ncol = ncol(mydistC))

C = SigmaT
R = SigmaN
G = SigGN

#svd C and R
svdC = svd(C)
svdR = svd(R)

Uc = svdC$u
Sc = svdC$d
Uc_t = t(svdC$v)

Ur = svdR$u
Sr = svdR$d
Ur_t = t(svdR$v)

# This is the very slow way of getting post conditional mean
p=nrow(C)
k=nrow(R)
KronCovmat = margvar * kronecker(SigmaT, SigmaN) + diag(nugget, nrow = p*k, ncol = p*k) 
mySim = rmultnorm(1, rep(0, p*k), KronCovmat)
postmeantruth = as.vector(margvar*kronecker(SigmaT, SigGN) %*% solve(KronCovmat) %*% t(mySim))

#inverse of kronecker of singular values of C and R
InvSSmat = diag(1/(margvar*as.vector(outer(Sr, Sc)) + nugget))

#reshape z to be p x k
zreshape = matrix(mySim, nrow=k, ncol=p)
dim(zreshape)

#res1 = InvSSmat %*% as.vector(solve(Ur) %*% zreshape %*% solve(t(Uc)))
res1 = InvSSmat %*% as.vector(solve(Ur) %*% zreshape %*% Uc)


resmat = matrix(res1, nrow = k, ncol = p)

PostMeanTest = margvar * as.vector(G%*%solve(Ur_t) %*% resmat %*% t(C%*%solve(Uc_t)))
testing = margvar * as.vector(G %*% Ur %*% resmat %*% Uc_t %*% t(C) ) # equivalent to above but simpler

head(postmeantruth)
head(PostMeanTest)
################################################################################################
# if we do sigmaT kronecker sigmaN instead, like Dave does in his example
################################################################################################
margvar=1
KronCovmat = margvar * kronecker(SigmaN, SigmaT) + diag(nugget, nrow = p*k, ncol = p*k) 
mySim = rmultnorm(1, rep(0, p*k), KronCovmat)
postmeantruth = as.vector(margvar*kronecker(SigGN, SigmaT) %*% solve(KronCovmat) %*% t(mySim))
dim(kronecker(SigGN, SigmaT))
dim(KronCovmat)
#inverse of kronecker of singular values of C and R
InvSSmat = diag(1/(margvar*as.vector(outer(Sc, Sr)) + nugget))

#reshape z to be p x k
zreshape = matrix(mySim, nrow=p, ncol=k)

#750x750 x vector(150x150 %*% 150x5 %*% 5x5)
res1 = InvSSmat %*% as.vector(solve(Uc) %*% zreshape %*% Ur)

resmat = matrix(res1, nrow = p, ncol = k)

result = margvar * as.vector(C %*% Uc %*% resmat %*% Ur_t %*% t(G) ) 

