library(matrixcalc)
library(MixMatrix)
?MixMatrix::dmatrixnorm
library(Matrix)
library(raster)
#function to construct distance matrix
doyDist <- function(t1,t2=tvals,daysInYear=365){
  # computes the distance in days between the days in the vector t1
  # and the days in the vector t2, using only the day of the year for this.
  # So Dec 30 and Jan 1 are only 2 days apart, even if they're in different
  # years.
  
  halfyear = round(daysInYear/2)
  t1c <- t1%%daysInYear
  t2c <- t2%%daysInYear
  dist1 <- outer(t1c,t2c,'-')
  #
  t1cs <- (t1+halfyear)%%daysInYear
  t2cs <- (t2+halfyear)%%daysInYear
  dist2 <- outer(t1cs,t2cs,'-')
  #
  dist3 <- ifelse(abs(dist1)<abs(dist2),dist1,dist2)
}

##########################################################
# -------  Dave's rmultnorm function
#Generates variates from multivariate normal distribution
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

covpow <- function(locs,pow=2,scale=5){
  # browser()
  d1 <- dist(locs)
  n <- dim(locs)[1]
  mat <- matrix(0,n,n)
  mat[lower.tri(mat)] <- d1
  mat <- mat+t(mat)
  cc <- exp(-(mat/scale)^pow)
  return(cc)
}
#ok, I think this is OK, could adjust scale to adjust strength of temporal autocorrelation.
#Y <- rmultnorm(1,rep(0,365),sigmaT)
#plot(1:365,Y,type = "l") #looks reasonable probably


locs = finegrid
SigmaGG = covpow(locs, pow = 2, scale = .25)
SigmaGGNear = nearPD(SigmaGG)
SigmaGG2 = SigmaGGNear$mat
SigmaGG2 = as.matrix(SigmaGG2)
is.positive.semi.definite(SigmaGG2)

##simulate from z ###########################################
#############################################################
#make spatio temporal matrix
SigGGT = kronecker(sigmaT, SigmaGG2) 


################################################################################
##conditional distribution
###############################################################################
#pick some locations--let's do 78
#############################################################################
testdays2 = rmultnorm(1, rep(0, nrow(SigGGT)), SigGGT)

samplestations = sample(1:nrow(finegrid), size = 18)
sampleLocs = finegrid[samplestations, ]

sampleLocs$Index = samplestations
start = 1
muvec =vector(mode = "numeric")
for (i in 1:nrow(sigmaT)){
  mystart = start
  end = mystart + (nrow(finegrid)-1) #no of grid points - 1
  d1 = testdays2[mystart:end]
  start = end + 1
  z <- matrix(d1, ncol=sqrt(nrow(SigmaGG)))
  obs = z[samplestations]
  muvec = c(muvec, obs)
}

sigmaT  = makeSigmaT(daysInYear = 11, scale = 5)
SigmaGrid_N = MakeSigGridN(sampleLocs, finegrid = scale = .25) #sigmaGrid,N
Zn = muvec
SigmaNN = covpow(sampleLocs[,1:2], pow = 2, scale = .25) #get sigmaN,N

#brute force method
dim(sigmaT)
BFresult = kronecker(sigmaT, SigmaGrid_N) %*% solve(kronecker(sigmaT, SigmaNN)) %*% Zn

#faster result
solveSigNN = solve(SigmaNN)

FastResult = SigmaGrid_N %*% solveSigNN %*% matrix(Zn, nrow = nrow(sampleLocs), ncol = ncol(sigmaT))
head(as.vector(FastResult))
dim(FastResult)
head(BFresult)
length(zn)




Zreshape = matrix(Zn, nrow = length(samplestations), ncol = 11)
dim(Zreshape)

product2 = SigmaGrid_N %*% solveSigNN
dim(product2)

test = as.vector(product2 %*% Zreshape)
head(test)
head(BFresult)




#sigma grid, n
m=4
n=3
p=2
sigGridN = matrix(rnorm(m*n), nrow = m, ncol = n) #grid points to stations
sigmann =matrix(rnorm(n*n), nrow = n, ncol = n) #stations to stations
sigmaT = matrix(runif(p*p), nrow = p, ncol = p) #temporal cov matrix
Zn = rnorm(n*p) 

sigGridN = SigmaGrid_N
sigmaT = sigmaT
sigmann = SigmaNN
p=11
n=18
Zn = muvec
#but we still need to invert...
Zreshape = matrix(Zn, nrow = n, ncol = p)
BFsolution = kronecker(sigmaT, sigGridN) %*% solve(kronecker(sigmaT, sigmann)) %*% Zn #solution is dim m x p
BFsolution
fastSolution = sigGridN %*% solve(sigmann) %*% Zreshape 

###################################################################################
#some useful functions
#####################################################################################
MakeSigGridN = function(sampleLocs, finegrid, pow = 2, scale) {
distStation2fineGrid = matrix(NA, nrow = nrow(sampleLocs), ncol = nrow(finegrid))
for (i in 1:nrow(sampleLocs)){
  x1 = as.numeric(sampleLocs[i, 1])
  y1= as.numeric(sampleLocs[i, 2])
  for (j in 1:nrow(finegrid)){
    x2 = finegrid$Var1[j]
    y2= finegrid$Var2[j]
    distStation2fineGrid[i,j] = sqrt(((x1-x2)^2 + (y1-y2)^2))
  }
}

SigmaGrid_N=t(distStation2fineGrid)

SigmaGrid_N = exp(-(SigmaGrid_N/scale)^pow)
return(SigmaGrid_N)
}


#daysInYear is however many days you want to model
makeSigmaT = function(daysInYear=11, pow = 2, scale){
  
  #daysInYear = 11
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
  dist3 <- ifelse(abs(dist1)<abs(dist2),dist1,dist2)
  dist3 = abs(dist3)
  mat = dist3

  sigttest <- exp(-(mat/scale)^pow)
  pdtest = nearPD(x = sigttest, keepDiag = T)
  pdtest = pdtest$mat
  sigmaT = as.matrix(pdtest)
  return(sigmaT)
}

CalcPostMean = function(SigmaGrid_N, sigmaNN, Zn, sigmaT){
  solveSigNN = solve(sigmaNN)
  Zrows = dim(SigmaGrid_N)[2]
  Zcols = dim(sigmaT)[1]
  Zreshape = matrix(Zn, nrow = Zrows, ncol = Zcols)
  CondGausMu = SigmaGrid_N %*% solveSigNN %*% Zreshape
  return(as.vector(CondGausMu))
  
}

BruteForceCondmean = function(SigmaGrid_N, SigmaNN, Zn, sigmaT){
  BFsolution= kronecker(sigmaT, sigGridN) %*% solve(kronecker(sigmaT, sigmann)) %*% Zn 
  return(BFsolution)
}

##alternative rmultnorm function when covariance matrix is kronecker product 
#finish this
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

#junk
distStation2fineGrid = matrix(NA, nrow = nrow(sampleLocs), ncol = nrow(finegrid))
for (i in 1:nrow(sampleLocs)){
  x1 = as.numeric(sampleLocs[i, 1])
  y1= as.numeric(sampleLocs[i, 2])
  for (j in 1:nrow(finegrid)){
    x2 = finegrid$Var1[j]
    y2= finegrid$Var2[j]
    distStation2fineGrid[i,j] = sqrt(((x1-x2)^2 + (y1-y2)^2))
  }
}

SigmaGrid_N=t(distStation2fineGrid)
SigmaGrid_N = exp(-(SigmaGrid_N/scale)^pow)


La = chol(matA, pivot=T)
Lb = chol(matB, pivot=T)
La = chol(matA)
Lb = chol(matB)

kronABChol = kronecker((t(La) %*% La) , (t(Lb) %*% Lb))
kronABChol = kronecker(CholFactorA %*% LaT, CholFactorB %*% LbT)

kronABChol = (kronecker(Lb, La)) %*% t(kronecker(Lb, La))
t(kronABChol)

pivB = attr(Lb, "pivot")
pivA = attr(La, "pivot")
#this is equivalent to kronab!!!!!!!!!!!!!!!!!!!!!!!!
kronecker(t(La),t(Lb)) %*% t(kronecker(t(La), t(Lb)))
n=2
m=3
matA = genPositiveDefMat("eigen",dim=n)
matA = matA$Sigma
matB = genPositiveDefMat("eigen",dim=m)
matB = matB$Sigma
kronAB = kronecker(matA, matB)
La = chol(matA, pivot = T)
Lb = chol(matB, pivot = T)
chol(kronAB, pivot = T)
kronecker(La, Lb)


t1 = rmultnormKRON(matA, matB, mu=rep(0,n*m), n=1)
t2 = rmultnorm(n=1, mu=rep(0,n*m), sigma = kronAB)
sort(t1)
sort(t2)

rmultnormKRON = function(matA, matB, mu, n=1){
  p = length(mu)
  Lb = chol(matB, pivot=T)
  La = chol(matA, pivot =T)
  pivA = attr(La, "pivot")
 print(pivA)
  pivB = attr(Lb, "pivot")
  print(pivB)
  ordA = order(pivA)
  ordB = order(pivB)
 
  L = kronecker(La, Lb)
  print(L)
  set.seed(44839)
  z = matrix(rnorm(n * p), nrow=n)
  
  zz = La %*% z %*% Lb
  realization = zz + matrix(mu, nrow = n, ncol = p, byrow = T)
  return(realization)
}
p


sigma=kronAB
rmultnorm <- function(n, mu, sigma){
  # returns n rows of multinormal mu,sigma random vectors
  # ie: returns a n x length(mu) matrix
  p <- length(mu)
  set.seed(44839)
  z <- matrix(rnorm(n * p), nrow=n)

  ch <- chol(sigma, pivot=T)
  print(ch)
  piv <- attr(ch, "pivot")
  print(piv)
  zz <- (z%*%ch)
  
  zzz <- 0*zz
  zzz[, piv] <- zz
  zzz + matrix(mu, nrow=n, ncol=p, byrow=T)
}
kronUV = kronecker(myU,myV)
sigma = kronUV
mysvd = svd(sigma)

SVD_U = mysvd$u
Droot = diag(sqrt(mysvd$d))
set.seed(44839)
z <- matrix(rnorm(2*400),ncol=1)
h = SVD_U %*% Droot %*% z
h

start = Sys.time()
testdays2 = rmultnorm(1, rep(0, nrow(SigGGT)), SigGGT)
end = Sys.time()
end - start


#junk

piv <- attr(ch, "pivot")
zz <- (z%*%ch)
zzz <- 0*zz
zzz[, piv] <- zz
TRUTH = zzz + matrix(mu, nrow=1, ncol=p, byrow=T)
plot(as.vector(TRUTH), type = "l")
.85*2.72


L = (kronecker(Lb, La))
zz <- (z%*%L)
# zzz <- 0*zz
# zzz[, c(2*pivA, 2*pivA-1)] <- zz
# pivB
# L = (kronecker(Lb, La))
# kronecker(La,Lb)
# dim(L)
# zz <- (z%*%L)
TEST = zz + matrix(mu, nrow=1, ncol=p, byrow=T)
TEST
TRUTH
SigGGT = kronecker(sigmaT, SigmaGG2) 
dim(SigGGT)
sigma = SigGGT




n=2
m=400
matA = genPositiveDefMat("eigen",dim=n)
matA = matA$Sigma
matB = genPositiveDefMat("eigen",dim=m)
matB = matB$Sigma


U <- matA
V <- matB
set.seed(44839)
z <- matrix(rnorm(n*m),nrow=nrow(U))

mu = matrix(rep(0,n*m),nrow=nrow(U))

Lu = chol(U)
Lu = chol(U, pivot = T) 
pivU = attr(Lu, "pivot")
ordU = order(pivU)
print(ordU)
Lu = Lu[,ordU]

Lv = chol(V)
Lv = chol(V, pivot = T)
pivV = attr(Lv, "pivot")
ordV = order(pivV)
print(ordV)
Lv = Lv[,ordV]

result = Lv %*% z %*% Lu + mu
sort(as.vector(t(z %*% Lu)  %*% Lv))
dim(Lv)

result
sort(result)
sort(truth)
kronUV = kronecker(U,V)
sigma = kronUV
p <- n*m
set.seed(44839)
z <- matrix(rnorm(1 * p), nrow=1)
ch=chol(sigma)
ch <- chol(sigma, pivot=T)

piv <- attr(ch, "pivot")
print(piv)
zz <- z %*% ch
zz
zzz <- 0*zz
zzz[, piv] <- zz
truth = zzz + matrix(rep(0,n*m), nrow=1, ncol=p, byrow=T)
sort(result)
sort(zz)
U
V
library(MixMatrix)


myU = U
myV = V

set.seed(44839)
testing <- rmatrixnorm(n = 1, mean=matrix(rep(0, n*m), nrow = 8), U=myU, V=myV)

mu + t(chol(U)) %*% y %*% chol(V)




set.seed(44839)
z <- matrix(rnorm(n*m),nrow=nrow(myU))
mu = matrix(rep(0,64),nrow=8)
testingagain = mu + t(chol(myU)) %*% z %*% chol(myV)

kronUV = kronecker(myU,myV)
sigma = kronUV
mysvd = svd(sigma)

SVD_U = mysvd$u
Droot = diag(sqrt(mysvd$d))
set.seed(44839)
z <- matrix(rnorm(2*400),ncol=1)
h = SVD_U %*% Droot %*% z
h
zz1 = zz[1:400]
h1 = h[1:400]
rz= raster(z)


z = matrix(h1, ncol = 20)

persp(z, theta = -60, phi=30) #ok, looks good
rz = raster(z)
plot(rz)

U
V
n*m

Utest <- 5 * diag(3) + 1
Vtest <- matrix(c(2,0,0,.1),nrow=2)
mu = matrix(1:6,nrow=3)
set.seed(20180203)
z <- rmatrixnorm(n = 1, mean=mu,U=Utest,V=Vtest)
Utest

###############################################################################

SigGGT = kronecker(sigmaT, SigmaGG2) 

rmultnorm <- function(n, mu, sigma){
  # returns n rows of multinormal mu,sigma random vectors
  # ie: returns a n x length(mu) matrix
  p <- length(mu)
  z <- matrix(rnorm(n * p), nrow=n)
  ch <- chol(sigma, pivot=T, tol = 1e-6)
  piv <- attr(ch, "pivot")
  zz <- (z%*%ch)
  zzz <- 0*zz
  zzz[, piv] <- zz
  zzz + matrix(mu, nrow=n, ncol=p, byrow=T)
}

?chol
testdays2=rmultnorm(1, mu = rep(0, nrow(SigGGT)), sigma = SigGGT )
d1 = testdays2[1:441]

z = matrix(testdays2, ncol = sqrt(nrow(SigmaGG2)))
persp(z, theta = -60, phi=30) 
rz= raster(z)
persp(z, theta = -60, phi=30) #ok, looks good
plot(rz)
dmatrixnorm
A <- rmatrixnorm(n=10,mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
                 L=matrix(c(2,1,0,.1),nrow=2),list=TRUE)
A[[1]]
# demonstrating the dmatrixnorm function
dmatrixnorm(A[[1]],mean=matrix(c(100,0,-100,0,25,-1000),nrow=2),
            L=matrix(c(2,1,0,.1),nrow=2),log=TRUE )

mykron = kronecker(cbm2, cga2)
mykron[mykron>0]
dim(mykron)
z <- rmultnorm(1,rep(0,25*25),mykron)
zmat = matrix(z, nrow=25, ncol = 25)
mymean=matrix(0, nrow=25, ncol = 25)
nearcbm = nearPD(cbm,keepDiag = T)
nearcga = nearPD(cga, keepDiag = T)
cbm2 = as.matrix(nearcbm$mat)
cga2 = as.matrix(nearcga$mat)
chol(cbm2, tol = 1e-100)
?chol
image(cbm2)
is.positive.semi.definite(cbm2)
dmatrixnorm(zmat, mean = mymean, U=cbm2, V=cga2, log=T)
cbm
