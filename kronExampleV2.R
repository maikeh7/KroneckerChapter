# show how a GP can be fit kronecker style

# make the 2-d x support
n1 = 11; n2 = 10        # good for plots
#n1 = 11*4; n2 = 10*4    # good for timing
n1 = 33; n2 = 30 
t1=seq(0,1,length=n1)
x1=seq(0,1,length=n2)
x=expand.grid(t1,x1)

# a deterministic function to be emulated
f = (x[,2]+1)*cos(pi*x[,1])

# look at function over a 2-d grid
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(1,1,1,1))
persp(t1,x1,matrix(f,nrow=n1),theta = 130-90, phi = 10,xlab='t',ylab='x',zlab='f',zlim=c(-2,2)) -> res
points(trans3d(x[,1], x[,2], f, pmat = res), col = 'black', pch = 16,cex=.7)

 # make a pretty fig for the chapter
PDF=FALSE
if(PDF) pdf('simplef.pdf',width=3,height=3)
par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(1,1,1,1))
persp(t1,x1,matrix(f,nrow=n1),theta = 40, phi = 10,xlab='t',ylab='x',zlab='f',zlim=c(-2.2,2.4)) -> res
points(trans3d(x[,1], x[,2], f, pmat = res), col = 'black', pch = 16,cex=.7)
if(PDF) dev.off()

###### let's use this example to show the knronecker approach
 # load is some spatial cov functions
#source('~/Dave/Courses/Spatial/Shortcourse/SampleCode/gensurf.r')
 # x1 and x2

# a function to make the Gaussian covariance matrix
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

 # produce the standard (n1*n2)^2 covariance matrix
CovFull = covpow(x,scale=1)

 # do the same thing with a kronecker trick
Covt = covpow(cbind(t1,0),scale=1)
Covx = covpow(cbind(0,x1),scale=1)
CovKr = kronecker(Covx,Covt)
 # check that these covariance matrices are equal
summary(as.vector((CovFull-CovKr)))

 # now the half quadratic summation - the long way
# sqrt{  f^T * Cov^-1 * f }
Cnug = diag(rep(1,nrow(CovFull)))*1e-6
t(f)%*%solve(CovFull+Cnug,f)

sum(solve(t(chol(CovFull+Cnug)),f)^2)
sum(forwardsolve(t(chol(CovFull+Cnug)),f)^2)

 # same thing with pivoting in the chol call
chCovFull = chol(CovFull+Cnug,pivot=T)
ipiv = attr(chCovFull,"pivot")
sum(forwardsolve(t(chCovFull),f[ipiv])^2)

# same solve using the eigen decomposition
C1Eigen = eigen(CovFull+Cnug,symmetric=TRUE)
 # look at first few eigen functions
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(1,1,1,1))
for(k in 1:4){
  persp(t1,x1,matrix(C1Eigen$vectors[,k]*sqrt((C1Eigen$values[k])),nrow=length(t1)),theta = -50, phi = 10,
        xlab='t',ylab='x',zlab='f',zlim=c(-1.1,1.1)/sqrt(1)) -> res
  mtext(paste('eigenfunction',k),side=3,line=-.1)
}

C1svd = svd(CovFull)  #svd and Eigen look about the same
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(1,1,1,1))
for(k in 1:4){
  persp(t1,x1,matrix(C1svd$u[,k]*sqrt((C1svd$d[k])),nrow=length(t1)),theta = -50, phi = 10,
        xlab='t',ylab='x',zlab='f',zlim=c(-1.1,1.1)/sqrt(1)) -> res
  mtext(paste('eof',k),side=3,line=-.1)
}

# f^T UD^{-1}U^T f - the half quadratic
# sqrt{  z^T * Cov^-1 * z } = sqrt(sum( solve(Cov^.5 * z)^2 ))
sum( (f %*% C1Eigen$vectors %*% diag(C1Eigen$values^(-.5)))^2 )
 # again, but with eigen applied to CovFull
C1Eigen0 = eigen(CovFull,symmetric=TRUE)
diagEV = C1Eigen0$values + diag(Cnug)
sum( (f %*% C1Eigen0$vectors %*% diag(diagEV^(-.5)) )^2 )
 # or a bit more computationally efficient by doing the first mult, and then a dot.
sum( ( as.vector(f %*% C1Eigen0$vectors)*diagEV^(-.5) )^2 )

## Next, the same solves with kroneckered eigen decompositions ##
 # recall in our construction of x1 and x2, x1 is varying faster so
fmat = matrix(f,nrow=n1)
persp(t1,x1,fmat,theta = 120, phi = 10, xlab='x1',ylab='x2',zlab='f',zlim=c(-2.3,2.3))
CovtE = eigen(Covt,symmetric=TRUE)
CovxE = eigen(Covx,symmetric=TRUE)

diagEVkron = as.vector(outer(CovtE$values,CovxE$values))+diag(Cnug)
 # kron version 1
sum( (f %*% kronecker(CovxE$vectors,CovtE$vectors) %*% diag(diagEVkron^(-.5)) )^2 )
 # now, do the solves separately

 # this is as far as I've gotten;  I'm hoping you can help out now!
# check 1: is eigen(C1) kron eigen(C2) = eigen(C1 kron C2)?
CovtE = eigen(Covt,symmetric=TRUE)
CovxE = eigen(Covx,symmetric=TRUE)
ec1 = CovxE$vectors %*% diag(CovxE$values^.5)
ec2 = CovtE$vectors %*% diag(CovtE$values^.5)
eKron = kronecker(ec1,ec2)
eKron
# 
eCov = C1Eigen0$vectors %*% diag(C1Eigen0$values^.5) #get some Nans
head(eCov)
#  # eigenvalues become negative or NAs with this poorly conditioned matrix
#  # Instead, try svd
# 
Cov2 = Covt
Cov1 = Covx
sCovFull = svd(CovFull); Cov2S = svd(Cov2); Cov1S = svd(Cov1)
sCov = sCovFull$u %*% diag(sqrt(sCovFull$d))
sc1 = Cov1S$u %*% diag(sqrt(Cov1S$d))
sc2 = Cov2S$u %*% diag(sqrt(Cov2S$d))
# 
image(sCov)
image(kronecker(sc1,sc2))
image(kronecker(sc2,sc1))

browser()

 # code to make speedier kronecker calculations
C = CovFull; C1 = Covx; C2 = Covt
Csvd = svd(C); C1svd = svd(C1); C2svd = svd(C2)
U = Csvd$u; d = Csvd$d
U1 = C1svd$u; d1 = C1svd$d
U2 = C2svd$u; d2 = C2svd$d
nugg = rep(1e-6,nrow(C))

 # compute some quadratic forms - first, the long way
f%*%kronecker(U1,U2)%*%solve(kronecker(diag(d1),diag(d2))+ diag(nugg))%*%t(kronecker(U1,U2))%*%f
t(f)%*%solve(C+diag(nugg),f)            

 # fast(er) computation of t(f)%*%solve(C+diag(nugg),f)
f12 = solve(U2, matrix(f,nrow=nrow(C2)))  # recall inv(U2) = t(U2)
f12 = t(U2)%*%matrix(f,nrow=nrow(C2))

f22 = solve(U1,matrix(t(f12),nrow=nrow(C1)))
f22 = t(U1)%*%matrix(t(f12),nrow=nrow(C1))


fkr = as.vector(t(f22))
svvec = kronecker(d1,d2)+nugg
fkr2 = fkr / sqrt(svvec)
sum(fkr2^2)


 # some timings
library(tictoc)
tic(); for(k in 1:10) t(f)%*%solve(C+diag(nugg),f); toc() # 9.5 sec

tic(); for(k in 1:10){                                     # 0.018 sec
 C1svd = svd(C1); C2svd = svd(C2)
  U = Csvd$u; d = Csvd$d
  U1 = C1svd$u; d1 = C1svd$d
  U2 = C2svd$u; d2 = C2svd$d
f12 = t(U2)%*%matrix(f,nrow=nrow(C2))
f22 = t(U1)%*%matrix(t(f12),nrow=nrow(C1))
fkr = as.vector(t(f22))
svvec = kronecker(d1,d2)+nugg
fkr2 = fkr / sqrt(svvec)
sum(fkr2^2)
}; toc()

##############################################################################
## Maike added this -- result for post cond mean
# Dave, you can check
##############################################################################

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


n1 = 11; n2 = 10        # good for plots
#n1 = 11*4; n2 = 10*4    # good for timing
t1=seq(0,1,length=n1)
x1=seq(0,1,length=n2)
x=expand.grid(t1,x1)
Covt = covpow(cbind(t1,0),scale=1)
Covx = covpow(cbind(0,x1),scale=1)
CovKr = kronecker(Covx,Covt)

pow=2
scale=1
nstar = 15
xstar = seq(0,1, length = nstar)
distmat = abs(outer(xstar, x1, "-"))
distmat2 = abs(outer(x1, xstar, "-"))
distmat3 = abs(outer(xstar, xstar, "-"))
# build cross cov mat
Cov_xstar_x = exp(-(distmat/scale)^pow) #nstar x n1
Cov_x_xstar = exp(-(distmat2/scale)^pow)

Cov_xstar_xstar = exp(- (distmat3/scale)^pow)
# rename stuff otherwise I get confused
C = Covt
R = Covx
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
p=nrow(C)
k=nrow(R)

# get conditional mean the fast way
margvar=1
# but first the slow way so we can check!

Cnug = 1e-6
KronCovmat = margvar * kronecker(Covx, Covt) + diag(Cnug, nrow = p*k, ncol = p*k) 
mySim = rmultnorm(1, rep(0, p*k), KronCovmat)

# this is the slow way
postmeantruth = as.vector(margvar*kronecker(Cov_xstar_x, Covt) %*% solve(KronCovmat) %*% t(mySim))

#inverse of kronecker of singular values of C and R
InvSSmat = diag(1/(margvar*as.vector(outer(Sc, Sr)) + Cnug))

#reshape z to be p x k
zreshape = matrix(mySim, nrow=p, ncol=k)

res1 = InvSSmat %*% as.vector(solve(Uc) %*% zreshape %*% Ur)

resmat = matrix(res1, nrow = p, ncol = k)

result = margvar * as.vector(C %*% Uc %*% resmat %*% Ur_t %*% t(G) ) 
# check if correct
all.equal(postmeantruth, result)

# deal w/ covariance
C = Covx
R = Covt

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

Sig11 = kronecker(Covx, Covt) + diag(Cnug, nrow = p*k, ncol = p*k) 
Sig12 = kronecker(Cov_x_xstar, Covt)
Sig21 = kronecker(Cov_xstar_x, Covt)
Sig22 = kronecker(Cov_xstar_xstar, Covt)
True_cond_cov = Sig22 - Sig21 %*% solve(Sig11) %*% Sig12
truth = Sig21 %*% solve(Sig11) %*% Sig12
X = diag(1/(margvar*as.vector(outer(Sr, Sc)) + Cnug))

Sig11_inv_truth = solve(Sig11)

Sig11_inv = kronecker(Uc, Ur) %*% X %*% solve(kronecker(Uc, Ur))

# this is not exactly the same as truth but I think that's bc of rounding error...hopefully
# I can't reduce this, so maybe we can work with this?
# can't cholesky decompose RHS or LHS! issue w/ positive definiteness
# this really isn't any faster
# we will have to simulate using Nychka's method 
test1 = Sig22 - Sig21 %*% Sig11_inv %*% kronecker(t(G), R)


