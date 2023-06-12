# show how a GP can be fit kronecker style

# make the 2-d x support
n1 = 11; n2 = 10        # good for plots
#n1 = 11*4; n2 = 10*4    # good for timing
t1=seq(0,1,length=n1)
x1=seq(0,1,length=n2)
x=expand.grid(t1,x1)

# a deterministic function to be emulated
f = (x[,2]+1)*cos(pi*x[,1])
f = (x[,2]+1)*cos(pi*x[,1]) + .03*(exp(x[,2]))

 # make a pretty fig for the chapter
PDF=FALSE
if(PDF) pdf('simplef.pdf',width=3,height=3)
par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(1,1,1,1))
persp(t1,x1,matrix(f,nrow=n1),theta = 40, phi = 10,xlab='t',ylab='x',zlab='f',zlim=c(-2.2,2.4)) -> res
points(trans3d(x[,1], x[,2], f, pmat = res), col = 'black', pch = 16,cex=.7)
if(PDF) dev.off()

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
Covtf = covpow(cbind(t1f,0),scale=1)
Covx = covpow(cbind(0,x1),scale=1)
CovKr = kronecker(Covx,Covt)
 # check that these covariance matrices are equal
summary(as.vector((CovFull-CovKr)))

 # use K-L decomp to come up with bases
eigt = eigen(Covt,symmetric=TRUE)

par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(4,4,1,1))
par(pty='s')
image(t1,t1,Covt)
par(pty='m')
plot(eigt$values^.5)
matplot(t1,eigt$vectors[,1:6],type='l')
 # grab the first q bases; n2 = number of multivariate simulations n_c
q=5
 # create basis matrix K; scale columns by sqrt(values) and length
 # This gives w-hat's that are ||k_i||^2 = 1.
K = eigt$vectors[,1:q] %*% diag(sqrt(eigt$values[1:q])) * sqrt(n1)
matplot(t1,K,type='l')
 # a slightly more efficient version
what = lsfit(K,fmat,int=FALSE)$coef
matplot(x1,t(what),type='p',pch=16)

 # note the crude variance (around 0) of each col of what not necessarily
 # 1 as it is in the svd basis case.  So we probably want to estimate
 # marginal variances for each w_k(x).
var0 <- function(x) mean(x^2)
apply(t(what),2,var0)

 # check: what could also be determined by projecting the sims onto K
In2 = diag(rep(1,n2))
Kbig = cbind(kronecker(In2,K[,1]),kronecker(In2,K[,2]),kronecker(In2,K[,3]),kronecker(In2,K[,4]),kronecker(In2,K[,5]))
# compute the w-hat's using lsfit - first with formula in chap
what2 = lsfit(Kbig,f,int=FALSE)$coef
 # check that this matches the other w-hat calucation
matplot(x1,matrix(what2,nrow=n2),type='p',pch=16)

 # computing the det needed in the likelihood calulation
 # det(t(Kbig)%*%Kbig) = det(t(K)%*%K)^n2
sum(log(diag(chol(t(Kbig)%*%Kbig))))  # slow, using Kbig
sum(log(diag(chol(t(K)%*%K))))*n2     # fast, requires orthogonal kj's

# make likelihood function
# find ML parameter settings
# make some "posterior" draws of the emulator at a new x*

