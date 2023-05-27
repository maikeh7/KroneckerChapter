# show how a GP can be fit kronecker style

# make the 2-d x support
n1 = 11; n2 = 10        # good for plots
#n1 = 11*4; n2 = 10*4    # good for timing
t1=seq(0,1,length=n1)
x1=seq(0,1,length=n2)
x=expand.grid(t1,x1)

# a deterministic function to be emulated
#f = (x[,2]+1)*cos(pi*x[,1])
f = (x[,2]+1)*cos(pi*x[,1]) + .03*(exp(x[,2]))

# look at function over a 2-d grid
par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(1,1,1,1))
persp(t1,x1,matrix(f,nrow=n1),theta = 130-90, phi = 10,xlab='t',ylab='x',zlab='f',zlim=c(-2.2,2.4)) -> res
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
Covx = covpow(cbind(0,x1),scale=1)

# recall in our construction of x1 and x2, x1 is varying faster so
# make the matrix of computer model output where each column is a run.
fmat = matrix(f,nrow=n1)

 # subtract out the mean of the columns of fmat
meanf = apply(fmat,1,mean)
fmat0 = fmat - meanf
fsvd = svd(fmat0)

par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(4,4,1,1))
matplot(t1,fmat0,type='l')
plot(fsvd$d)
matplot(t1,fsvd$u,type='l')
 # grab the first two bases; n2 = number of multivariate simulations n_c
q=2
K = fsvd$u[,1:q] %*% diag(fsvd$d[1:q]) / sqrt(n2)
matplot(t1,K,type='l')
 # compute the w-hat's
what = fsvd$v[,1:q] * sqrt(n2)
matplot(x1,what,type='p',pch=16)
 # note the crude variance (around 0) of each col of what is 1
 # So modeling each col of what with a GP(0,1*R(x,x')) isn't too bad.
var0 <- function(x) mean(x^2)
apply(what,2,var0)

 # check: what could also be determined by projecting the sims onto K
In2 = diag(rep(1,n2))
Kbig = cbind(kronecker(In2,K[,1]),kronecker(In2,K[,2]))
#Kbig = kronecker(In2,K[,1])

what2 = solve(t(Kbig)%*%Kbig,t(Kbig)%*%as.vector(fmat0))
solve(t(K)%*%K,t(K)%*%fmat0)
#what2 = solve(t(Kbig)%*%Kbig+diag(rep(1e-10,q*n2)),t(Kbig)%*%f)
# note, since this system is so singular (1 singular value explains everything)
# what2 blows up for the second component.  The svd-based what calculation above is more stable.

 # det(t(Kbig)%*%Kbig) = det(t(K)%*%K)^n2
sum(log(diag(chol(t(Kbig)%*%Kbig))))  # slow, using Kbig
sum(log(diag(chol(t(K)%*%K))))*n2     # fast, requires orthogonal kj's

# make likelihood function
# find ML parameter settings
# make some "posterior" draws of the emulator at a new x*






