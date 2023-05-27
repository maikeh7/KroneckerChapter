# show how a GP can be fit kronecker style

# make the 2-d x support
n1 = 11; n2 = 10        # good for plots
#n1 = 11*4; n2 = 10*4    # good for timing
t1=seq(0,1,length=n1); t1f = seq(0,1,length=101)
x1=seq(0,1,length=n2)
x=expand.grid(t1,x1)

# a deterministic function to be emulated
f = (x[,2]+1)*cos(pi*x[,1])
f = (x[,2]+1)*cos(pi*x[,1]) + .03*(x[,2]-.5)^2 -.003056
f = (x[,2]+1)*cos(pi*x[,1]) + .03*(exp(x[,2]))

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
matplot(t1f,eigt$vectors[,1:6],type='l')
 # grab the first q bases; n2 = number of multivariate simulations n_c
q=5
K = eigt$vectors[,1:q] %*% diag(sqrt(eigt$values[1:q])) * sqrt(n1)
#K = eigt$vectors[,1:q] * sqrt(n1)
matplot(t1,K,type='l')
 # compute the w-hat's [can't: not svd on t support]
what = lsfit(Kbig,f,int=FALSE)$coef
what = lsfit(K,fmat,int=FALSE)$coef
matplot(x1,matrix(what2,nrow=n2),type='p',pch=16)
matplot(x1,t(what),type='p',pch=16)

 # note the crude variance (around 0) of each col of what is 1
 # So modeling each col of what with a GP(0,1*R(x,x')) isn't too bad.
var0 <- function(x) mean(x^2)
apply(what,2,var0)

 # a pdf for the paper
PDF=FALSE
if(PDF) pdf('klSimple.pdf',width=10.5,height=2.8)
par(mfrow=c(1,4),oma=c(0,0,0,0),mar=c(4,4,1.7,1))
par(pty='s')
image(t1,t1,Covt,xlab='',ylab='')
mtext('output support t',side=1,line=2.2)
mtext('output support',side=2,line=2.2)
mtext('correlation R(t)',side=3,line=.2)
par(pty='m')
plot(sqrt(eigt$values),pch=16,col=c(1:5,rep(8,10)),xlab='',ylab='')
mtext('basis number',side=1,line=2.2)
mtext(expression(sqrt(lambda[j])),side=2,line=2.2)
mtext('eigenvalues',side=3,line=.2)
matplot(t1,K,type='l',lty=1,xlab='',ylab='',ylim=c(-3.5,3.5))
matlines(t1,K[,c(2,4)],type='l',lty=1,lwd=2,col=c(2,4))
mtext('output support t',side=1,line=2.2)
mtext(expression(k[j]),side=2,line=2.2)
mtext('q=5 KL bases',side=3,line=.2)
matplot(x1,t(what),pch=16,xlab='',ylab='',ylim=c(-2,2))
mtext('input x',side=1,line=2.2)
mtext(expression(hat(w)[j](x)),side=2,line=2.1)
mtext('weights',side=3,line=.2)
if(PDF) dev.off()

 # check: what could also be determined by projecting the sims onto K
In2 = diag(rep(1,n2))
Kbig = cbind(kronecker(In2,K[,1]),kronecker(In2,K[,2]),kronecker(In2,K[,3]),kronecker(In2,K[,4]),kronecker(In2,K[,5]))
#Kbig = kronecker(In2,K[,1])

what2 = solve(t(Kbig)%*%Kbig,t(Kbig)%*%f)
solve(t(K)%*%K,t(K)%*%fmat)
#what2 = solve(t(Kbig)%*%Kbig+diag(rep(1e-10,q*n2)),t(Kbig)%*%f)
# note, since this system is so singular (1 singular value explains everything)
# what2 blows up for the second component.  The svd-based what calculation above is more stable.

 # det(t(Kbig)%*%Kbig) = det(t(K)%*%K)^n2
sum(log(diag(chol(t(Kbig)%*%Kbig))))  # slow, using Kbig
sum(log(diag(chol(t(K)%*%K))))*n2     # fast, requires orthogonal kj's

# make likelihood function
# find ML parameter settings
# make some "posterior" draws of the emulator at a new x*






