# Read in the cosmology design 5x37
x <- matrix(scan('design.txt'),ncol=5,byrow=T)
p <- ncol(x); nc = nrow(x)
deslabs <- c('n','h',expression(sigma[8]),expression(Omega[D]),
             expression(Omega[B]))
pairs(x,labels=deslabs)

# read in the dark matter simulations
y <- matrix(scan('smoothspectra.txt'),byrow=T,ncol=223)
logkNative <- y[,1]
 # grab the simulations corresponding to a redshift of z=0
ycNative = y[,seq(7,by=6,length=37)]
matplot(logkNative,ycNative,type='l',lty=1,col='gray',
        xlab='', ylab='')
mtext(expression(log[10](frac(Delta,k^{1.5}))),side=2,line=1.8)
mtext(expression(log[10](k)),side=1,at=-.6,line=2.4)
mtext('wavenumber',side=1,at=-1.3,line=2.4)

 # interpolate evenly over logk
logk = seq(min(logkNative),max(logkNative),length=400)
 # to use the language of the book chapter
# support t of the simulations
t1 = logk
# matrix of simulation output, put on logk support
interpLogk <- function(y,x=logkNative,xout=logk) approx(x,y,xout)$y
yc = apply(ycNative,2,interpLogk)

# design x which is p = 5 parameters X nc = 37 runs
# standardize the design so that it resides on [0,1]^5

to01 <- function(x) (x-min(x))/diff(range(x))
x01 = apply(x,2,to01)
pairs(x01,labels=deslabs)

 # plot balanced data
matplot(logk,yc,type='l',lty=1,col='gray',
        xlab='', ylab='')
mtext(expression(log[10](frac(Delta,k^{1.5}))),side=2,line=1.8)
mtext(expression(log[10](k)),side=1,at=-.6,line=2.4)
mtext('wavenumber',side=1,at=-1.3,line=2.4)

write(t(x01),'c3design01-37x5.txt',ncol=5)
write(t(yc),'c3simsz0-300x37.txt',ncol=37)
write(logk,'c3logk-300.txt')


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

#---- stop here -----#
fmat = matrix(f,nrow=n1)
fsvd = svd(fmat)

par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(4,4,1,1))
matplot(t1,fmat,type='l')
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
what2 = solve(t(Kbig)%*%Kbig,t(Kbig)%*%f)
solve(t(K)%*%K,t(K)%*%fmat)
#what2 = solve(t(Kbig)%*%Kbig+diag(rep(1e-10,q*n2)),t(Kbig)%*%f)
# note, since this system is so singular (1 singular value explains everything)
# what2 blows up for the second component.  The svd-based what calculation above is more stable.

 # det(t(Kbig)%*%Kbig) = det(t(K)%*%K)^n2
sum(log(diag(chol(t(Kbig)%*%Kbig))))
sum(log(diag(chol(t(K)%*%K))))*n2

# make likelihood function
# find ML parameter settings
# make some "posterior" draws of the emulator at a new x*






