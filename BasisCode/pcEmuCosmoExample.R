# Read in the cosmology design 5x37
x <- matrix(scan('c3design01-37x5.txt'),ncol=5,byrow=T)
p <- ncol(x); nc = nrow(x)
deslabs <- c('n','h',expression(sigma[8]),expression(Omega[D]),
             expression(Omega[B]))
pairs(x,labels=deslabs)

# read in the dark matter simulations
y <- matrix(scan('c3simsz0-400x37.txt'),byrow=T,ncol=37)
logk <- scan('c3logk-400.txt')

# n_c number of model runs
nc = nrow(x)
s = length(logk)

 # plot the simulations
matplot(logk,y,type='l',lty=1,col='gray',
        xlab='', ylab='')
mtext(expression(log[10](frac(Delta,k^{1.5}))),side=2,line=1.8)
mtext(expression(log[10](k)),side=1,at=-.6,line=2.4)
mtext('wavenumber',side=1,at=-1.3,line=2.4)

# support t of the simulations
t1 = logk

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

 # subtract out the mean of the simulations - it can be added for predictions later
meany = apply(y,1,mean)
y0 = y - meany
#---- get pc/svd basis -----#
ysvd = svd(y0)

par(mfrow=c(2,2),oma=c(0,0,0,0),mar=c(4,4,1,1))
matplot(t1,y,type='l')
plot(ysvd$d)
matplot(t1,ysvd$u,type='l')
 # grab the first 5 bases; n2 = number of multivariate simulations n_c
q=5
K = ysvd$u[,1:q] %*% diag(ysvd$d[1:q]) / sqrt(nc)
matplot(t1,K,type='l')
 # compute the w-hat's
what = ysvd$v[,1:q] * sqrt(nc)
matplot(x[,4],what,type='p',pch=16)
 # note the crude variance (around 0) of each col of what is 1
 # So modeling each col of what with a GP(0,1*R(x,x')) isn't too bad.
var0 <- function(x) mean(x^2)
apply(what,2,var0)

 # check: what could also be determined by projecting the sims onto K
Inc = diag(rep(1,nc))
Kbig = cbind(kronecker(Inc,K[,1]),kronecker(Inc,K[,2]),kronecker(Inc,K[,3]),kronecker(Inc,K[,4]),kronecker(Inc,K[,5]))
what2 = solve(t(Kbig)%*%Kbig,t(Kbig)%*%as.vector(y))
solve(t(K)%*%K,t(K)%*%y)

 # det(t(Kbig)%*%Kbig) = det(t(K)%*%K)^n2
sum(log(diag(chol(t(Kbig)%*%Kbig))))
sum(log(diag(chol(t(K)%*%K))))*nc

# make likelihood function
# find ML parameter settings
# make some "posterior" draws of the emulator at a new x*






