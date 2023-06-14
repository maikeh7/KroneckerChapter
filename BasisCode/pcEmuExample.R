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
dim(fmat)
nc=ncol(fmat)
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
 # compute the w-hat's -- see page 10 of chapter
what = fsvd$v[,1:q] * sqrt(n2)
matplot(x1,what,type='p',pch=16)
dim(K)
dim(what) 
# note the crude variance (around 0) of each col of what is 1
 # So modeling each col of what with a GP(0,1*R(x,x')) isn't too bad.
var0 <- function(x) mean(x^2)
apply(what,2,var0)

 # check: what could also be determined by projecting the sims onto K
In2 = diag(rep(1,n2))
Kbig = cbind(kronecker(In2,K[,1]),kronecker(In2,K[,2]))
#Kbig = kronecker(In2,K[,1])
# this is the calculation on page 8
what2 = solve(t(Kbig)%*%Kbig,t(Kbig)%*%as.vector(fmat0))
solve(t(K)%*%K,t(K)%*%fmat0)
dim(Kbig)

fvec %*% (diag(1, nrow = 110, ncol = 110) - Kbig %*% solve(t(Kbig) %*% Kbig) %*% t(Kbig) ) %*% t(fvec)
fvec = matrix(f, nrow=1, ncol = 110) 

testing=t(fmat) %*% (diag(1, nrow = 11, ncol=11)- K %*% solve(t(K) %*% K) %*% t(K)) %*% fmat0 
sum(testing)
dim(fmat0)

#what2 = solve(t(Kbig)%*%Kbig+diag(rep(1e-10,q*n2)),t(Kbig)%*%f)
# note, since this system is so singular (1 singular value explains everything)
# what2 blows up for the second component. The svd-based what calculation above is more stable.

 # det(t(Kbig)%*%Kbig) = det(t(K)%*%K)^n2 # I don't get this--ask Dave! How are these equal??
sum(log(diag(chol(t(Kbig)%*%Kbig))))  # slow, using Kbig
sum(log(diag(chol(t(K)%*%K))))*n2     # fast, requires orthogonal kj's

# make likelihood function
# find ML parameter settings
# make some "posterior" draws of the emulator at a new x*
j=1
nrow(K)
covpow
length(f)
K %*% (solve)
# what is s???
# note n2 = nc
nc = n2
parameters = list(1, c(2,2), c(1,2))
Estimate_params = function(parameters, f, n1, x1){
  nu = parameters[1]
  phi = parameters[2:3]
  sig2 = parameters[4:5]
  # make f into matrix of appropriate dims
  fmat = matrix(f, ncol = n1)
  nc = ncol(fmat)
  # subtract means of rows
  meanf = apply(fmat,1,mean)
  fmat0 = fmat - meanf
  
  # do svd
  fsvd = svd(fmat0)
  
  # number of bases to include
  q = length(phi)
  n2 = nc
  n1 = nrow(fmat0)
  # make K
  K = fsvd$u[,1:q] %*% diag(fsvd$d[1:q]) / sqrt(n2)
  
  s = nrow(K)
  
  # also stupid!!! we don't want to do this!!! bad bad bad
  In2 = diag(rep(1,n2))
  Kbig = cbind(kronecker(In2,K[,1]),kronecker(In2,K[,2]))
  fvec = matrix(fmat0, nrow=1, ncol = length(f))
  
  # this is stupid but wtf can't quickly see another way to get this bad bad bad
  fcalc = fvec %*% (diag(1, nrow = length(f), ncol = length(f)) - 
                      Kbig %*% solve(t(Kbig) %*% Kbig) %*% t(Kbig) ) %*% t(fvec)
  
  # or nc?? not n2? i guess they are equal
  part_two = ((s-q)*n2 / 2) * log(nu) - 1/2 * 1 / nu * fcalc
  
  # get what...this would only work for svd though
  # actually, don't need this
  #what = fsvd$v[,1:q] * sqrt(n2)
  
  results_vec = vector(length = q)
  Xdist = make_distmat(cbind(0, x1))
  for(j in q){
    kj = K[, j]
    covmat_j  = sig2[j]*exp(-(phi[j]*Xdist)^2)
    part_one = nu / t(kj) %*% (kj) # this is a scalar
    # so...this does not make sense to me....can't add scalar to matrix in R....
    # ask Dave...not sure about this!!
    part_one = matrix(rep(part_one, nrow(covmat_j)*nrow(covmat_j)), nrow=nrow(covmat_j))
    chunk = part_one + covmat_j
    
    LogDet = 2*sum(log(diag(chol(chunk)))) # why is this off by a factor of 2?
    LogDet = determinant(chunk, logarithm = TRUE)$modulus
    what_j = fsvd$v[ ,q] * sqrt(n2)
    # faster
    ssq = sum(solve(t(chol(chunk)),what_j)^2)
    loglike = - .5*LogDet - .5*ssq - 0.5 * part_two
    
    results_vec[j] = loglike
  } 
  
  MylogLik = sum(results_vec)
  return(-MylogLik) #return negative log likelihood
}

params = nmkb(parameters = c(1,1,1,1,1), 
              Estimate_params, f=myF, n1=myn, x1=myx,
              lower = c(.01,.01,.01,.01,.01), upper = rep(20, 5))

lower = .01, upper = 100)
myF = f
myn = n1
myx = x1
params = nmkb(parameters = 1, 
              Estimate_params, f=myF, n1=myn, x1=myx),
              lower = .01, upper = 100)
dim(fmat)
nc
Estimate_params = function(parameters, f, n1, x1){
  nu =1
  phi = c(.001, 1)
  sig2 = c(2,2)
  # make f into matrix of appropriate dims
  fmat = matrix(f, ncol = n1)
  nc = ncol(fmat)
  # subtract means of rows
  meanf = apply(fmat,1,mean)
  fmat0 = fmat - meanf
  
  # do svd
  fsvd = svd(fmat0)
  
  # number of bases to include
  q = length(phi)
  n2 = nc
  n1 = nrow(fmat0)
  # make K
  K = fsvd$u[,1:q] %*% diag(fsvd$d[1:q]) / sqrt(n2)
  
  s = nrow(K)
  
  # also stupid!!! we don't want to do this!!! bad bad bad
  In2 = diag(rep(1,n2))
  Kbig = cbind(kronecker(In2,K[,1]),kronecker(In2,K[,2]))
  fvec = matrix(fmat0, nrow=1, ncol = length(f))
  
  # this is stupid but wtf can't quickly see another way to get this bad bad bad
  fcalc = fvec %*% (diag(1, nrow = length(f), ncol = length(f)) - 
                      Kbig %*% solve(t(Kbig) %*% Kbig) %*% t(Kbig) ) %*% t(fvec)
  
  # or nc?? not n2? i guess they are equal
  part_two = ((s-q)*n2 / 2) * log(nu) - 1/2 * 1 / nu * fcalc
  
  # get what...this would only work for svd though
  # actually, don't need this
  #what = fsvd$v[,1:q] * sqrt(n2)
  
  results_vec = vector(length = q)
  Xdist = make_distmat(cbind(0, x1))
  j=1
  for(j in q){
    kj = K[, j]
    covmat_j  = sig2[j]*exp(-(phi[j]*Xdist)^2)
    part_one = nu / t(kj) %*% (kj) # this is a scalar
    # so...this does not make sense to me....can't add scalar to matrix in R....
    # ask Dave...not sure about this!!
    part_one = matrix(rep(part_one, nrow(covmat_j)*nrow(covmat_j)), nrow=nrow(covmat_j))
    chunk = part_one + covmat_j
    
    #LogDet = 2*sum(log(diag(chol(chunk)))) # why is this off by a factor of 2?
    LogDet = determinant(chunk, logarithm = TRUE)$modulus
    what_j = fsvd$v[ ,q] * sqrt(n2)
    # faster
    
    cholChunk = chol(Ch+ diag(0.00001, nrow=nrow(chunk)))
    ssq = sum(solve(cholChunk, what_j^2))
    loglike = - .5*LogDet - .5*ssq - 0.5 * part_two
    
    results_vec[j] = loglike
  } 
  
  MylogLik = sum(results_vec)
  return(-MylogLik) #return negative log likelihood
}





q=1
dim(kj)
nu=10
kj=rnorm(11)
kj = rep(5, 11)
nrow(K)
distmat = make_distmat(cbind(0, kj))
covmat_j  = 2*exp(-(.9*distmat)^2)
part_one = nu / (kj) %*% t(kj)
chunk = part_one + covmat_j
what_j = fsvd$v[ ,q] * sqrt(n2)

t(what_j) %*% solve(chunk) %*% what_j
sum((solve(chol(chunk)) %*% what_j)^2)
LogDet = determinant(chunk, logarithm = TRUE)
what_j = fsvd$v[ ,q] * sqrt(n2)
# this is slow--fix
ssq = t(what) %*% solve(chunk) %*% what_j
dim(solve(chol(chunk)))
dim(what_j)
dim(chunk)
make_distmat = function(locs){
  d1 <- dist(locs)
  n <- dim(locs)[1]
  mat <- matrix(0,n,n)
  mat[lower.tri(mat)] <- d1
  mat <- mat+t(mat)
  return(mat)
}

get_log_det = function(part_one, covmat_j){
  
}

Estimate_params(list(1, c(2,2), c(1,2)))
params = list(1, c(2,2), c(1,2))
params[[1]]
params[[2]]

loglike = - .5*LogDet -n/2*log(2 *pi) - .5*ssqKronecker



Estimate_params = function(parameters, distC, distR, z){
  # parameters we want to estimate
  scaleC = parameters[1]
  scaleR = parameters[2]
  nugget = parameters[3]
  marginal_variance = parameters[4]
  
  #calculate eigendecomp of C and R 
  C = exp(-(distC / scaleC)^2) 
  R = exp(-(distR / scaleR)^2) 
  
  # do svd (eigen) decomp of both cov mats
  svdC = svd(C)
  Uc = svdC$u
  Sc = svdC$d
  Uc_t = t(svdC$v)
  
  R = exp(-(distR / scaleR)^2) 
  svdR = svd(R)
  Ur = svdR$u
  Sr = svdR$d
  Ur_t = t(svdR$v)
  
  MylogLik = LogLik_MVN(C, R, Uc, Sc, Uc_t, Ur, Sr, Ur_t, z, nugget, marginal_variance)
  # print(MylogLik)
  return(-MylogLik) #return negative log likelihood
}


