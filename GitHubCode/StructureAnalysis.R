### FINAL STRUCTURE ANALYSIS SCRIPT

#Dependencies and data
library(pracma)
library(rstudioapi) # Kræver pakken rstudioapi
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("PointsFun.R")
source("gen.incX.R")
#------------------

## Note that the following only works for bcicov-dataset.

##Determining k
IRFOrder <- function(q){
  points <- bcicov[,1:2]
  l <- nrow(points)
  errmat <- matrix(NA, nrow = 3*l, ncol = 2); rownames(errmat) <- 1:(3*l)
  o <- 0
  for (k in 0:2) {
    m <- factorial(k+2)/(factorial(2)*factorial(k))+1
    for (i in 1:(l-m)){
      o <- o+1
      lincomp.points <- points[i+(1:m),]
      lincomp.val <- bcicov[i+(1:m),q]
      val <- bcicov[i,q]
      if(k==0){
        lambda <- c(1,-1)
        n <- 1
        fun <- function(par){
          lma <- par*lambda
          as.numeric((val - lincomp.val%*%lma)^2)
        }
        par1 <- optimize(f = fun, interval = c(-500.1,500))$minimum
        weight <- par1*lambda
        errmat[o,] <- c(lincomp.val%*%weight - val,k)
      }else{
        lambda <- ALC(lincomp.points,k, All = T)
        n <- ncol(lambda)
        lsq.par <- rep(0,n)
        lsq.val <- rep(0,n)
        for (c in 1:n) {
          fun <- function(par){
            lma <- par*lambda[,c]
            as.numeric((lincomp.val%*%lma-val)^2)
          }
          temp <- optimize(f = fun, interval = c(-500.1,500))
          lsq.val[c] <- temp$objective
          lsq.par[c] <- temp$minimum
        }
        lsq.lambda <- which.min(lsq.val)
        par <- lsq.par[lsq.lambda]
        weight <- par*lambda[,lsq.lambda]
        errmat[o,] <- c(lincomp.val%*%weight-val,k)
      }
    }
  }
  perm <- order(errmat[,1])
  errorsort <- errmat[perm,]
  out <- c(mean(which(errorsort[,2]==0)),
           mean(which(errorsort[,2]==1)),
           mean(which(errorsort[,2]==2)))
  return(out)
}

IRFOrder(3)


#####Estimate GC##### Using all data points, even non-grid##
#Ends line 163

k <- 1
q <- 3
n <- 100

ALCPoints.allK <- Points(k,n)

## Computing corresponding ALC
ALCs.allK <- mclapply(ALCPoints.allK, 
                               function(x){ALC(x,k,F)}, mc.cores = 4)
## Elementary model definition
K <- function(p,h){
  (-1)^(p+1)* Norm(h)^(2*p+1)
}

## Data transformation
G <- function(i){
  ind <- as.numeric(rownames(ALCPoints.allK[[i]]))
  d <- t(as.matrix(bcicov[ind, q]))
  res <- (d%*%ALCs.allK[[i]])^2
  return(res)
}


v <- rep(0,n)
for (i in 1:n) {
  v[i] <- G(i)
}
res <- v

#Generating Covariates

covar <- as.data.frame(matrix(0, n, k+1))
for (j in 1:n) {
  for (p in 0:(k)) {
      k.p <- c()
      l <- nrow(ALCPoints.allK[[j]])
      for (m in 1:l) {
        k.p <- c(k.p, K(p,as.numeric(ALCPoints.allK[[j]][m,])))
      }
      covar[j,p+1] <- k.p%*%ALCs.allK[[j]]
  }
}
reg.data <- as.data.frame(cbind(res, covar))

data.row <- nrow(reg.data)

coefmat <- matrix(NA, ncol = k+1, nrow=2000)
v <- 1
for (i in 3:data.row) {
  o <- floor(data.row/i)
  idxs <- 1:data.row
  for (j in 1:o){
    sample.ind <- sample(idxs, i)
    idxs <- idxs[-sample.ind]
    sam.data <- reg.data[sample.ind,]
    fit <- lm(res~.-1, data = sam.data, weights = rep(var(res),i))
    if((k<=1)&(all(fit$coefficients>=0))){
      coefmat[v,] <- c(fit$coefficients)
      v <- v+1
    }
    if((k==2)){
      b <- fit$coefficients
      b02.cond<- (b[1]>=0) & (b[3]>=0)
      if(b02.cond){
        b1.cond <- b[2] >= -(10/3)*sqrt(b[1]*b[3])
        if(b1.cond){
          coefmat[v,] <- b
          v <- 1+v
        }
      }
      
    }
  }
    
}

#Model fit criteria
Q <- function(b){
  d <- as.matrix(reg.data[,-1])
  sum(1/var(res)^2*(res-d%*%b)^2)
}

#Calculating criteria for all estimated coefficients.
ind <- sum(!is.na(coefmat[,1]))
critmat <- rep(0,ind)
for (i in 1:ind) {
  critmat[i] <- Q(coefmat[i,])/Q(rep(0,k+1))
}

best <- which.min(critmat)

coefmat[best,]; critmat[best]
  p1 <- coefmat[best,]

#######Estimation af GC####### using only grid points ####
# ends line 254

k=1
q=3

K <- function(p,h){
  (-1)^(p+1)* Norm(h)^(2*p+1)
}
h <- c(50,0)
ALC55 <- gen.incX(h,k,q)


respons <- c(ALC55$Val)
#points <- rbind(ALC55$points, ALC50$points, ALC05$points)

l <- length(ALC55$Val)
points <- ALC55$points

covar <- rep(0,l)

difalc <- function(k){
  vec <- rep(0,k+2)
  for (i in 0:(k+1)) {
    vec[i+1] <- (-1)^(i)*nchoosek(k+1,i)
  }
  return((-1)^(k+1)*vec)
}


mat <- matrix(0, nrow = l, ncol = k+1)
for (p in 0:k) {
  for (i in 1:l) {
    kCov <- rep(0,p+1)
    for (j in 0:(p+1)) {
      kCov[j+1] <- K(p,points[i,]+j*h)
    }
    mat[i,p+1] <- as.numeric(difalc(p)%*%kCov)
  }
}

reg.data <- as.data.frame(cbind(respons, mat))

data.row <- nrow(reg.data)

coefmat <- matrix(NA, ncol = k+1, nrow=2000)
v <- 1
for (i in 2:data.row) {
  o <- floor(data.row/i)
  idxs <- 1:data.row
  for (j in 1:o){
    sample.ind <- sample(idxs, i)
    idxs <- idxs[-sample.ind]
    sam.data <- reg.data[sample.ind,]
    fit <- lm(respons~.-1, data = sam.data, weights = rep(var(respons),i))
    if((k<=1)&(all(fit$coefficients>=0))){
      coefmat[v,] <- c(fit$coefficients)
      v <- v+1
    }
    if((k==2)){
      b <- fit$coefficients
      b02.cond<- (b[1]>=0) & (b[3]>=0)
      if(b02.cond){
        b1.cond <- b[2] >= -(10/3)*sqrt(b[1]*b[3])
        if(b1.cond){
          coefmat[v,] <- b
          v <- 1+v
        }
      }
      
    }
    
  }
}
Q <- function(b){
  d <- as.matrix(reg.data[,-1])
  sum(1/var(reg.data[,1])^2*(reg.data[,1]-d%*%b)^2)
}


ind <- sum(!is.na(coefmat[,1]))
critmat <- rep(0,ind)
for (i in 1:ind) {
  critmat[i] <- Q(coefmat[i,])/Q(rep(0,k+1))
}

best <- which.min(critmat)

coefmat[best,]; critmat[best]
p2 <- coefmat[best,]


######Estimation of GC ved brug af GV
#ends line 271

k=1 #Order
q=3 #Column in BCICOV

emp <- IRFVariogram(k,q)

optim.func <- function(parm){
  temp <- unlist(lapply(emp$h, function(x){model(emp$k, x, parm)}))
  return(mse(emp$gamma, temp))
}

parameter <- optim(par = rep(0,k+1), fn = optim.func)$par
p3 <- parameter


#Plotting (may not work)
l <- floor(length(unique(grid$gx))/(k+1))
h.x <- seq(0, l-1, by=1)*50
h <- lapply(1:l, function(x){c(h.x[x],0)})
EGV <- unlist(lapply(h, function(x){GenVar(x,k,q)}))
plot(h.x, EGV, type = 'l', main = paste("Aluminium", "k = ", k),ylim = c(0,180000))
model.p3 <- unlist(lapply(h, function(x){model(k,x,p3)}))
lines(h.x,model.p3, col = "red")
model.p2 <- unlist(lapply(h, function(x){model(k,x,c(1.441206e+00, 2.429460e-07))}))
lines(h.x,model.p2, col = "blue")
model.p1 <- unlist(lapply(h, function(x){model(k,x,c(6.011120e+05, 6.391678e-01))}))
lines(h.x,model.p1, col = "blue")

