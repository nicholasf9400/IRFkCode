### FINAL STRUCTURE ANALYSIS SCRIPT

#Dependencies and data
library(pracma)
library(rstudioapi) # Kræver pakken rstudioapi
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("PointsFun.R")
source("gen.incX.R")
#------------------


##Determining k
IRFOrder <- function(obs, obs.points){
  points <- obs.points
  l <- nrow(points)
  errmat <- matrix(NA, nrow = 3*l, ncol = 2); rownames(errmat) <- 1:(3*l)
  o <- 0
  for (k in 0:2) {
    m <- factorial(k+2)/(factorial(2)*factorial(k))+1
    for (i in 1:(l-m)){
      o <- o+1
      lincomp.points <- points[i+(1:m),]
      lincomp.val <- obs[i+(1:m)]
      val <- obs[i]
      #K <- polynomial_matrix(lincomp.points,k)
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
        fun <- function(par){
          lma <- lambda%*%par
          as.numeric((lincomp.val%*%lma-val)^2)
        }
        if(n==1){
          parameter <- optimize(fun, interval = c(-500.1,500))$minimum
        }else{
          parameter <- optim(par = rep(0,n), fn = fun)$par
        }
        lambda.opt <- lambda%*%parameter
        errmat[o,] <- c(lincomp.val%*%lambda.opt-val,k)
      }
    }
  }
  perm <- order(errmat[,1])
  errorsort <- errmat[perm,]
  out <- c(mean(which(errorsort[,2]==0)),
           mean(which(errorsort[,2]==1)),
           mean(which(errorsort[,2]==2)))
  return(which.min(out)-1)
}

EstimateGC <- function(obs,k,n, Points){
  
  #Choosing points used for computing ALC
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
    d <- t(as.matrix(obs[ind]))
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
        for (w in 1:l) {
          k.p <- c(k.p, ALCs.allK[[j]][m]*ALCs.allK[[j]][w]*K(p,as.numeric(ALCPoints.allK[[j]][m,] - ALCPoints.allK[[j]][w,])))
        }
      }
      covar[j,p+1] <- sum(k.p)
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
      if(any(is.na(fit$coefficients))){}else{
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
        
      }}
    }
    
  }
  
  
  Q <- function(b){
    d <- as.matrix(reg.data[,-1])
    sum(1/var(res)^2*(res-d%*%b)^2)
  }
  
  ind <- sum(!is.na(coefmat[,1]))
  critmat <- rep(0,ind)
  for (i in 1:ind) {
    critmat[i] <- Q(coefmat[i,])/Q(rep(0,k+1))
  }
  
  best <- which.min(critmat)
  
  return(list(Coefficients = coefmat[best,], FitCriteria = critmat[best]))
  
}

######Estimation of GC ved brug af GV

GVFitGC <- function(k,obs, obs.points, h){
  emp <- IRFVariogram(k,obs, obs.points, h)
  
  optim.func <- function(parm){
    temp <- unlist(lapply(emp$h, function(x){model(emp$k, x, parm)}))
    return(mse(emp$gamma, temp))
  }
  
  parameter <- optim(par = rep(0,k+1), fn = optim.func, lower = rep(0,k+1),
                     method = "L-BFGS-B")$par
  return(parameter)
}
