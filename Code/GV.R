# Estimation of Generalised Variogram (GV).

#Packages and Data import.
#-------------------------------
library(pracma)
library(hydroGOF)
library(parallel)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("Generalised\ increment.R")
#-------------------------------

M <- function(k){
  nchoosek(2*k+2,k+1)
}

#Define grid
v <- (bcicov$gx-25)/50
keep.x <- which(v-floor(v)==0)
grid <- bcicov[keep.x,]

w <- (grid$gy-25)/50
keep.y <- which(w-floor(w)==0)
grid <- grid[keep.y,]
rm(v,keep.x,w,keep.y)

# Function determines if the matrix x is made up entirely
# from rows in M. requires ncol(x)==nrow(M).
isin2 <- function(x,M){
  l <- nrow(x); n <- nrow(M)
  logi <- logical(length = l)
  for (j in 1:l) {
    for (i in 1:n) {
      if(all(x[j,]==M[i,])){
        logi[j] <- T
      }
    }
  }
  if(all(logi)){
    return(T)
  }
  return(F)
}

#Returns the k+1st diff operator of bcicov[,q], at x with translations h.
gen.inc <- function(x,h,k,q){
  V <- function(x){
    x1 <- x[1]
    x2 <- x[2]
    Z(x1,x2,n=q)
  }
  if(any(mod(h,50)!=0)){
    stop("h must be multiple of 50")
  }
  if(k==0){
    o <- V(x+h) - V(x)
  }
  if(k==1){
    o <- V(x+2*h) - 2*V(x+h) + V(x)
  }
  if(k==2){
    o <- V(x+3*h) - 3*V(x+2*h) + 3*V(x+h) - V(x)
  }
  if(length(o)==0){return("hej")}
  return(o)
}
#-------------------------------

# h <- c(50,50)
# k <- 2
# q <- 3

#The kth order empirical generalised variogram at distance h, of bcicov[,q].
GenVar <- function(h,k=0,q=3){
  if (!any(mod(h,50)==0)) {
    stop("H must be multiple of 50")
  }
  points <- as.matrix(grid[,1:2])
  l <- nrow(points)
  trans <- matrix(NA,2,k+2); trans[,1] <- c(0,0)
  for (i in 1:(k+1)){
    trans[,i+1] <- i*h
  }
  
  S.h <- matrix(NA,2, l)
  logic <- logical(length = l)
  for (i in 1:l){
    p <- as.numeric(points[i,])
    transed.p <- t(as.matrix(p + trans))
    logic[i] <- isin2(transed.p, points) 
  }
  
  ind <- 1:nrow(points); ind <- ind[logic]
  S.h <- points[ind,]
  m <- length(ind)
  if(m==0){
    return(NA)
  }
  if(m==1){
    points.list <- list(S.h)
  }else{
    points.list <- lapply(1:m, function(x){S.h[x,]})
  }
  
  sum.list <- lapply(points.list, function(x){gen.inc(x,h,k,q)^2})
  
  (out <- 1/(M(k)*m)*sum(unlist(sum.list)))
}

# l <- seq(0,6, by = 1)*50
# eval.list <- lapply(1:7, function(x){c(l[x],0)})
# 
# y <- unlist(lapply(eval.list, function(x){GenVar(x,k=0,q=10)}))
# 
# plot(l,y, type = 'l')

## Covariance models
KMod <- function(h,b,k){
  if (k==0) {
    return(b*(-1)*Norm(h))
  }
  if(k==1){
    return(b[1]*(-1)*Norm(h)+b[2]*(-1)^2*Norm(h)^3)
  }
  if(k==2){
    return(b[1]*(-1)*Norm(h)+b[2]*(-1)^2*Norm(h)^3 + b[3]*(-1)*Norm(h)^5)
  }
}
KMod1 <- function(h,b,k){
  red <- k
  return(b[1]*exp(-Norm(h)/b[2]))
}

#Relation between GV and GC.
model <- function(k,h,b){
  if(k==0){
    o <- 1/2*(-KMod(-h,b,k)+2*KMod(0,b,k)-KMod(2*h,b,k))
  }
  if(k==1){
    o <- 1/6*(2*KMod(2*h,b,k) - 8*KMod(h,b,k) + 6*KMod(0,b,k))
  }
  if(k==2){
    o <- 1/20*(20*KMod(0,b,k)-30*KMod(h,b,k)+12*KMod(2*h,b,k)-2*KMod(3*h,b,k))
  }
  return(o)
}

#Empirical GV for all h
IRFVariogram <- function(k,q){
  l <- floor(length(unique(grid$gx))/(k+1))
  m <- floor(length(unique(grid$gy))/(k+1))
  h.x <- seq(0, l-1, by=1)*50
  h.y <- seq(0, m-1, by=1)*50
  h.all <- as.matrix(expand.grid(h.x,h.y))
  h.l <- lapply(1:nrow(h.all), function(x){h.all[x,]})
  GV.emp <- unlist(mclapply(h.l, function(x){GenVar(x,k,q)},mc.cores = 4))
  out <- list(h = h.l, gamma = GV.emp, k=k)
  return(out)
}

# #Fitting GC to EGV
# 
# k=0 #Order
# q=3 #Column in BCICOV
# 
# emp <- IRFVariogram(k,q)
# 
# optim.func <- function(parm){
#   temp <- unlist(lapply(emp$h, function(x){model(emp$k, x, parm)}))
#   return(mse(emp$gamma, temp))
# }
# 
# paramter <- optim(par = rep(0,k+1), fn = optim.func)$par
# 
# 
# #dadsadsa
# l <- floor(length(unique(grid$gx))/(k+1))
# h.x <- seq(0, l-1, by=1)*50
# h <- lapply(1:l, function(x){c(h.x[x],0)})
# EGV <- unlist(lapply(h, function(x){GenVar(x,k,q)}))
# plot(h.x, EGV, type = 'l', main = paste("Ressource", q, "Order", k))
# modelGV <- unlist(lapply(h, function(x){model(k,x,paramter)}))
# lines(h.x,modelGV, col = "red")

