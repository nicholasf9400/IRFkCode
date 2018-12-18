# Estimation of Generalised Variogram (GV).

#Packages and Data import.
#-------------------------------
library(pracma)
library(hydroGOF)
library(parallel)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("Generalised\ increment.R")
source("isin.R")
#-------------------------------

M <- function(k){
  nchoosek(2*k+2,k+1)
}


gen.inc <- function(x,h,k,obs, obs.points){
  V <- function(x){
    temp <- which(obs.points[,1]==x[1] & obs.points[,2]==x[2])
    return(obs[temp])
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
  if(length(o)==0){return("Argument obs.point should be a regular grid")}
  return(o)
}
#-------------------------------

GenVar <- function(h,k=0,obs, obs.points){
  if (!any(mod(h,50)==0)) {
    stop("H must be multiple of 50")
  }
  points <- as.matrix(obs.points)
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
    logic[i] <- isin(transed.p, points) 
  }
  
  ind <- 1:nrow(points); ind <- ind[logic]
  S.h <- points[ind,]
  m <- length(ind)
  if(m==0){
    stop("Argument obs.points should be a regular grid")
  }
  if(m==1){
    points.list <- list(S.h)
  }else{
    points.list <- lapply(1:m, function(x){S.h[x,]})
  }
  
  sum.list <- lapply(points.list, function(x){gen.inc(x,h,k,obs,obs.points)^2})
  
  (out <- 1/(M(k)*m)*sum(unlist(sum.list)))
}

## Covariance mod
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

#Empirical GV for all h Works only for BCICOV dataset 
# 
IRFVariogram <- function(k,obs, obs.points,h){
  init_check_temp <- as.matrix(obs.points - apply(obs.points,2,min))
  if (any(mod(init_check_temp,h)!=0)) {
    warning("Either grid or grid distance (h) wrongly specified")
  }
  l <- floor(length(unique(obs.points[,1]))/(k+1))
  m <- floor(length(unique(obs.points[,2]))/(k+1))
  h.x <- seq(0, l-1, by=1)*h
  h.y <- seq(0, m-1, by=1)*h
  h.all <- as.matrix(expand.grid(h.x,h.y))
  h.l <- lapply(1:nrow(h.all), function(x){h.all[x,]})
  GV.emp <- unlist(mclapply(h.l, function(x){GenVar(x,k,obs, obs.points)},mc.cores = 4))
  out <- list(h = h.l, gamma = GV.emp, k=k)
  return(out)
}
