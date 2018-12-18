setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("Generalised\ increment.R")
source("GV.R")

gen.incX <- function(h,k,q){
  points <- as.matrix(bcicov[,1:2])
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
  out <- lapply(points.list, function(x){gen.inc(x,h,k,q)})
  return(list(Val = unlist(out),points = S.h))
}
