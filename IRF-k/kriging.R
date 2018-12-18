setwd(dirname(rstudioapi::getSourceEditorContext()$path))
Rcpp::sourceCpp('GenCovMatrix.cpp')
source("Generalised\ increment.R")
library(limSolve)
library(parallel)

IntrinsicKriging <- function(points,data,inter.points,b,...){
  points <- as.matrix(points); inter.points <- as.matrix(inter.points); data <- as.numeric(data)
  k <- length(b)-1
  kn <- factorial(k+2)/(factorial(2)*factorial(k))
  n <- nrow(points);  q <- nrow(inter.points)
  
  H <- GenCovMatrix(as.matrix(points),b)
  K0 <- genCov(c(0,0),b)
  diag(H) <- K0
  H <- H+t(H)
  P <- polynomial_matrix(points,k)
  
  Big.matrix <- matrix(0,n+kn,n+kn)
  Big.matrix[1:n,1:n] <- H
  Big.matrix[1:n,(n+1):(n+kn)] <- P
  Big.matrix[(n+1):(n+kn),1:n] <- t(P)
  
  h <- lapply(as.list(1:q), function(i) GenCovVector(points, inter.points[i,], b))
  
  inter.points.list <- split(inter.points, rep(1:nrow(inter.points), each = ncol(inter.points)))
  p <- mclapply(inter.points.list, function(x) as.numeric(polynomial_matrix(x,k)),...)
  
  Big.vec <- lapply(as.list(1:q), function(i) c(h[[i]],p[[i]]))
  lambdamu <- mclapply(Big.vec, function(x) Solve(Big.matrix,x),...)
  lambda <- lapply(lambdamu, function(x) x[1:n])
  krig <- lapply(lambda, function(x) as.numeric(x %*% data))
  z <- unlist(krig)
  return(z)
}
