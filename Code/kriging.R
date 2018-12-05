setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# bcicov <- read.csv("Data/bci.dataJul05.csv")
# source("StructureAnal.R")
Rcpp::sourceCpp('GenCovMatrix.cpp')
source("Generalised\ increment.R")
library(limSolve)
library(parallel)
# library(hydroGOF)

# For debugging
# points <- bcicov[1:3,1:2]
# inter.points <- bcicov[4,1:2]
# data <- bcicov$Zn[1:3]
# IntrinsicKriging(points,data,inter.points,b)

#Returns values intrinsic kriging estimator at "inter.points", 
# computed using observations "data",
# at points "points", using polynomial GC with coefficients "b".
#### "inter.points" may be matrix
IntrinsicKriging <- function(points,data,inter.points,b){
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
  
  C <- min(detectCores(), 20); C <- min(C, q)
  inter.points.list <- split(inter.points, rep(1:nrow(inter.points), each = ncol(inter.points)))
  p <- mclapply(inter.points.list, function(x) as.numeric(polynomial_matrix(x,k)),mc.cores = C )
  
  Big.vec <- lapply(as.list(1:q), function(i) c(h[[i]],p[[i]]))
  lambdamu <- lapply(Big.vec, function(x) Solve(Big.matrix,x))
  lambda <- mclapply(lambdamu, function(x) x[1:n], mc.cores = C)
  krig <- lapply(lambda, function(x) as.numeric(x %*% data))
  z <- unlist(krig)
  return(z)
}

# Cross-validation
# cv.list <- list()
# cv.data <- bcicov[,c(1,2,14)]
# cv.index <- 1:300
# rear <- NULL
# for(i in 1:10){
#   id <- sample(cv.index,30)
#   cv.list[[i]] <- id
#   cv.index <- cv.index[!(cv.index %in% id)]
#   rear <- c(rear,id)
# }
# cv.cal <- lapply(cv.list, function(x) IntrinsicKriging(cv.data[-x,1:2],cv.data[-x,3],cv.data[x,1:2],b))
# comp <- as.numeric(cv.data[rear,3])
# rmse(unlist(cv.cal),comp)
# mae(unlist(cv.cal),comp)
