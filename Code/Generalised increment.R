setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(MASS)

#Data Import
bcicov <- read.csv("bci.dataJul05.csv")

# Generalized increments

#Function returns values of the nth column of bcicov at 
# input points x and y.

Z <- function(x,y,n=3){
  temp <- which(bcicov$gx==x & bcicov$gy==y)
  return(bcicov[temp,n])
}

# Function constructs a matrix of all monomial of deg =< k,
# evaluated at points
polynomial_matrix <- function(points,k){
  # Converts a numeric vector into an appropriate matrix
  if(class(points)=="numeric"){
    points <- t(as.matrix(points))
  }
  # Computes the number of two dimensional polynomials up to degree k
  # and sets l to the number of points
  kn <- factorial(k+2)/(factorial(2)*factorial(k))
  l <- nrow(points)
  
  # Constructing all multi-indices with sum 0,1,...,k
  t <- 1
  L <- list()
  for(i in 0:k){
    for(j in 0:i){
      L[[t]] <- c(j,i-j)
      t <- t+1
    }
  }
  
  X <- matrix(1,l,kn)
  for(i in 1:l){
    x <- as.numeric(points[i,])
    p <- lapply(L, function(n) prod(x^n))
    X[i,] <- unlist(p)
  }
  return(X)
}


# Computes allowable measure using "points".
# Returns 1 allowable measure if "ALL=F" and all allowable
# measures if "All=T".
ALC <- function(points,k,All=TRUE){
  # Computes the number of two dimensional polynomials up to degree k
  # and sets l to the number of points
  kn <- factorial(k+2)/(factorial(2)*factorial(k))
  l <- nrow(points)
  
  # Sets up the polynomial matrix
  X <- polynomial_matrix(points,k)
  
  r <- qr(X)$rank
  if(r==l){
    stop("Too few points chosen")
  }
  
  All_lambda <- Null(X)
  lambda <- All_lambda[,1]
  
  if(All){
    return(All_lambda) 
  }
  if(!All){
    return(lambda)
  }
}


#Uses allowable measures from "ALC"-function, to compute ALCs
 # from the bcicov dataset.
bcicov_grid_ALC <- function(lambda, points){
  l <- nrow(points)
  lambda <- as.numeric(lambda)
  # Remove all the points that are not on the grid
  v <- (bcicov$gx-25)/50
  keep.x <- which(v-floor(v)==0)
  grid <- bcicov[keep.x,]
  
  w <- (grid$gy-25)/50
  keep.y <- which(w-floor(w)==0)
  grid <- grid[keep.y,]
  
  # Computes the number of steps that can be taken in the x- an 
  # y-direction without stepping outside the grid
  xstep.neg <- (25-min(points[,1]))/50
  xstep.pos <- (975-max(points[,1]))/50
  ystep.neg <- (25-min(points[,2]))/50
  ystep.pos <- (475-max(points[,2]))/50
  n.x <- xstep.pos+abs(xstep.neg)+1
  n.y <- ystep.pos+abs(ystep.neg)+1
  
  # constructs the dataset and put the valid coordinates in the first two columns
  gen.inc <- data.frame(matrix(0,n.x*n.y,21))
  colnames(gen.inc) <- colnames(bcicov)
  gen.inc$gx <- 50*seq(xstep.neg,xstep.pos,1)
  gen.inc$gy <- 50*seq(ystep.neg,ystep.pos,1); gen.inc$gy <- sort(gen.inc$gy)
  
  for(n in 3:21){
    for(i in 1:(n.x*n.y)){
      geninc <- rep(0,l)
      for(j in 1:l){
        x <- points[j,1]+gen.inc$gx[i]
        y <- points[j,2]+gen.inc$gy[i]
        geninc[j] <- Z(x,y,n)
      }
      gen.inc[i,n] <- lambda %*% geninc
    }
  }
  return(list(Inc = gen.inc, Coef = lambda))
}

# For testing purposes
#L2 <- matrix(0,7,2); L2[,1] <- c(25,75,125,25,75,125,75)+50; L2[,2] <- c(25,25,25,75,75,75,125)+50
# ALC(L,2)

