
irf_sim_1d <- function(n,b,distance){
  # n is number of points, b is the vector of coefficients,
  # distance is the distance between the points.
  points <- seq(0,n*distance,distance)
  k <- length(b)-1
  X <- rnorm(n-1, sd=sqrt(distance)); W0 <- c(0,cumsum(X))
  if(k==0){
    Z <- 2*sqrt(b)*W0
    xw <- matrix(0,nrow = length(Z)+1,ncol = 2)
    xw[,1] <- points; xw[,2] <- c(0,Z)
    return( xw ) 
  }
  W <- matrix(0,n,k+1); W[,1] <- W0
  for(i in 1:k){
    W[,(i+1)] <- distance*cumsum(W[,i])
  }
  a <- rep(0,k+1)
  if(k==1){
    a[1] <- -sqrt(2*b[1])
    a[2] <- sqrt(2*b[2])
  }
  if(k==2){
    a[1] <- -sqrt(2*b[1])
    a[3] <- -sqrt(2*b[3])
    a[2] <- sqrt(2*(b[2] + 2*a[1]*a[3] ))
  }
  Z <- W %*% matrix(a)
  xw <- matrix(0,nrow = length(Z)+1,ncol = 2)
  xw[,1] <- points; xw[,2] <- c(0,Z)
  return(xw)
}



