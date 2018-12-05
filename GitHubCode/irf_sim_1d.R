
irf_sim_1d <- function(n,b,distance){
  # n is number of points, b is the vector of coefficients,
  # distance is the distance between the points.
  points <- seq(0,n*distance,distance)
  k <- length(b)-1
  U <- runif(n,-1,1); R <- sign(U)
  W0 <- cumsum(R)
  if(k==0){
    Z <- -2*distance*sqrt(b)*W0
    xw <- matrix(0,nrow = length(Z)+1,ncol = 2)
    xw[,1] <- points; xw[,2] <- c(0,Z)
    return( xw ) 
    }
  W <- matrix(0,n,k+1); W[,1] <- W0
  for(i in 1:k){
    W[,(i+1)] <- distance*cumsum(W[,1])
  }
  a <- rep(0,k+1)
  if(k==1){
    a[1] <- -2*distance*sqrt(b[1])
    a[2] <- 2*distance*sqrt(b[2])
  }
  if(k==2){
    a[1] <- -2*distance*sqrt(b[1])
    a[3] <- 2*distance*sqrt(b[3])
    a[2] <- -2*distance*sqrt(b[2] + 2*b[1]*b[3] )
  }
  Z <- W %*% matrix(a)
  xw <- matrix(0,nrow = length(Z)+1,ncol = 2)
  xw[,1] <- points; xw[,2] <- c(0,Z)
  return(xw)
}

# Z <- irf_sim_1d(n=10000,b=c(568,0.0014),distance = 0.1)
# plot(Z,type = 'l')



