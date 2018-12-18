##Turning Bands simulation for an IRF-k in R^2 with desired
# polynomial GCF with parameter vector b. b should be of length k+1.
## Furthermore a list of grid points in R^2 should be constructed. These should
# be equidistant.

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("irf_sim_1d.R")
source("decimalplaces.R")

B <- function(n=2,p){
  factorial(p)*gamma(n/2)*(sqrt(pi)*gamma(p+(1+n)/2))
}

tbsimx <- function(theta,dist,d, grid.list){
  y.eval <- unlist(lapply(grid.list, 
                          function(x){abs(as.numeric(x%*%c(cos(theta),sin(theta))))}))
  n.dig <- decimalplaces(dist)
  line.point <- round(y.eval, digits = n.dig)
  endpoint <- ceiling(max(line.point)/dist)
  line.sim <- irf_sim_1d(n=endpoint,d,dist)
  sim.ind <- line.point/dist+1
  out <- line.sim[sim.ind,2]
  return(out)
}

tbsim <- function(b,n.line, grid.list){
  theta <- runif(n.line, 0, pi)
  m <- length(b)
  d <- rep(0,m)
  for(i in 0:(m-1)){
    d[i+1] <- b[i+1]*factorial(2*i+1)/B(2,i)
  }
  sims <- lapply(theta, function(i){tbsimx(i, 0.01,d, grid.list)})
  out <- 1/sqrt(n.line)*Reduce("+",sims)
  grid.sim <- matrix(unlist(grid.list),ncol=2, byrow=T)
  return(cbind(grid.sim,out))
}
