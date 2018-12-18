setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("kriging.R")
source("tbsim_2.R")
source("Generalised\ Increment.R")
source("isin.R")

Condsim <- function(obs, obs.points, b, n, h){
  #Forming grid for sim
  x <- seq(0,1000, by = h)
  y <- seq(0, 500, by=h)
  grid.sim <- as.matrix(expand.grid(x,y)); k <- nrow(grid.sim)
  grid.list <- lapply(1:k, function(x){c(grid.sim[x,1], grid.sim[x,2])})
  #Which points of the simulation grid are observed.
  krig.ind <- isin(obs.points, grid.sim, bool = F)
  #Kriged val of Z at x
  print("Kriging of Z: Started")
  obs_ind <- isin(grid.sim[krig.ind,], obs.points, bool = F)
  Z.krig <- IntrinsicKriging(grid.sim[krig.ind,], obs[obs_ind], grid.sim, b)
  print("Kriging of Z: Done")
  #Sim S(x)
  s <- tbsim(b,n, grid.list)
  S.x <- cbind(grid.sim,s)
   #Kriged val of S at x
  print("Kriging of S: Started")
  S.krig <- IntrinsicKriging(S.x[krig.ind,1:2],S.x[krig.ind,3], grid.sim, b)
  print("Kriging of S: Done")
  return(Z.krig + (S.x[,3] - S.krig))
}
