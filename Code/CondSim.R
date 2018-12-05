setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("kriging.R")
source("tbsim_2.R")
source("Generalised\ Increment.R")
source("GV.R")

#Returns coditional simulation of IRF with polynomial GC with coefficients b,
# conditional on points where bcicov[,q] is known. 
# n is number of lines in turning bands used,
# h is distance in grid of points for eval.

Condsim <- function(q,b,n,h){
  #Forming grid for sim
  x <- seq(0,1000, by = h)
  y <- seq(0, 500, by=h)
  grid.sim <- as.matrix(expand.grid(x,y)); k <- nrow(grid.sim)
  grid.list <- lapply(1:k, function(x){c(grid.sim[x,1], grid.sim[x,2])})
  krig.ind <- which(unlist(mclapply(grid.list, 
                                      function(x){isin2(matrix(x,1,2),grid[,1:2])}, mc.cores = 4)))
  #Kriged val of Z at x
  print("Kriging of Z: Started")
  Z.krig <- IntrinsicKriging(grid.sim[krig.ind,], grid[,3], grid.sim, b)
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

# Testing 

# q <- 3
# b <- c(2.89)
# n <- 100
# h <- 5 #Afstand af punkter i grid
# 
# sims <- Condsim(q,b,n,h)  


# Heatmap

# sim.matrix <- matrix(sims, ncol = 101, nrow= 201)
# 
# colours <- colorRampPalette(c("blue","cyan","yellow","red"))(100)
# 
# levelplot(sim.matrix,
#           row.values = seq(0,1000,by = 5),
#           column.values = seq(0,500,by = 5),
#           xlab = "x", ylab="y",
#           col.regions = colours)
# 
# krig <- IntrinsicKriging(bcicov[,1:2],data = bcicov$Al,inter.points = bcicov[,1:2],b)
# all.equal(krig,bcicov$Al)
