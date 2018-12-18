#Points Function

library(pracma)
library(parallel)
library(rstudioapi) # Kræver pakken rstudioapi
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("Generalised\ increment.R")

Points <- function(k,n){
  locations=bcicov[,1:2]
  l <- nrow(locations)
  kn <- factorial(k+2)/(factorial(2)*factorial(k))
  points.list <- list()
  border_ind <- which(bcicov$gx<75 | bcicov$gx>925 | bcicov$gy<75 | bcicov$gy>425)
  no_border <- locations[-border_ind,]
  for(i in 1:n){
    m <- 0
    while(m<kn){ # Usually only one iteration is neccesary
      s <- sample((1:300)[-border_ind],1);
      x <- locations[s,1]; y <- locations[s,2]
      ind <- which(bcicov$gx >= x-50 & bcicov$gx <= x+50
                   & bcicov$gy >= y-50 & bcicov$gy <= y+50
                   & bcicov$gx != x & bcicov$gx != x)
      if(length(ind) >= kn){
        pts_ind <- c(s,sample(ind,kn))
        points.list[[i]] <- locations[pts_ind,1:2]
        m <- Inf # To break the loop
      }
    }
  }
  return(points.list)
}