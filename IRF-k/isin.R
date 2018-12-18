# Is in 

isin <- function(X,M, bool = T){
  l <- nrow(X); n <- nrow(M)
  logi <- logical(length = l)
  ind <- c()
  for (j in 1:l) {
    for (i in 1:n) {
      if(all(X[j,]==M[i,])){
        logi[j] <- T
        ind <- c(ind,i)
        break
      }
    }
    
  }
  if(bool){
    if(all(logi)){
      return(T)
    }
    return(F)
  }else{
    return(ind)
  }
}
