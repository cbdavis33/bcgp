getCorMatR <- function(x, rho){
  
  # # If I choose not to export this function, then I'll skip the error-checking
  # # since the only time this function would be called is if everything is correct.
  # stopifnot(is.matrix(x), is.numeric(x), (0 <= rho && rho <= 1),
  #           ncol(x) == length(rho))
  
  n <- nrow(x)
  d <- ncol(x)
  
  R <- matrix(0, nrow = n, ncol = n)
  
  for(i in 1:(n-1)){
    R[i, i] <- 1.0
    for(j in (i+1):n){
      R[i, j] <- 1.0
      dist <- x[i, ] - x[j, ]
      R[i, j] <- prod(rho^(16*dist^2))
      R[j, i] <- R[i, j]
    }
    
  }
  R[n, n] <- 1.0
  return(R)
}


getCovMatR <- function(V, R, sig2){
  rootV <- sqrt(V)
  C <- t(rootV*R) * rootV + diag(sig2, length(V))
  return(C)
}