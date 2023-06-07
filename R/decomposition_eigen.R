decomposition_eigen <- function(A) {
  
  eigA <- eigen(A)
  D <- eigA$values
  V <- eigA$vectors
  expA <- V %*% diag(exp(D)) %*% solve(V)
  
  return(expA)
}

