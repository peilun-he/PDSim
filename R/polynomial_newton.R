polynomial_newton <- function(A) {
  
  eigA <- eigen(A)
  lambda <- eigA$values
  n <- dim(A)[1]
  expA <- exp(lambda[1]) * diag(n)
  
  for (j in 2: n) {
    prod <- diag(n)
    for (k in 1: (j-1)) {
      prod <- prod %*% (A - lambda[k] * diag(n))
    }
    expA <- expA + divided_difference(lambda[1: j]) * prod
  }
  
  return(expA)
}

divided_difference <- function(lambda) {
  n <- length(lambda)
  if (n == 2) {
    l1 <- lambda[1]
    l2 <- lambda[2]
    dd <- ( exp(l1) - exp(l2) ) / (l1 - l2)
  } else {
    dd <- ( divided_difference(lambda[1: (n-1)]) - divided_difference(lambda[2: n]) ) / ( lambda[1] - lambda[n] )
  }
  return(dd)
}



