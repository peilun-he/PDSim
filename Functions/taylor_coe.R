taylor_coe <- function(order) {
  # Coefficients of Taylor series of f(chi, xi) = exp(chi+xi) at (0, 0)
  # Inputs: 
  #   Order: the order of Taylor series
  # Outputs:
  #   coe: a vector of coefficients
  
  coe <- 1
  
  if (order == 0) {
    return(0)
  } else {
    for (s in 1: order) {
      for (i in seq(from = s, to = 0, by = -1)) {
        j <- s - i
        coe <- c( coe, 1 / (factorial(i) * factorial(j)) )
      }
    }
  }
  
  return(coe)
}