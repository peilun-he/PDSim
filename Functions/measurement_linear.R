measurement_linear <- function(x, par, mats) {
  # Linear measurement equation
  # Inputs: 
  #   x: x_t
  #   par: a vector of parameters and model coefficients
  #   mats: maturities
  # Outputs:
  #   y: y_t
  #   y_jacobian: Jacobian of y 
  
  n_obs <- dim(mats)[1]
  n_contract <- dim(mats)[2]
  
  y <- matrix(0, nrow = n_obs, ncol = n_contract)
  y_jacobian <- array(0, dim = c(2, n_contract, n_obs)) 
  
  for (i in 1: n_obs) {
    y_jacobian[, , i] <- rbind( exp(-par[1]*mats[i, ]), exp(-par[2]*mats[i, ]) )
    y[i, ] <- t(x[, i]) %*% y_jacobian[, , i] + AofT(par, mats[i,])
  }

  return(list(y = y, y_jacobian = y_jacobian))
}