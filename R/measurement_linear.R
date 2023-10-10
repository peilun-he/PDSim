#' Linear measurement equation
#' 
#' Measurement equation in a linear form. This function can be passed to [simulate_data]
#' as the argument of `func_g`. 
#' 
#' @param x State variables \eqn{x_t}. 
#' @param par A vector of parameters and model coefficients. 
#' @param mats Time to maturity. 
#' 
#' @return This function returns a list with components: 
#' \item{y}{The futures price \eqn{y_t}. }
#' \item{y_jacobian}{The jacobian of \eqn{y_t}. }
#' 
#' @export
#' @seealso [measurement_polynomial], [state_linear] for other forms of state and 
#' measurement equations. 

measurement_linear <- function(x, par, mats) {
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