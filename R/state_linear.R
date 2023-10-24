#' Linear state equation
#' 
#' State equation in a linear form \eqn{x_t = A + B x_{t-1}}. This function can be passed to [simulate_data], 
#' [EKF] and [UKF] as the argument of `func_f`. 
#' 
#' @param x State variables \eqn{x_{t-1}} at previous time point \eqn{t-1}. 
#' @param par A vector of parameters and model coefficients. 
#' @param dt \eqn{\Delta t}. The interval between two consecutive time points. 
#' 
#' @return This function returns a list with components: 
#' \item{y}{The state variables \eqn{x_t} at current time point \eqn{t}. }
#' \item{y_jacobian}{The jacobian of \eqn{x_t}. }
#' 
#' @export
#' @seealso [measurement_linear], [measurement_polynomial] for other forms of state and 
#' measurement equations. 

state_linear <- function(x, par, dt) {
  
  kappa_chi  <- par[1]
  kappa_xi   <- par[2]
  mu_xi      <- par[3]
  sigma_chi  <- par[4]
  sigma_xi   <- par[5]
  rho        <- par[6]
  lambda_chi <- par[7]
  lambda_xi  <- par[8]
  
  A <- c( 0, mu_xi/kappa_xi*(1-exp(-kappa_xi*dt)) )
  B <- matrix(c( exp(-kappa_chi*dt), 0, 0, exp(-kappa_xi*dt) ), nrow = 2, byrow = TRUE)

  y <- A + B %*% x
  y_jacobian <- B
  
  return(list(y = y, y_jacobian = y_jacobian))
}

