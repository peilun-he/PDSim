#' Simulate commodity futures data
#' 
#' Simulate commodity futures price, time to maturity, and hidden factors based on state-space 
#' model. See `Details` for more information about the model. 
#' 
#' The state-space model is give by
#' \deqn{y_t = g(x_t, m_t) + u_t}
#' \deqn{x_t = f(x_{t-1}, m_t) + v_t}
#' where \eqn{y_t} is the futures price, \eqn{x_t} is the hidden state variable, \eqn{m_t} is 
#' the time to maturity, \eqn{u_t} and \eqn{v_t} are noises with mean 0. 
#' 
#' @param par A vector of parameters. 
#' @param x0 Initial values of state variables. 
#' @param n_obs The number of observations. 
#' @param n_contract The number of contracts. 
#' @param func_f Function `f`, which should take two arguments, xt and a vector of parameters, 
#' and return two values, the function value f(x) and the gredient f'(x). ( f'(x) is useless in 
#' this function, just make it consistent as used in other functions )
#' @param func_g Function `g`, which should take three arguments, xt, a vector of parameters 
#' and maturities, and return two values, g(x) and g'(x). ( g'(x) is useless in this function, 
#' just make it consistent as used in other functions )
#' @param n_coe The number of model coefficients.  
#' @param noise The distribution of noise, currently only "Gaussian" works. 
#' @param seed Integer. Seed for random values. 
#' 
#' @return This function returns a \code{list} with components: 
#' \item{yt}{A data frame. Commodity futures price.}
#' \item{mats}{A data frame. Time to maturity. It has the same dimension as `yt`.} 
#' \item{xt}{A data frame. Hidden state variables. }
#' 
#' @import MASS
#' @export
#' @seealso [state_linear], [measurement_linear], [measurement_polynomial] for examples of function 
#' `f` and `g`. 
#' @examples
#' n_obs <- 100
#' n_contract <- 10
#' par <- c(0.5, 0.3, 1, 1.5, 1.3, -0.3, 0.5, 0.3, seq(from = 0.1, to = 0.01, length.out = n_contract))
#' x0 <- c(0, 1/0.3)
#' dt <- 1/360 # daily data
#' n_coe <- 0
#' func_f <- function(xt, par) state_linear(xt, par, dt) # state equation
#' dat <- simulate_data(par, x0, n_obs, n_contract, func_f, measurement_linear, n_coe, "Gaussian", 1234)

simulate_data <- function(par, x0, n_obs, n_contract, func_f, func_g, n_coe, noise, seed) {

  # Parameters
  kappa_chi  <- par[1]
  kappa_xi   <- par[2]
  mu_xi      <- par[3]
  sigma_chi  <- par[4]
  sigma_xi   <- par[5]
  rho        <- par[6]
  lambda_chi <- par[7]
  lambda_xi  <- par[8]
  
  monthdays <- 30 # number of days per month
  yeardays <- 360 # number of days per year
  dt <- 1 / yeardays # delta_t
  n_state <- length(x0) # number of state variables
  n_para <- length(par) # number of parameters
  
  # Generate random noises for xt and yt
  set.seed(seed) # fix random seed
  
  if (noise == "Gaussian") {
    if (length(par) - 8 - n_coe != n_contract) {
      stop("Incorrect number of contracts or parameters. ")
    }
    if (n_contract == 1) {
      V <- par[9]^2
    } else {
      V <- diag( par[9: (n_para-n_coe)] ^ 2 )
    }
    W = matrix(c(sigma_chi^2/(2*kappa_chi) * ( 1-exp(-2*kappa_chi*dt) ), 
                 rho*sigma_chi*sigma_xi/(kappa_chi+kappa_xi) * ( 1-exp(-(kappa_chi+kappa_xi)*dt) ),  
                 rho*sigma_chi*sigma_xi/(kappa_chi+kappa_xi) * ( 1-exp(-(kappa_chi+kappa_xi)*dt) ), 
                 sigma_xi^2/(2*kappa_xi) * ( 1-exp(-2*kappa_xi*dt) )), 
               nrow = 2, byrow = TRUE)
    noise_xt = mvrnorm(n_obs, c(0, 0), W)
    noise_yt = mvrnorm(n_obs, rep(0, n_contract), V)
  } else {
    stop("Incorrect distribution of noises. ")
  }
  
  # Simulate xt
  xt <- matrix(0, nrow = n_obs+1, ncol = n_state)
  xt[1, ] <- x0
  for (j in 2: (n_obs+1)) {
    x_temp <- func_f(xt[j-1, ], par)
    xt[j, ] <- x_temp$y + noise_xt[j-1, ]
  }
  xt <- xt[-1, ]
  
  # Simulate time to maturities 
  TT <- seq(from = monthdays, to = n_contract*monthdays, by = monthdays)
  TT <- TT + 1
  mats <- matrix(0, nrow = n_obs, ncol = n_contract)
  for (j in 1: n_obs) {
    if ((j-1) %% monthdays == 0 & j != 1) {
      TT <- TT + monthdays
    }
    mats[j, ] <- (TT - j) / yeardays
  }
  
  # Simulate yt
  y_temp <- func_g(t(xt), par, mats)
  yt <- y_temp$y + noise_yt
  
  return(list(yt = yt, mats = mats, xt = xt))
}


