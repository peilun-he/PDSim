simulate_data <- function(par, x0, n_obs, n_contract, func_f, func_g, n_coe, noise, seed) {
  # Simulate data, state variables and time to maturities.
  # Inputs:
  #   par: vector of parameters
  #   x0: initial values of state variables
  #   n_obs: number of observations
  #   n_contract: number of contracts
  #   func_f: function f(x), which should take two arguments, xt and a vector of parameters, 
  #     and return two values, f(x) and f'(x) 
  #     ( f'(x) is useless in this function, just make it consistent as used in other functions )
  #   func_g: function g(x), which should take three arguments, xt, a vector
  #       of parameters and maturities, and return two values, g(x) and g'(x) 
  #       ( g'(x) is useless in this function, just make it consistent as used in other functions ) 
  #   noise: Gaussian -> Gaussian noise for both state and measurement equations
  #   seed: seed for random values
  # Outputs: 
  #   yt: data
  #   mats: time to maturities
  #   xt: state variables
  
  require(MASS)
  
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
    
    V <- diag( par[9: (n_para-n_coe)] ^ 2 )
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


