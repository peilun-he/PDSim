#' Unscented Kalman Filter
#' 
#' The Unscented Kalman Filter (UKF) is a variant of the Kalman Filter that is used
#' for estimating the hidden state variables of a non-linear dynamic system. 
#' 
#' Details: Eric Wan & Rudolph van der Merwe (2000). Non-linear state-space model: 
#' X_t = f(X_{t-1}) + epsilon_t, epsilon_t ~ N(0, W)
#' Y_t = g(X_t) + eta_t, eta_t ~ N(0, V)
#' 
#' @param par A vector of parameters. 
#' @param yt Futures price. 
#' @param mats Time to maturity. 
#' @param func_f Function `f`, which should take two arguments, xt and a vector of parameters, 
#' and return two values, the function value f(x) and the gredient f'(x). 
#' @param func_g Function `g`, which should take three arguments, xt, a vector of parameters 
#' and maturities, and return two values, g(x) and g'(x). 
#' @param dt \eqn{\Delta t}. The interval between two consecutive time points. 
#' @param n_coe The number of model coefficient to be estimated. 
#' @param noise The distribution of noise, `currently only "Gaussian" works.`Gaussian` or `Gamma`.  
#' 
#' @return This function returns a list with components: 
#' \item{nll}{Negative log-likelihood. }
#' \item{ll_table}{A vector to store cumulative log-likelihood at each time point. 
#' Used to calculate Sandwich variance. }
#' \item{table_xt_filter}{Filtered state variable. }
#' \item{table_xt_prediction}{Predicted state variable. }
#' 
#' @export
#' @seealso [KF], [EKF] for other filtering methods. 
#' @examples
#' ######################################
#' ##### Polynomial diffusion model #####
#' ######################################
#' n_obs <- 100
#' n_contract <- 10
#' par <- c(0.5, 0.3, 1, 1.5, 1.3, -0.3, 0.5, 0.3, seq(from = 0.1, to = 0.01, length.out = n_contract))
#' x0 <- c(0, 1/0.3)
#' dt <- 1/360 # daily data
#' n_coe <- 6 # polynomial with order 2`
#' par_coe <- c(1, 1, 1, 1, 1, 1)
#' func_f <- function(xt, par) state_linear(xt, par, dt) # state equation
#' func_g <- function(xt, par, mats) measurement_polynomial(xt, par, mats, 2, n_coe) # measurement equation 
#' dat <- simulate_data(c(par, par_coe), x0, n_obs, n_contract, func_f, func_g, n_coe, "Gaussian", 1234)
#' price <- dat$yt # measurement_polynomial function returns the futures price
#' mats <- dat$mats
#' est <- UKF(c(par, par_coe, x0), price, mats, func_f, func_g, dt, n_coe, "Gaussian")

UKF <- function(par, yt, mats, func_f, func_g, dt, n_coe, noise) {
  
  x0 <- c( par[length(par)-1], par[length(par)] )
  par <- par[1: (length(par)-2)]
  par_all <- par # model parameters and model coefficients
  par <- par[1: (length(par) - n_coe)] # model parameters
  
  kappa_chi  <- par[1]
  kappa_xi   <- par[2]
  mu_xi      <- par[3]
  sigma_chi  <- par[4]
  sigma_xi   <- par[5]
  rho        <- par[6]
  lambda_chi <- par[7]
  lambda_xi  <- par[8]
  
  n_obs <- dim(yt)[1]
  n_contract <- dim(yt)[2]
  n_state <- 2
  table_xt_filter <- matrix(0, nrow = n_obs, ncol = n_state)  # a_t|t
  table_Pt_filter <- array(0, dim = c(n_state, n_state, n_obs)) # P_t|t
  table_xt_prediction <- matrix(0, nrow = n_obs, ncol = n_state) # a_t|t-1
  table_Pt_prediction <- array(0, dim = c(n_state, n_state, n_obs)) # P_t|t-1
  table_Pyy <- array(0, dim = c(n_contract, n_contract, n_obs)) # Covariance of y
  nll <- 0 # negative log-likelihood
  ll_table <- matrix(0, nrow = 1, ncol = n_obs) 
  
  # Covariance matrices
  if (noise == "Gaussian") {
    if (n_contract == 1) {
      V <- par[length(par)]^2
    } else {
      V <- diag( par[9: length(par)]^2 )
    }
    
    if (kappa_xi != 0) {
      W <- matrix(c(sigma_chi^2/(2*kappa_chi) * ( 1-exp(-2*kappa_chi*dt) ), 
                    rho*sigma_chi*sigma_xi/(kappa_chi+kappa_xi) * ( 1-exp(-(kappa_chi+kappa_xi)*dt) ), 
                    rho*sigma_chi*sigma_xi/(kappa_chi+kappa_xi) * ( 1-exp(-(kappa_chi+kappa_xi)*dt) ), 
                    sigma_xi^2/(2*kappa_xi) * ( 1-exp(-2*kappa_xi*dt) )), 
                  nrow = 2, byrow = TRUE)
    } else if (kappa_xi == 0) {
      W <- matrix(c(sigma_chi^2/(2*kappa_chi) * ( 1-exp(-2*kappa_chi*dt) ), 
                    rho*sigma_chi*sigma_xi/(kappa_chi) * ( 1-exp(-(kappa_chi)*dt) ), 
                    rho*sigma_chi*sigma_xi/(kappa_chi) * ( 1-exp(-(kappa_chi)*dt) ), 
                    sigma_xi^2*dt), 
                  nrow = 2, byrow = TRUE)
    }
    
  } else if (noise == "Gamma") {
    s_sq <- par(9)
    V <- diag(rep(s_sq, n_contract))
    if (kappa_xi != 0) {
      W <- matrix(c(sigma_chi^2/(2*kappa_chi) * ( 1-exp(-2*kappa_chi*dt) ),
                    rho*sigma_chi*sigma_xi/(kappa_chi+kappa_xi) * ( 1-exp(-(kappa_chi+kappa_xi)*dt) ), 
                    rho*sigma_chi*sigma_xi/(kappa_chi+kappa_xi) * ( 1-exp(-(kappa_chi+kappa_xi)*dt) ), 
                    sigma_xi^2/(2*kappa_xi) * ( 1-exp(-2*kappa_xi*dt) )), 
                  nrow = 2, byrow = TRUE)
    } else if (kappa_xi == 0) {
      W <- matrix(c(sigma_chi^2/(2*kappa_chi) * ( 1-exp(-2*kappa_chi*dt) ), 
                    rho*sigma_chi*sigma_xi/(kappa_chi) * ( 1-exp(-(kappa_chi)*dt) ), 
                    rho*sigma_chi*sigma_xi/(kappa_chi) * ( 1-exp(-(kappa_chi)*dt) ), 
                    sigma_xi^2*dt), 
                  nrow = 2, byrow = TRUE)
    }
    
  } else {
    stop("Incorrect distribution of noises. ")
  }
  
  # Initialization
  # xt_filter <- c( 0, mu_xi / kappa_xi ) # x_0|0
  xt_filter <- x0
  if (kappa_xi != 0) {
    Pt_filter <- matrix(c(sigma_chi^2 / (2*kappa_chi), 
                          sigma_chi*sigma_xi*rho / (kappa_chi + kappa_xi), 
                          sigma_chi*sigma_xi*rho / (kappa_chi + kappa_xi), 
                          sigma_xi^2 / (2*kappa_xi)), 
                        nrow = 2, byrow = TRUE) # P_0|0
  } else if (kappa_xi == 0) {
    Pt_filter <- matrix(c(sigma_chi^2 / (2*kappa_chi), 
                          sigma_chi*sigma_xi*rho / (kappa_chi), 
                          sigma_chi*sigma_xi*rho / (kappa_chi), 
                          sigma_xi^2), 
                        nrow = 2, byrow = TRUE) # P_0|0
  }
  
  
  # Calculate weights of sigma points
  alpha <- 1
  beta <- 0
  kappa <- 0
  lambda <- alpha^2 * (n_state + kappa) - n_state
  weight_mean <- matrix(0, nrow = 1, ncol = 2*n_state+1) 
  weight_mean[1] <- lambda / (n_state + lambda)
  weight_mean[2: length(weight_mean)] <- 1 / ( 2*(n_state + lambda) )
  weight_cov <- weight_mean
  weight_cov[1] <- lambda / (n_state + lambda) + (1 - alpha^2 + beta)
  
  for (i in 1: n_obs) {
    # Calculate sigma points
    SP <- cbind( xt_filter, c(xt_filter) + t( chol((n_state+lambda) * Pt_filter) ), c(xt_filter) - t( chol((n_state+lambda) * Pt_filter) ) ) # each column forms a sigma point   
    
    # Prediction step
    temp <- func_f(SP, par) # 2*(2*n_state+1) matrix
    SPx_prediction <- temp$y
    xt_prediction <- SPx_prediction %*% t(weight_mean) # 2*1 matrix
    Pt_prediction <- W # 2*2 matrix
    for (j in 1: (2*n_state+1)) {
        Pt_prediction <- Pt_prediction + weight_cov[j] * (SPx_prediction[, j] - xt_prediction) %*% t(SPx_prediction[, j] - xt_prediction) 
    }
    temp <- func_g(SPx_prediction, par_all, matrix(mats[i, ], nrow = 1)) # (2*n_state+1)*n_contract matrix
    SPy_prediction <- temp$y
    yt_prediction <- weight_mean %*% SPy_prediction # 1*n_contract matrix

    # Filter step
    Pyy <- V # n_contract*n_contract matrix
    Pxy <- matrix(0, nrow = 2, ncol = n_contract) # 2*n_contract
    for (j in 1: (2*n_state+1)) {
      Pyy <- Pyy + weight_cov[j] * t(SPy_prediction[j, ] - yt_prediction) %*% (SPy_prediction[j, ] - yt_prediction)
      Pxy <- Pxy + weight_cov[j] * (SPx_prediction[, j] - xt_prediction) %*% (SPy_prediction[j, ] - yt_prediction)
    }
    K <- Pxy %*% solve(Pyy) # n_state * n_contract matrix
    et <- yt[i, ] - yt_prediction
    xt_filter <- xt_prediction + K %*% t(et)
    Pt_filter <- Pt_prediction - K %*% Pyy %*% t(K)
    #Pt_filter <- Pt_prediction - Pxy %*% t(K) - K %*% t(Pxy) + K %*% Pyy %*% t(K) # Joseph covariance update

    # Update tables
    table_xt_filter[i, ] <- xt_filter
    table_xt_prediction[i, ] <- xt_prediction
    table_Pt_filter[, , i] <- Pt_filter
    table_Pt_prediction[, , i] <- Pt_prediction
    table_Pyy[, , i] <- Pyy
    
    if (det(Pyy)<0) {
      message("matrix is not semi positive definite (EKF)")
    }
  
    # Update likelihood 
    nll <- nll + 0.5*length(yt[i, ])*log(2*pi) + 0.5*log(det(Pyy)) + 0.5*et %*% solve(Pyy) %*% t(et)
    ll_table[i] <- -(0.5*length(yt[i, ])*log(2*pi) + 0.5*log(det(Pyy)) + 0.5*et %*% solve(Pyy) %*% t(et))
  }
  
  return(list(nll = nll, 
              ll_table = ll_table, 
              xt_filter = table_xt_filter, 
              Pt_filter = table_Pt_filter, 
              xt_prediction = table_xt_prediction, 
              cov_y = table_Pyy))
}


