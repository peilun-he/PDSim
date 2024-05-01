#' Standard Kalman Filter and Smoother
#' 
#' An algorithm for estimating the hidden state variables of a linear dynamic system. 
#' 
#' Details. 
#' 
#' @param par A vector of parameters. 
#' @param yt The logarihm of futures prices. 
#' @param T Time to maturity. 
#' @param delivery_time A vector of date, which is necessary if `seasonality` is `Constant`.
#' @param dt \eqn{\Delta t}. The interval between two consecutive time points.  
#' @param smoothing Boolean. Indicate if Kalman Smoothing is required. 
#' @param seasonality `Constant` or `None`. 
#' 
#' @return This function returns a list with components: 
#' \item{nll}{The negative log likelihood. }
#' \item{ll_table}{A vector to store cumulative log-likelihood at each time point. 
#' Used to calculate Sandwich variance. }
#' \item{table_at_filter}{A matrix gives the filtered values of state variables. }
#' \item{table_at_prediction}{A matrix gives the predicted values of state variables. }
#' \item{table_at_smoother}{A matrix gives the smoothed values of state variables. 
#' The algorithm of Kalman Smoother is given by Bierman (1973) and De Jong (1989). }
#' \item{ft}{Seasonal effect. }
#' 
#' @import lubridate
#' @export
#' @seealso [EKF], [UKF] for other filtering methods. 
#' @examples
#' ###############################################
#' ##### Schwartz and Smith two-factor model #####
#' ###############################################
#' n_obs <- 100
#' n_contract <- 10
#' par <- c(0.5, 0.3, 1, 1.5, 1.3, -0.3, 0.5, 0.3, seq(from = 0.1, to = 0.01, length.out = n_contract))
#' x0 <- c(0, 1/0.3)
#' dt <- 1/360 # daily data
#' n_coe <- 0
#' func_f <- function(xt, par) state_linear(xt, par, dt) # state equation
#' dat <- simulate_data(par, x0, n_obs, n_contract, func_f, measurement_linear, n_coe, "Gaussian", 1234)
#' log_price <- dat$yt # measurement_linear function returns the logarithm of futures price
#' mats <- dat$mats
#' est <- KF(par = c(par, x0), yt = log_price, mats = mats, delivery_time = 0, # delivery_time is unnecessary 
#'           dt = dt, smoothing = FALSE, seasonality = "None")


KF <- function(par, yt, mats, delivery_time, dt, smoothing, seasonality) {
  
  x0 <- c( par[length(par)-1], par[length(par)] )
  par <- par[1: (length(par)-2)]
  
  n_obs <- dim(yt)[1]
  n_contract <- dim(yt)[2]
  
  table_xt_filter <- matrix(0, nrow = n_obs, ncol = 2) # a_t|t
  table_Pt_filter <- array(0, dim = c(2, 2, n_obs)) # P_t|t
  table_xt_prediction <- matrix(0, nrow = n_obs+1, ncol = 2) # a_t|t-1
  table_Pt_prediction <- array(0, dim = c(2, 2, n_obs+1)) # P_t|t-1
  table_xt_smoother <- matrix(0, nrow = n_obs, ncol = 2) # a_t|s, s > t
  table_Pt_smoother <- array(0, dim = c(2, 2, n_obs)) # P_t|s, s > t
  
  nll <- 0 # negative log-likelihood 
  ll_table <- matrix(0, nrow = 1, ncol = n_obs) # table of log-likelihood
  
  table_et <- matrix(0, nrow = n_obs, ncol = n_contract) # e_t
  table_invL <- array(0, dim = c(n_contract, n_contract, n_obs)) # inverse of L_t|t-1
  table_K <- array(0, dim = c(2, n_contract, n_obs)) # Kalman gain matrix
  table_D <- matrix(0, nrow = n_obs, ncol = n_contract) # d_t
  table_E <- array(0, dim = c(2, n_contract, n_obs)) # F_t
  table_L <- array(0, dim = c(n_contract, n_contract, n_obs)) # Covariance of y
  
  # Seasonal component
  if (seasonality == "Constant") {
    par_seasonal <- par[(length(par)-11): length(par)]
    ft <- matrix(0, nrow = n_obs, ncol = n_contract)
    for (i in 1: n_obs) {
      ft[i, ] <- par_seasonal * t( (month( t(delivery_time[i]) ) == 1: 12) )
    }
    par <- par[1: (length(par)-12)] 
  } else if (seasonality == "None") {
    ft <- rep(0, n_obs)
  } else {
    stop("The seasonal component must be 'Constant' or 'None'. ")
  }
  
  # Parameters
  kappa_chi  <- par[1]
  kappa_xi   <- par[2]
  mu_xi      <- par[3]
  sigma_chi  <- par[4]
  sigma_xi   <- par[5]
  rho        <- par[6]
  lambda_chi <- par[7]
  lambda_xi  <- par[8]
  
  # Covariance matrices
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
  
  
  # Initialization
  # xt_filter <- c( 0, mu_xi / kappa_xi ) # x_0|0
  xt_filter <- x0
  if (kappa_xi != 0) {
    Pt_filter <- matrix(c(sigma_chi^2 / (2*kappa_chi), 
                          sigma_chi*sigma_xi*rho / (kappa_chi + kappa_xi), 
                          sigma_chi*sigma_xi*rho / (kappa_chi + kappa_xi), 
                          sigma_xi^2 / (2*kappa_xi)), 
                        nrow = 2, byrow = TRUE) # P_0|0
    
    C <- c( 0, mu_xi/kappa_xi*(1-exp(-kappa_xi*dt)) )
    G <- matrix(c( exp(-kappa_chi*dt), 0, 0, exp(-kappa_xi*dt) ), 
                nrow = 2, byrow = TRUE)
  } else if (kappa_xi == 0) {
    Pt_filter <- matrix(c(sigma_chi^2 / (2*kappa_chi), 
                          sigma_chi*sigma_xi*rho / (kappa_chi), 
                          sigma_chi*sigma_xi*rho / (kappa_chi), 
                          sigma_xi^2), 
                        nrow = 2, byrow = TRUE) # P_0|0
    
    C <- c( 0, mu_xi*dt )
    G <- matrix(c( exp(-kappa_chi*dt), 0, 0, 1 ), nrow = 2, byrow = TRUE)
  }
  
  
  
  # Kalman Filter
  for (i in 1:n_obs) {   
    D <- AofT(par, mats[i,]) + ft[i] # d_t + f_t
    E <- t (rbind( exp(-kappa_chi*mats[i, ]), exp(-kappa_xi*mats[i, ]) ) ) # E_t
    
    # Prediction step
    xt_prediction  <- C + G %*% xt_filter # a_t+1|t 
    Pt_prediction <- G %*% Pt_filter %*% t(G) + W # P_t+1|t
    y_prediction <- D + E %*% xt_prediction # ytilde_t|t-1 = d_t + F_t a_t|t-1
    
    # Filter step
    et <- yt[i, ] - t(y_prediction) # e_t = y_t - ytilde_t|t-1
    L <- E %*% Pt_prediction %*% t(E) + V # Covariance matrix of et
    invL <- solve(L) # inverse of L 
    K <- Pt_prediction %*% t(E) %*% invL # Kalman gain matrix: K_t
    
    xt_filter <- xt_prediction + K %*% t(et) # a_t
    Pt_filter <- (diag(2) - K %*% E) %*% Pt_prediction # P_t
    #Pt_filter <- (diag(2) - K %*% E) %*% Pt_prediction %*% t(diag(2) - K %*% E) + K %*% V %*% t(K) 
           
    # Update tables
    table_xt_filter[i, ] <- t(xt_filter)
    table_Pt_filter[, , i] <- Pt_filter
    table_xt_prediction[i+1, ] <- t(xt_prediction)
    table_Pt_prediction[, , i+1] <- Pt_prediction
    table_et[i, ] <- et
    table_invL[, , i] <- invL 
    table_K[, , i] <- K
    table_D[i, ] <- t(D)
    table_E[, , i] <- t(E)
    table_L[, , i] <- L
    
    if (det(L)<0) {
      message("matrix is not semi positive definite (KF)")
    }
    
    # Update likelihood 
    nll <- nll + 0.5*length(yt[i, ])*log(2*pi) + 0.5*log(det(L)) + 0.5*et %*% solve(L) %*% t(et)
    ll_table[i] <- -(0.5*length(yt[i, ])*log(2*pi) + 0.5*log(det(L)) + 0.5*et %*% solve(L) %*% t(et))
  }
  
  # Kalman Smoother
  if (smoothing) {
    for (t in seq(from = n_obs, to = 1, by = -1)) {
      E <- t(table_E[, , t])
      D <- t(table_D[t, ])
      K <- table_K[, , t]
      invL <- table_invL[, , t]
      et <- t(table_et[t, ])
            
      xt_prediction <- t(table_xt_prediction[t, ])
      Pt_prediction <- table_Pt_prediction[, , t]
      
      if (t == n_obs) {
        rt <- matrix(0, nrow = 2, ncol = 1) 
        Rt <- matrix(0, nrow = 2, ncol = 2)
      }
      
      rt <- t(E) %*% invL %*% et + t(G - G %*% K %*% E) %*% rt # 2 * 1 matrix 
      Rt <- t(E) %*% invL %*% E + t(G - G %*% K %*% E) %*% Rt %*% (G - G %*% K %*% E) # 2 * 2 matrix
      
      xt_smoother <- xt_prediction + Pt_prediction %*% rt # a_t|n
      Pt_smoother <- Pt_prediction - Pt_prediction %*% Rt %*% Pt_prediction # P_t|n
      
      # Update tables
      table_xt_smoother[t, ] <- xt_smoother
      table_Pt_smoother[, , t] <- Pt_smoother
    }    
  } else {
    table_xt_smoother <- 0
  }
  
  return(list(nll = nll, 
              ll_table = ll_table, 
              xt_filter = table_xt_filter,
              Pt_filter = table_Pt_filter, 
              xt_prediction = table_xt_prediction, 
              cov_y = table_L, 
              xt_smoother = table_xt_smoother, 
              seasonal = ft))
}