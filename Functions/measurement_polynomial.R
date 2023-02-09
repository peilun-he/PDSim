measurement_polynomial <- function(x, par, mats, degree, n_coe) {
  # Polynomial measurement equation
  # Inputs: 
  #   x: x_t
  #   par: a vector of parameters and model coefficients
  #   mats: maturities
  #   degree: degree of polynomial
  #   n_coe: the number of model coefficients
  # Outputs:
  #   y: y_t
  #   y_jacobian: Jacobian of y
  
  n_para <- length(par)
  
  if (n_coe != 0) {
    par_coe <- par[(n_para-n_coe+1): n_para] # model coefficients
    par <- par[1: (n_para-n_coe)] # model parameters
  }
  
  n_contract <- dim(mats)[2]
  n_point <- dim(x)[2] # number of points
  
  if (n_coe == (degree+1)*(degree+2)/2) {
    p_coordinate <- matrix(c(par_coe), ncol = 1)
  } else if (n_coe == 0) {
    p_coordinate <- matrix(taylor_coe(degree), ncol = 1)
  } else {
    stop("Incorrect number of coefficients. ")
  }
  
  G <- G_matrix(par, degree)
  y <- matrix(0, nrow = n_point, ncol = n_contract) 
  
  chi <- x[1, ]
  xi <- x[2, ]
  Hx <- rep(1, n_point)
  for (s in 1: degree) {
    for (i in seq(from = s, to = 0, by = -1)) {
      j <- s - i
      Hx <- rbind(Hx, chi^i * xi^j)
    }
  }
  Hx <- t(Hx)

  if (dim(mats)[1] == 1 & n_point == 1) {
    y_jacobian <- matrix(0, nrow = n_contract, ncol = 2) 
    for (j in 1: n_contract) {
      exp_matG <- decomposition_eigen(mats[, j] * G)
      exp_matG_p <- exp_matG %*% p_coordinate
      y[, j] <- Hx %*% exp_matG %*% p_coordinate
      for (s in 1: degree) {
        for (i in seq(from = s, to = 0, by = -1)) {
          k <- s - i
          if (i != 0) {
            y_jacobian[j, 1] <- y_jacobian[j, 1] + i * exp_matG_p[s*(s+1)/2+k+1] * chi^(i-1) * xi^k
          }
          if (k != 0) {
            y_jacobian[j, 2] <- y_jacobian[j, 2] + k * exp_matG_p[s*(s+1)/2+k+1] * chi^i * xi^(k-1)   
          }
        }
      }
    }
  } else if (dim(mats)[1] == 1 && n_point > 1) {
    y_jacobian <- 0
    for (j in 1: n_contract) {
      exp_matG <- decomposition_eigen(mats[, j] * G)
      y[, j] <- Hx %*% exp_matG %*% p_coordinate 
    }
  } else {
    y_jacobian <- 0
    for (i in 1: n_point) {
      for (j in 1: n_contract) {
        exp_matG <- decomposition_eigen(mats[i, j]*G)
        y[i, j] <- Hx[i, ] %*% exp_matG %*% p_coordinate
      }
    }
  }

  return(list(y = y, y_jacobian = y_jacobian))
}