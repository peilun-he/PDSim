n_obs <- 10 # number of observations
n_contract <- 2 # number of contracts
dt <- 1/360  # interval between two consecutive time points,

par <- c(0.5, 0.3, 1, 1.5, 1.3, -0.3, 0.5, 0.3,
         seq(from = 0.1, to = 0.01, length.out = n_contract)) # set of parameters
x0 <- c(0, 1/0.3) # initial values of state variables
n_coe <- 6 # number of model coefficient
par_coe <- c(1, 1, 1, 1, 1, 1) # model coefficients

# state equation
func_f <- function(xt, par) state_linear(xt, par, dt)
# measurement equation
func_g <- function(xt, par, mats) measurement_polynomial(xt, par, mats, 2, n_coe)

dat <- simulate_data(c(par, par_coe), x0, n_obs, n_contract,
                     func_f, func_g, n_coe, "Gaussian", 1234)
price <- dat$yt # measurement_polynomial function returns the futures price
mats <- dat$mats # time to maturity
xt <- dat$xt # state variables

est_UKF <- UKF(c(par, par_coe, x0), price, mats, func_f, func_g, dt, n_coe, "Gaussian")

results <- list(nll = as.numeric(round(est_UKF$nll, 4)), 
                xt_filter = round(est_UKF$xt_filter, 4))

exp_results <- list(nll = 16.8462, 
                    xt_filter = matrix(c(-0.6301, 3.1493, 
                                         -1.7531, 4.0916, 
                                         -1.8459, 4.1428, 
                                         -1.8514, 4.1444,
                                         -1.8576, 4.0808,
                                         -1.8903, 4.0938, 
                                         -1.9062, 4.1325,
                                         -1.9540, 4.2030, 
                                         -2.0007, 4.2702, 
                                         -2.0106, 4.1014), 
                                       nrow = 10, byrow = TRUE)) # expected results 

expect_equal(results, exp_results)
