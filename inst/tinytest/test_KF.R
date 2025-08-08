n_obs <- 10 # number of observations
n_contract <- 2 # number of contracts
dt <- 1/360 

par <- c(0.5, 0.3, 1, 1.5, 1.3, -0.3, 0.5, 0.3,
         seq(from = 0.1, to = 0.01, length.out = n_contract)) # set of parameters
x0 <- c(0, 1/0.3) # initial values of state variables
n_coe <- 0 # number of model coefficient

# state equation
func_f <- function(xt, par) state_linear(xt, par, dt)
# measurement equation
func_g <- function(xt, par, mats) measurement_linear(xt, par, mats)

dat <- simulate_data(par, x0, n_obs, n_contract,
                     func_f, func_g, n_coe, "Gaussian", 1234)
log_price <- dat$yt # logarithm of futures price
mats <- dat$mats # time to maturity
xt <- dat$xt # state variables

est <- KF(par = c(par, x0), yt = log_price, mats = mats,
          delivery_time = 0, dt = dt, smoothing = FALSE,
          seasonality = "None")

results <- list(nll = as.numeric(round(est$nll, 4)), 
               xt_filter = round(est$xt_filter, 4))

exp_results <- list(nll = -19.6466, 
                   xt_filter = matrix(c(0.0214, 3.3697, 
                                        0.1832, 3.3044, 
                                        0.2870, 3.2443, 
                                        0.1984, 3.3786,
                                        0.2675, 3.2363,
                                        0.5356, 2.9677, 
                                        0.3988, 3.1637,
                                        0.6151, 3.0347, 
                                        0.6439, 3.0774, 
                                        0.6983, 2.8524), 
                                      nrow = 10, byrow = TRUE)) # expected results 

expect_equal(results, exp_results)

