n_obs <- 10 # number of observations
n_contract <- 2 # number of contracts
dt <- 1/360  

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

est_EKF <- EKF(c(par, par_coe, x0), price, mats, func_f, func_g, dt, n_coe, "Gaussian")

results <- list(nll = as.numeric(round(est_EKF$nll, 4)), 
                xt_filter = round(est_EKF$xt_filter, 4))

exp_results <- list(nll = 0.7012, 
                    xt_filter = matrix(c(0.0676, 3.3209, 
                                         0.3172, 3.2577, 
                                         0.4222, 3.2507, 
                                         0.2630, 3.3284,
                                         0.3354, 3.2214,
                                         0.4729, 3.1481, 
                                         0.3595, 3.2505,
                                         0.4387, 3.2624, 
                                         0.4499, 3.3103, 
                                         0.4464, 3.1427), 
                                       nrow = 10, byrow = TRUE)) # expected results 

expect_equal(results, exp_results)

