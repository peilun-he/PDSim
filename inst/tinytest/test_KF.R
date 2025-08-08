n_obs <- 100 # number of observations
n_contract <- 10 # number of contracts
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

nll <- as.numeric(round(est$nll, 4))

expect_equal(nll, -1459.0542)

