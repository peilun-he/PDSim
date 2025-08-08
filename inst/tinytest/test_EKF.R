n_obs <- 100 # number of observations
n_contract <- 10 # number of contracts
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

nll <- as.numeric(round(est_EKF$nll, 4))

expect_equal(nll, -1203.2297)

