# Polynomial Diffusion Model Simulation and Estimation (V2.1.2)

## Introduction

This package allows users to simulate commodity futures data from two models,
Schwartz and Smith two-factor model (Schwartz & Smith, 2000) and polynomial
diffusion model (Filipovic & Larsson, 2016), through both GUI and R scripts.
Additionally, it gives state variables and contract estimations through Kalman
Filter (KF), Extended Kalman Filter (EKF) or Unscented Kalman Filter (UKF).

Plans:
- Add decomposition of data through the "seasonal" package
- Add forecasting and smoothing

## Installation

PDSim can be accessed in two ways:

1. You can use PDSim on the Shiny server. This way, you don't need to have R
installed on your computer. Just go to <https://peilunhe.shinyapps.io/pdsim/>
and use it there.

2. Additionally, you can download and run PDSim locally, by running the
following R code:

```r
# install.packages("devtools") # uncomment if you do not have devtools installed
devtools::install_github("peilun-he/PDSim", build_vignettes = TRUE)
PDSim::run_app()
```

A tutorial of how to use this app is available by running the following code and select "PDSim app tutorial":
```r
browseVignettes("PDSim")
```

## How to Use PDSim (GUI)

The graphical user interface (GUI) is a easy way for everyone to use PDSim package, even though you have no knowledge of programming. Just enter all necessary parameters, it will simulate data, and provide well-designed interactive visualisations. Currently, PDSim can simulate data from two models, Schwartz and Smith two-factor model (Schwartz & Smith, 2000), and polynomial diffusion model (Filipovic & Larsson, 2016). In this section, we will explain how to use GUI to simulate data. A detailed description of two models are available in [Model Description](#model-description).

### Schwartz-Smith Model

Firstly, we establish certain global configurations, such as defining the number of observations (trading days) and contracts. Furthermore, we make a selection regarding the model from which the simulated data is generated.

<img src="figures/SS1.png" alt="drawing" width="400"/>

For Schwartz-Smith model (Schwartz & Smith, 2000), we assume the logarithm of spot price $S_t$, is the sum of two hidden factors:
$$\log{(S_t)} = \chi_t + \xi_t, $$
where $\chi_t$ represent the short term fluctuation and $\xi_t$ is the long term equilibrium price level. We assume both $\chi_t$ and $\xi_t$ follow a risk-neutral mean-reverting process:
$$d\chi_t = (- \kappa \chi_t - \lambda_{\chi}) dt + \sigma_{\chi} dW_t^{\chi*}, $$
$$d\xi_t = (\mu_{\xi} - \gamma \xi_t - \lambda_{\xi}) dt + \sigma_{\xi} dW_t^{\xi*}. $$
$\kappa, \gamma \in \mathbb{R}^+$ are called the speed of mean-reversion parameters, which controls how fast those two latent factors converge to their mean levels. Most of experiments suggest that $\kappa, \gamma \in (0, 3]$. In addition, we recommend that users limit $\kappa$ to be greater than $\gamma$, which means that the short-term fluctuation factor converges faster than the long-term factor. $\mu_{\xi} \in \mathbb{R}$ is the mean level of the long-term factor $\xi_t$. Here we assume that the short-term factor converges to 0. $\sigma_{\chi}, \sigma_{\xi} \in \mathbb{R}^+$ are volatilities, which represent the degree of variation of a price series over time. $\lambda_{\chi}, \lambda_{\xi} \in \mathbb{R}$ are risk premiums. We price commodity under arbitrage-free assumption: "the price of the derivative is set at the same level as the value of the replicating portfolio, so that no trader can make a risk-free profit by buying one and selling the other. If any arbitrage opportunities do arise, they quickly disappear as traders taking advantage of the arbitrage push the derivative’s price until it equals the value of replicating portfolios" ([Risk.Net](https://www.risk.net/definition/no-arbitrage-pricing), n.d., para. 2). In reality, we need to introduce the mean terms corrections - which are $\lambda_{\chi}$ and $\lambda_{\xi}$. $W_t^{\chi*}$ and $W_t^{\xi*}$ are correlated standard Brownian Motions with correlation coefficient $\rho$. In discrete time, these processes are discretised to Gaussian noises. All of them are parameters need to be specified. Additionally, for simplicity, we assume the standard errors $\sigma_i, i = 1, \dots, m$ for futures contracts are evenly spaced, i.e., $\sigma_1 - \sigma_2 = \sigma_2 - \sigma_3 = \dots = \sigma_{m-1} - \sigma_m$. If users don't know the value of a parameter, we recommend using the default value.

If users have special needs for the standard errors, please use R script.

<img src="figures/SS2.png" alt="drawing" width="400"/>

Finally, all the simulated data are downloadable. Please click `Download prices` and `Download maturities` buttons to download futures price and maturities data. Please note, even though Schwartz and Smith (2000) models the logarithm of spot price, **all data downloaded or plotted are real price, they have been exponentiated**. The other button `Generate new data` is designed for users who want to simulate multiple realisations from the same set of parameters. Once clicking it, PDSim will get another set of random noises, so the futures price will change as well. This button is not compulsory if users only need one realisations. The data will updated automatically when you change any parameters.

<img src="figures/SS3.png" alt="drawing" width="400"/>

### Polynomial Diffusion Model

The procedure for simulating data from the polynomial diffusion model (Filipovic & Larsson, 2016) closely resembles that of the Schwartz and Smith model (Schwartz & Smith, 2000). Nevertheless, it involves the specification of additional parameters.

Firstly, let's look at the difference between these two models. Both the polynomial diffusion model (Filipovic & Larsson, 2016) and the Schwartz and Smith model (Schwartz & Smith, 2000) assume that the spot price $S_t$ is influenced by two latent factors, $\chi_t$ and $\xi_t$. However, the Schwartz and Smith model (Schwartz & Smith, 2000) assumes that the logarithm of $S_t$ is the sum of two factors, whereas the polynomial diffusion model (Filipovic & Larsson, 2016) posits that $S_t$ takes on a polynomial form. Currently, PDSim GUI can only deal with polynomials with degree 2, i.e.,
$$S_t = \alpha_1 + \alpha_2 \chi_t + \alpha_3 \xi_t + \alpha_4 \chi_t^2 + \alpha_5 \chi_t \xi_t + \alpha_6 \xi_t^2. $$
$\alpha_i, i = 1, \dots, 6$ are extra parameters to the polynomial diffusion model (Filipovic & Larsson, 2016). If users want to specify a polynomial with degree 1, just set $\alpha_4 = \alpha_5 = \alpha_6 = 0$. Additionally, users are required to specify one non-linear filtering methods, Extended Kalman Filter (EKF) or Unscented Kalman Filter (UKF).

<img src="figures/PD1.png" alt="drawing" width="400"/>

All other procedures are the same as the Schwartz and Smith model (Schwartz & Smith, 2000).

### Some Other Hints

1. Once users enter all parameters, the data will be generated automatically. Users do NOT need to click any buttons. However, if users wish to generate more realisations under the same set of parameters, please click the 'Generate new data' button.

2. The seed to generate random numbers is fixed, i.e., for the same set of parameters, users will get exactly the same data every time they use PDSim.

3. Futures prices in all tables / plots are REAL prices (NOT the logarithm), no matter which model is used.

4. The 95% confidence interval is shown as a grey ribbon on each plot.

5. Because of the limitation of filtering methods, the standard error of the estimated futures price on the first day is extremely large. All plots of contracts estimation start from the second day.

## How to Use PDSim (R Script)

The GUI should be suffice. However, if you want to have more control of the data simulated, you can use R script. In this section, we will discuss how to use exported functions from this package to simulate data, as well as how to use Kalman Filter (KF), Extended Kalman Filter (EKF) and Unscented Kalman Filter (UKF) to estimate the hidden state variables.

Firstly, load the package:
```r
library(PDSim)
```
If you don't have PDSim installed, please refer [Installation](#installation).

Next, we specify the necessary global setups:
```r
n_obs <- 100 # number of observations
n_contract <- 10 # number of contracts
dt <- 1/360 # interval between two consecutive time points, 1/360 represents daily data
```

### Schwartz-Smith Model

Next, we specify parameters. For the Schwartz-Smith model (Schwartz & Smith, 2000), there is no model coefficients.
```r
par <- c(0.5, 0.3, 1, 1.5, 1.3, -0.3, 0.5, 0.3, seq(from = 0.1, to = 0.01, length.out = n_contract)) # set of parameters
x0 <- c(0, 1/0.3) # initial values of state variables
n_coe <- 0 # number of model coefficient
```
The set of parameters are in the order of: $\kappa, \gamma, \mu_{\xi}, \sigma_{\chi}, \sigma_{\xi}, \rho, \lambda_{\chi}, \lambda_{\xi}$. The long sequence is the standard errors of measurement equation. We assume all futures curves are uncorrelated and standard errors are evenly decreasing. You can have your own assumptions on standard errors, but the independence of curves must be hold.

Then, we specify the measurement and state equations. You can use the exported functions `measurement_linear` and `state_linear` directly, or write you own functions.
```r
func_f <- function(xt, par) state_linear(xt, par, dt) # state equation
func_g <- function(xt, par, mats) measurement_linear(xt, par, mats) # measurement equation
```

Finally, we can simulate the futures price, time to maturity, and hidden state variables:
```r
dat <- simulate_data(par, x0, n_obs, n_contract, func_f, measurement_linear, n_coe, "Gaussian", 1234)
log_price <- dat$yt # measurement_linear function returns the logarithm of futures price
mats <- dat$mats # time to maturity
xt <- dat$xt # state variables
```
Please note, `measurement_linear` returns the logarithm of futures price (which is required by the Schwartz and Smith model), so the data simulated is also the logarithm.

Additionally, we can estimate the hidden state variables through Kalman Filter (KF):
```r
est <- KF(par = c(par, x0), yt = log_price, mats = mats, delivery_time = 0, dt = dt, smoothing = FALSE, seasonality = "None") # delivery_time is unnecessary as we don't have seasonality
```

### Polynomial Diffusion Model

For the polynomial diffusion model (Filipovic & Larsson, 2016), we have to specify both parameters and model coefficients:
```r
par <- c(0.5, 0.3, 1, 1.5, 1.3, -0.3, 0.5, 0.3, seq(from = 0.1, to = 0.01, length.out = n_contract)) # set of parameters
x0 <- c(0, 1/0.3) # initial values of state variables
n_coe <- 6 # number of model coefficient
par_coe <- c(1, 1, 1, 1, 1, 1) # model coefficients
```
Currently, PDSim can deal with a polynomial with order 2, i.e., 6 model coefficients.

Then, we specify the measurement and state equations. Again, you can use the exported functions `state_linear` and `measurement_polynomial`.
```r
func_f <- function(xt, par) state_linear(xt, par, dt) # state equation
func_g <- function(xt, par, mats) measurement_polynomial(xt, par, mats, 2, n_coe) # measurement equation
```

Finally, simulate the data:
```r
dat <- simulate_data(c(par, par_coe), x0, n_obs, n_contract, func_f, func_g, n_coe, "Gaussian", 1234)
price <- dat$yt # measurement_polynomial function returns the futures price
mats <- dat$mats # time to maturity
xt <- dat$xt # state variables
```
`measurement_polynomial` returns the actual price, rather than the logarithm.

We can also estimate the hidden state variables through Extended Kalman Filter (EKF) or Unscented Kalman Filter (UKF):
```r
est_EKF <- EKF(c(par, par_coe, x0), price, mats, func_f, func_g, dt, n_coe, "Gaussian")
est_UKF <- UKF(c(par, par_coe, x0), price, mats, func_f, func_g, dt, n_coe, "Gaussian")
```

## Model Description

### Schwartz-Smith Model

The spot price $S_t$ is modelled as
$$\log(S_t) = \chi_t + \xi_t,$$
where $\chi_t$ represents the short-term fluctuation and $\xi_t$ is the long-term equilibrium price level. We assume both $\chi_t$ and $\xi_t$ follow an Ornstein–Uhlenbeck process
$$d\chi_t = - \kappa \chi_t dt + \sigma_{\chi} dW_t^{\chi}$$
$$d\xi_t = (\mu_{\xi} - \gamma \xi_t) dt + \sigma_{\xi} dW_t^{\xi}$$
for real processes and
$$d\chi_t = (- \kappa \chi_t - \lambda_{\chi}) dt + \sigma_{\chi} dW_t^{\chi*}$$
$$d\xi_t = (\mu_{\xi} - \gamma \xi_t - \lambda_{\xi}) dt + \sigma_{\xi} dW_t^{\xi*}$$
for risk-neutral processes. $\kappa, \gamma \in \mathbb{R}^+$ are speed of mean-reversion parameters; $\mu_{\xi} \in \mathbb{R}$ is the mean level of the long-term factor; $\sigma_{\chi}, \sigma_{\xi} \in \mathbb{R}^+$ are volatilities; $\lambda_{\chi}, \lambda_{\xi} \in \mathbb{R}$ are risk premiums; $W_t^{\chi*}$ and $W_t^{\xi*}$ are correlated standard Brownian Motions with correlation coefficient $\rho $. In the original Schwartz-Smith model (Schwartz & Smith, 2000), the parameter $\gamma$ is set to zero. However, in our extended model, we introduce the flexibility for this mean-reversion parameter associated with the long-term factor to take on non-zero values.

Under the arbitrage-free assumption, the futures price $F_{t,T}$ at current time $t$ with maturity time $T$ must be equal to the expected value of spot price at maturity time, i.e.,
$$\log(F_{t,T}) = \log(\mathbb{E}^\*(S_T \mathcal{F}_t | \mathcal{F}_t)),$$

where $F_t$ is the filtration until time $t$. In discrete time, we have the linear Gaussian state-space model:
$$x_t = c + E x_{t-1} + w_t,$$
$$y_t = d_t + F_t x_t + v_t,$$
where

$$
x_t = \left[ \begin{matrix}
\chi_t\\
\xi_t
\end{matrix} \right],
c = \left[ \begin{matrix} 0\\
\frac{\mu_{\xi}}{\gamma} (1 - e^{-\gamma \Delta t}) \end{matrix} \right],
E = \left[ \begin{matrix} e^{-\kappa \Delta t} & 0\\
0 & e^{-\gamma \Delta t} \end{matrix} \right],
F_t = \left[ \begin{matrix} e^{-\kappa (T_1 - t)}, \dots, e^{-\kappa (T_m - t)} \\
e^{-\gamma (T_1 - t)}, \dots, e^{-\gamma (T_m - t)} \end{matrix} \right]^\top,
$$

$$
y_t = \left( \log{(F_{t,T_1})}, \dots, \log{(F_{t,T_m})} \right)^\top,
d_t = \left( A(T_1 - t), \dots, A(T_m - t) \right)^\top,
$$

and $m$ is the number of contracts. The function $A(\cdot)$ is given by
$$A(t) = -\frac{\lambda_{\chi}}{\kappa}(1 - e^{-\kappa t}) + \frac{\mu_{\xi} - \lambda_{\xi}}{\gamma}(1 - e^{-\gamma t}) + \frac{1}{2} \left(\frac{1 - e^{-2\kappa t}}{2\kappa}\sigma_{\chi}^2 + \frac{1 - e^{-2\gamma t}}{2\gamma}\sigma_{\xi}^2 + 2\frac{1 - e^{-(\kappa + \gamma)t}}{\kappa + \gamma}\sigma_{\chi}\sigma_{\xi}\rho \right).$$

$w_t$ and $v_t$ are multivariate Gaussian noises with mean 0 and covariance matrix $\Sigma_w$ and $\Sigma_v$ respectively, where

$$
\Sigma_w = \left[ \begin{matrix}
\frac{1 - e^{-2\kappa \Delta t}}{2\kappa}\sigma_{\chi}^2 & \frac{1 - e^{-(\kappa + \gamma) \Delta t}}{\kappa + \gamma} \sigma_{\chi}\sigma_{\xi}\rho \\
\frac{1 - e^{-(\kappa + \gamma) \Delta t}}{\kappa + \gamma} \sigma_{\chi}\sigma_{\xi}\rho & \frac{1 - e^{-2\gamma \Delta t}}{2\gamma}\sigma_{\xi}^2
\end{matrix} \right],
$$

and we assume $\Sigma_v$ is diagonal, i.e.,

$$
\Sigma_v = \left[ \begin{matrix}
\sigma_1^2 & 0 & \dots & 0\\
0 & \sigma_2^2 & \dots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \dots & \sigma_m^2
\end{matrix} \right].
$$

Moreover, we assume $\sigma_1, \dots \sigma_m$ are evenly spaced, i.e., $\sigma_1 - \sigma_2 = \sigma_2 - \sigma_3 = \dots = \sigma_{m-1} - \sigma_m$.

### Polynomial Diffusion Model

Consider the stochastic differential equation
$$dX_t = b(X_t)dt + \sigma(X_t)dW_t,$$
where $W_t$ is a $d$-dimensional standard Brownian Motion and map $\sigma: \mathbb{R}^d \to \mathbb{R}^{d \times d}$ is continuous. Define $a := \sigma \sigma^\top$. For maps $a: \mathbb{R}^d \to \mathbb{S}^{d}$ and $b: \mathbb{R}^d \to \mathbb{R}^d$, suppose we have $a_{ij} \in Pol_2$ and $b_i \in Pol_1$. $\mathbb{S}^d$ is the set of all real symmetric $d \times d$ matrices and $Pol_n$ is the set of all polynomials of degree at most $n$. Then the solution of the SDE is a polynomial diffusion. Moreover, we define the generator $\mathcal{G}$ associated to the polynomial diffusion $X_t$ as
$$\mathcal{G}f(x) = \frac{1}{2} Tr\left( a(x) \nabla^2 f(x)\right) + b(x)^\top \nabla f(x)$$
for $x \in \mathbb{R}^d$ and any $C^2$ function $f$. Let $N$ be the dimension of $Pol_n$, and $H: \mathbb{R}^d \to \mathbb{R}^N$ be a function whose components form a basis of $Pol_n$. Then for any $p \in Pol_n$, there exists a unique vector $\vec{p} \in \mathbb{R}^N$ such that
$$p(x) = H(x)^\top \vec{p}$$
and $\vec{p}$ is the coordinate representation of $p(x)$. Moreover, there exists a unique matrix representation $G \in \mathbb{R}^{N \times N}$ of the generator $\mathcal{G}$, such that $G \vec{p}$ is the coordinate vector of $\mathcal{G} p$. So we have
$$\mathcal{G} p(x) = H(x)^\top G \vec{p}.$$
<br/><br/>
**Theorem 1**: Let $p(x) \in Pol_n$ be a polynomial with coordinate representation $\vec{p} \in \mathbb{R}^N$, $G \in \mathbb{R}^{N \times N}$ be a matrix representation of generator $\mathcal{G}$, and $X_t \in \mathbb{R}^d$ satisfies the SDE. Then for $0 \le t \le T$, we have
$$\mathbb{E} \left[ p(X_T) | \mathcal{F}_t \right] = H(X_t)^\top e^{(T-t)G} \vec{p},$$
where $\mathcal{F}_t$ represents all information available until time $t$.
<br/><br/><br/>

Next, we apply this theorem to the two-factor model. Assume the spot price $S_t$ is modelled as
$$S_t = p_n(x_t) = \alpha_1 + \alpha_2 \chi_t + \alpha_3 \xi_t + \alpha_4 \chi_t^2 + \alpha_5 \chi_t \xi_t + \alpha_6 \xi_t^2,$$
which is a polynomial with order $n = 2$. $x_t = (\chi_t, \xi_t)^\top$ is a vector of state variables. Obviously, $x_t$ satisfies the SDE with

$$
b(x_t) = \left[ \begin{matrix}
-\kappa \chi_t - \lambda_{\chi} \\
\mu_{\xi} - \gamma \xi_t - \lambda_{\xi}
\end{matrix} \right],
\sigma(x_t) = \left[ \begin{matrix}
\sigma_{\chi} & 0 \\
0 & \sigma_{\xi}
\end{matrix} \right],
a(x_t) = \sigma(x_t) \sigma(x_t)^\top = \left[ \begin{matrix}
\sigma_{\chi}^2 & 0 \\
0 & \sigma_{\xi}^2
\end{matrix} \right].
$$

The basis
$$H(x_t) = (1, \chi_t, \xi_t, \chi_t^2, \chi_t \xi_t, \xi_t^2)^\top$$
has a dimension $N = 6$. The coordinate representation is
$$\vec{p} = (\alpha_1, \alpha_2, \alpha_3, \alpha_4, \alpha_5, \alpha_6)^\top.$$
By applying $\mathcal{G}$ to each element of $H(x_t)$, we get the matrix representation

$$
G = \left[ \begin{matrix}
0 & -\lambda_{\chi} & \mu_{\xi} - \lambda_{\xi} & \sigma_{\chi}^2 & 0 & \sigma_{\xi}^2 \\
0 & -\kappa & 0 & -2 \lambda_{\chi} & \mu_{\xi} - \lambda_{\xi} & 0 \\
0 & 0 & -\gamma & 0 & -\lambda_{\chi} & 2\mu_{\xi} - 2\lambda_{\xi} \\
0 & 0 & 0 & -2\kappa & 0 & 0 \\
0 & 0 & 0 & 0 & -\kappa - \gamma & 0 \\
0 & 0 & 0 & 0 & 0 & -2\gamma
\end{matrix} \right].
$$

Then, by Theorem 1, the futures price $F_{t,T}$ is given by
$$F_{t,T} = \mathbb{E}^*(S_T | \mathcal{F}_t) = H(x_t)^\top e^{(T-t)G} \vec{p}.$$

Therefore, we have the non-linear state-space model
$$x_t = c + E x_{t-1} + w_t,$$
$$y_t = H(x_t)^\top e^{(T-t)G} \vec{p} + v_t.$$

## Contributions and Supports

If you find any bugs or want to make a contribution to this package, please create a GitHub issue at: <https://github.com/peilun-he/PDSim/issues>.

Additionally, you are very welcome to provide any kind of feedback and comments. Please send me an email at: peilun.he93\@gmail.com.

If you have questions about how to use this package, please also send me an email. I will get back to you as soon as possible.

## Acknowledgements

We would like to thank Sam Forbes, Blake Rayfield and Mark Van de Vyver for testing PDSim and providing valuable feedback and suggestions.

## Version History

**Version 2.1.2** (current version):
- Main functions are exported, with short executable examples.
- Add Contributions and Supports section.

**Version 2.1.1**:
- Add a vignette.

**Version 2.1**:
- PDSim is packaged into an R package. Some structures is changed to achieve this.
- A exported function "run_app" is added to run PDSim.
- Add some documentation.  

**Version 2.0**:
- Add navigation bar: welcome page, app, user guide, team members.
- Descriptions of models and some hints are added to the user guide page.
- Allow users to download simulated data as csv files.
- Add 95% confidence intervals to all estimations.  
- Add a 3D surface of data.
- Allow users to generate new realisations of data using same set of parameters.
- Bugs fixed.

**Version 1.0**: basic functions and UI

## References

Aspinall, T., Gepp, A., Harris, G., Kelly, S., Southam, C., & Vanstone, B. (2022). NFCP: N-factor commodity pricing through term structure estimation. *The Comprehensive R Archive Network*. [https://cran.rstudio.com/web/packages/NFCP/index.html](https://cran.rstudio.com/web/packages/NFCP/index.html).

Filipovic, D., & Larsson, M. (2016). Polynomial diffusions and applications in finance. *Finance and Stochastics*, 20(4), 931–972.

Harvey, A. C. (1990). Forecasting, structural time series models and the kalman filter. *Cambridge University Press*.

Julier, S. J., & Uhlmann, J. K. (1997). New extension of the kalman filter to nonlinear systems. *Signal Processing, Sensor Fusion, and Target Recognition VI*, 3068, 182–193.

Julier, S. J., & Uhlmann, J. K. (2004). Unscented filtering and nonlinear estimation. *Proceedings of the IEEE*, 92(3), 401–422.

Kleisinger-Yu, X., Komaric, V., Larsson, M., & Regez, M. (2020). A multifactor polynomial framework for long-term electricity forwards with delivery period. *SIAM Journal on Financial Mathematics*, 11(3), 928–957.

Risk.net. (n.d.). No arbitrage pricing. Retrieved from <https://www.risk.net/definition/no-arbitrage-pricing>.

Schwartz, E. S., & Smith, J. E. (2000). Short-term variations and long-term dynamics in commodity prices. *Management Science*, 46(7), 893–911.

Wan, E. A., & Van Der Merwe, R. (2000). The unscented kalman filter for nonlinear estimation. *Proceedings of the IEEE 2000 Adaptive Systems for Signal Processing, Communications, and Control Symposium (Cat. No. 00EX373)*, 153–158.
