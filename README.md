# Polynomial Diffusion Model Simulation and Estimation (V2.1)
## Introduction
Web application for the polynomial diffusion model. This app generates a futures price data by providing all parameters. Also, it gives state variables and contracts estimations through Extended Kalman Filter (EKF) or Unscented Kalman Filter (UKF). The Schwartz and Smith's two-factor model is also provided for comparison. 

Future plan: 
- Add decomposition of data through "seasonal" package 
- Add forecasting and smoothing 

## Installation
PDSim can be accessed in two ways: 

1. You can use PDSim on the Shiny server. This way, you don't need to have R installed on your computer. Just go to https://peilunhe.shinyapps.io/pdsim/ and use it there. 

2. Additionally, you can download and run PDSim locally, by running the following R code: 

```r
# install.packages("devtools") # uncomment if you do not have devtools installed
devtools::install_github("peilun-he/PDSim")
PDSim::run_app()
```

## User Guide

### Schwartz-Smith Model
The spot price $S_t$ is modelled as
$$\log(S_t) = \chi_t + \xi_t,$$
where $\chi_t$ represents the short-term fluctuation and $\xi_t$ is the long-term equilibrium price level. We assume both $\chi_t$ and $\xi_t$ follow an Ornstein–Uhlenbeck process
$$d\chi_t = - \kappa \chi_t dt + \sigma_{\chi} dW_t^{\chi}$$
$$d\xi_t = (\mu_{\xi} - \gamma \xi_t) dt + \sigma_{\xi} dW_t^{\xi}$$
for real processes and
$$d\chi_t = (- \kappa \chi_t - \lambda_{\chi}) dt + \sigma_{\chi} dW_t^{\chi*}$$
$$d\xi_t = (\mu_{\xi} - \gamma \xi_t - \lambda_{\xi}) dt + \sigma_{\xi} dW_t^{\xi*}$$
for risk-neutral processes. $\kappa, \gamma \in \mathbb{R}^+$ are speed of mean-reversion parameters; $\mu_{\xi} \in \mathbb{R}$ is the mean level of the long-term factor; $\sigma_{\chi}, \sigma_{\xi} \in \mathbb{R}^+$ are volatilities; $\lambda_{\chi}, \lambda_{\xi} \in \mathbb{R}$ are risk premiums; $W_t^{\chi*}$ and $W_t^{\xi*}$ are correlated standard Brownian Motions with correlation coefficient $\rho $. Under the arbitrage-free assumption, the futures price $F_{t,T}$ at current time $t$ with maturity time $T$ must be equal to the expected value of spot price at maturity time, i.e.,
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

### Some Other Hints
1. Once users enter all parameters, the data will be generated automatically. Users do NOT need to click any buttons. However, if users wish to generate more realisations under the same set of parameters, please click the 'Generate new data' button.

2. The seed to generate random numbers is fixed, i.e., for the same set of parameters, users will get exactly the same data every time they use PDSim.

3. Futures prices in all tables / plots are REAL prices (NOT the logarithm), no matter which model is used.

4. The 95% confidence interval is shown as a grey ribbon on each plot.

5. Because of the limitation of filtering methods, the standard error of the estimated futures price on the first day is extremely large. All plots of contracts estimation start from the second day.

## Version History 
**Version 2.1** (current version): 
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
