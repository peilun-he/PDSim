---
title: 'PDSim: A Shiny App for Polynomial Diffusion Model Simulation and Estimation'
tags:
  - Shiny
  - finance
  - dynamics
  - futures pricing
authors:
  - name: Peilun He
    orcid: 0000-0002-2740-3390
    equal-contrib: true
    affiliation: 1
  - name: Nino Kordzakhia
    orcid: 0000-0002-7853-4550
    equal-contrib: true
    affiliation: 2 
  - name: Gareth W. Peters
    orcid: 0000-0003-2768-8979
    equal-contrib: true
    affiliation: 3
  - name: Pavel V. Shevchenko
    orcid: 0000-0001-8104-8716
    equal-contrib: true
    affiliation: 1
affiliations:
 - name: Department of Actuarial Studies and Business Analytics, Macquarie University, Australia
   index: 1
 - name: School of Mathematical and Physical Sciences, Macquarie University, Australia
   index: 2
 - name: Department of Statistics and Applied Probability, University of California Santa Barbara, USA
   index: 3
date: 09 February 2023
bibliography: paper.bib

---

# Summary

The Schwartz-Smith two-factor model [@schwartz:2000] was commonly used in the pricing of commodity futures in the last two decades. In 2016, [@filipovic:2016] introduced a polynomial diffusion framework which allows a more complex structure of spot price. This framework has been applied to electricity forwards [@kleisinger:2020], in which the spot price is modelled in a quadratic form of two factors. PDSim aims to estimate futures prices as well as the latent state variables, and provides well-designed visualisations. This application is available at [https://github.com/peilun-he/PDSim](https://github.com/peilun-he/PDSim). 

## Schwartz-Smith two-factor model

Under the Schwartz-Smith framework, the logarithm of spot price $S_t$ is modelled as the sum of two factors $\chi_t$ and $\xi_t$, 
\begin{equation}\label{eq:st}
\log{(S_t)} = \chi_t + \xi_t,
\end{equation}
where $\chi_t$ represents the short-term fluctuation and $\xi_t$ is the long-term equilibrium price level. Additionally, we assume both $\chi_t$ and $\xi_t$ follow a risk-neutral Ornstein-Uhlenbeck process, 
\begin{equation}
d\chi_t = (-\kappa \chi_t - \lambda_{\chi}) dt + \sigma_{\chi} d W_t^{\chi}, 
\label{eq:SS_rn_chi}
\end{equation}
and 
\begin{equation}
d\xi_t = (\mu_{\xi} - \gamma \xi_t - \lambda_{\xi})dt + \sigma_{\xi} dW_t^{\xi},
\label{eq:SS_rn_xi}
\end{equation}
where $\kappa, \gamma \in \mathbb{R}^+$ are the speed of mean-reversion parameters, $\mu_{\xi} \in \mathbb{R}$ is the mean level of the long-term factor, $\sigma_{\chi}, \sigma_{\xi} \in \mathbb{R}^+$ are the volatility parameters, and $\lambda_{\chi}, \lambda_{\xi} \in \mathbb{R}$ are risk premiums. The processes $(W_t^{\chi})_{t \ge 0}$ and $(W_t^{\xi})_{t \ge 0}$ are correlated standard Brownian Motions with 
$$\mathbb{E} \left(dW_t^{\chi} dW_t^{\xi}\right) = \rho dt. $$ 
We set $\lambda_{\chi} = \lambda_{\xi} = 0$ in \autoref{eq:SS_rn_chi} and \autoref{eq:SS_rn_xi} to get the real processes. We use the risk-neutral processes for futures pricing, and real processes for modelling state variables. 

In discrete time, $\chi_t$ and $\xi_t$ are jointly normally distributed. Therefore, the spot price is log-normally distributed. Moreover, under the arbitrage-free assumption, the futures price ($F_{t,T}$) at current time $t$ must be equal to the expected value of spot price at maturity time $T$, 
$$F_{t,T} = \mathbb{E}^*(S_T | \mathcal{F}_t), $$
where $\mathcal{F}_t$ is a natural filtration and $\mathbb{E}^*(\cdot)$ is the expectation under the risk-neutral processes from \autoref{eq:SS_rn_chi} and \autoref{eq:SS_rn_xi}. Then we can get the linear Gaussian state space model: 
\begin{equation}
x_t = c + Ex_{t-1} + w_t, 
\end{equation}
\begin{equation}
y_t = d_t + F_t x_t + v_t, 
\end{equation}
where $x_t = \left[ \begin{matrix} \chi_t \\ \xi_t \end{matrix} \right], c = \left[ \begin{matrix} 0 \\ \frac{\mu_{\xi}}{\gamma} \left(1 - e^{-\gamma \Delta t} \right) \end{matrix} \right], E = \left[ \begin{matrix} e^{-\kappa \Delta t} & 0 \\ 0 & e^{-\gamma \Delta t}\end{matrix} \right], y_t = \left( \log{(F_{t,T_1})}, \dots, \log{(F_{t,T_m})} \right)^\top, d_t = \left( A(T_1 - t), \dots, A(T_m - t) \right)^\top, F_t = \left[ \begin{matrix} e^{-\kappa (T_1 - t)}, \dots, e^{-\kappa (T_m - t)} \\ e^{-\gamma (T_1 - t)}, \dots, e^{-\gamma (T_m - t)} \end{matrix} \right]^\top$ and $m$ is the number of futures contracts. The function $A(\cdot)$ is given by 
\begin{align}
A(t) =& -\frac{\lambda_{\chi}}{\kappa}(1 - e^{-\kappa t}) + \frac{\mu_{\xi} - \lambda_{\xi}}{\gamma}(1 - e^{-\gamma t}) \nonumber \\
&+ \frac{1}{2}\left(\frac{1 - e^{-2\kappa t}}{2\kappa}\sigma_{\chi}^2 + \frac{1 - e^{-2\gamma t}}{2\gamma}\sigma_{\xi}^2 + 2\frac{1 - e^{-(\kappa + \gamma)t}}{\kappa + \gamma}\sigma_{\chi}\sigma_{\xi}\rho \right). \nonumber
\end{align}
$w_t$ and $v_t$ are multivariate Gaussian noises with mean $\textbf{0}$ and covariance matrix $\Sigma_w$ and $\Sigma_v$ respectively, where 
$$\Sigma_w = \left[ \begin{matrix}
\frac{1 - e^{-2\kappa \Delta t}}{2\kappa}\sigma_{\chi}^2 & \frac{1 - e^{-(\kappa + \gamma) \Delta t}}{\kappa + \gamma}\sigma_{\chi}\sigma_{\xi}\rho \\
\frac{1 - e^{-(\kappa + \gamma) \Delta t}}{\kappa + \gamma}\sigma_{\chi}\sigma_{\xi}\rho & \frac{1 - e^{-2\gamma \Delta t}}{2\gamma}\sigma_{\xi}^2
\end{matrix} \right],$$ 
and we assume $\Sigma_v$ is diagonal, $\Sigma_v = diag(\sigma_1^2, \sigma_2^2, \dots, \sigma_m^2)$. 
Under this framework, $c, E, \Sigma_w$ and $\Sigma_v$ are deterministic but $d_t$ and $F_t$ are time-variant. 

## Polynomial diffusion model

In this section, we present a general framework of the polynomial diffusion model first, and then we give the application in the two-factor model. The mathematical foundations and applications of the polynomial diffusion model in finance are provided in [@filipovic:2016].

Consider the stochastic differential equation
\begin{equation}
dX_t = b(X_t)dt + \sigma(X_t)dW_t, 
\label{eq:sde}
\end{equation}
where $W_t$ is a $d$-dimensional standard Brownian motion and map $\sigma: \mathbb{R}^d \to \mathbb{R}^{d \times d}$ is continuous. Define $a := \sigma \sigma^\top$. For maps $a: \mathbb{R}^d \to \mathbb{S}^{d}$ and $b: \mathbb{R}^d \to \mathbb{R}^d$, suppose we have 
$a_{ij} \in Pol_2$ and $b_i \in Pol_1$. $\mathbb{S}^d$ is the set of all real symmetric $d \times d$ matrices and $Pol_n$ is the set of all polynomials of degree at most $n$. Then the solution of \autoref{eq:sde} is a polynomial diffusion. 

Moreover, we define the generator $\mathcal{G}$ associated to the polynomial diffusion $X_t$ as
\begin{equation}
\mathcal{G}f(x) = \frac{1}{2} Tr\left( a(x) \nabla^2 f(x)\right) + b(x)^\top \nabla f(x)
\label{eq:generator}
\end{equation}
for $x \in \mathbb{R}^d$ and any $f \in C^2$, twice continuously differentiable functions. Let $N$ be the dimension of $Pol_n$, and $H: \mathbb{R}^d \to \mathbb{R}^N$ be a function whose components form a basis of $Pol_n$. Then for any $p \in Pol_n$, there exists a unique vector $\vec{p} \in \mathbb{R}^N$ such that 
\begin{equation}
p(x) = H(x)^\top \vec{p}
\label{eq:vec_p}
\end{equation}
and $\vec{p}$ is the coordinate representation of $p(x)$. Moreover, there exists a unique matrix representation $G \in \mathbb{R}^{N \times N}$ of the generator $\mathcal{G}$, such that $G \vec{p}$ is the coordinate vector of $\mathcal{G} p$. Hence, we have
\begin{equation}
\mathcal{G} p(x) = H(x)^\top G \vec{p}. 
\label{eq:G}
\end{equation}

Theorem 1: Let $p(x) \in Pol_n$ be a polynomial with coordinate representation $\vec{p} \in \mathbb{R}^N$ satisfying \autoref{eq:vec_p}, $G \in \mathbb{R}^{N \times N}$ be a matrix representation of generator $\mathcal{G}$ satisfying \autoref{eq:G}, and $X_t \in \mathbb{R}^d$ satisfies \autoref{eq:sde}. Then for $0 \le t \le T$, we have
	$$\mathbb{E} \left(p(X_T)|\mathcal{F}_t \right) = H(X_t)^\top e^{(T-t)G} \vec{p}, $$
	where $\mathcal{F}_t$ is a natural $\sigma$-algebra generated up tp time $t$. 
	\label{th:pd}

The proof of Theorem 1 is given in [@filipovic:2016]. 

Next, we apply this theorem to the two-factor model. Assume the spot price $S_t$ is modelled as
\begin{equation}
S_t = \alpha_1 + \alpha_2 \chi_t + \alpha_3 \xi_t + \alpha_4 \chi_t^2 + \alpha_5 \chi_t \xi_t + \alpha_6 \xi_t^2, 
\label{eq:poly_st}
\end{equation}
where $S_t$ is a polynomial function with a degree $n = 2$ and $x_t = (\chi_t, \xi_t)^\top$ is a vector of state variables with $\chi_t$ and $\xi_t$ are the short-term and long-term factors defined in \autoref{eq:SS_rn_chi} and \autoref{eq:SS_rn_xi} for risk-neutral processes. Then $x_t$ satisfies the stochastic differential equation \autoref{eq:sde}, with 
$$b(x_t) = \left[ \begin{matrix} -\kappa \chi_t - \lambda_{\chi} \\ \mu_{\xi} - \gamma \xi_t - \lambda_{\xi} \end{matrix} \right], \sigma(x_t) = \left[ \begin{matrix} \sigma_{\chi} & 0 \\ 0 & \sigma_{\xi} \end{matrix} \right], a(x_t) = \sigma(x_t) \sigma(x_t)^\top = \left[ \begin{matrix} \sigma_{\chi}^2 & 0 \\ 0 & \sigma_{\xi}^2 \end{matrix} \right]. $$
The basis $H(x_t) = (1, \chi_t, \xi_t, \chi_t^2, \chi_t \xi_t, \xi_t^2)^\top$ has a dimension $N = 6$. The coordinate representation $\vec{p} = (\alpha_1, \alpha_2, \alpha_3, \alpha_4, \alpha_5, \alpha_6)^\top$. 
The generator $\mathcal{G}$ is given by 
$$\mathcal{G}f(x) = \frac{1}{2} Tr \left( \left[ \begin{matrix} \sigma_{\chi}^2 & 0 \\ 0 & \sigma_{\xi}^2 \end{matrix} \right] \nabla^2 f(x) \right) + \left[ \begin{matrix} -\kappa \chi_t - \lambda_{\chi} \\ \mu_{\xi} - \gamma \xi_t - \lambda_{\xi} \end{matrix} \right]^\top \nabla f(x). $$
By applying $\mathcal{G}$ to each element of $H_n(x_t)$, we obtain the matrix representation 
$$G = \left[ \begin{matrix} 
             0 & -\lambda_{\chi} & \mu_{\xi} - \lambda_{\xi} & \sigma_{\chi}^2 & 0 & \sigma_{\xi}^2 \\
             0 & -\kappa & 0 & -2 \lambda_{\chi} & \mu_{\xi} - \lambda_{\xi} & 0 \\
             0 & 0 & -\gamma & 0 & -\lambda_{\chi} & 2\mu_{\xi} - 2\lambda_{\xi} \\
             0 & 0 & 0 & -2\kappa & 0 & 0 \\ 
             0 & 0 & 0 & 0 & -\kappa - \gamma & 0 \\ 
             0 & 0 & 0 & 0 & 0 & -2\gamma
             \end{matrix} \right].$$ 
Then, by Theorem 1, the futures price $F_{t,T}$ is given by
\begin{equation}
F_{t,T} = \mathbb{E}^*(S_T | \mathcal{F}_t) = H(x_t)^\top e^{(T-t)G} \vec{p}. 
\label{eq:qua_ftt}
\end{equation}
Therefore, we have the non-linear state-space model 
\begin{equation}
x_t = c + E x_{t-1} + w_t, w_t \sim N(\textbf{0}, \Sigma_w), 
\label{eq:qua_xt}
\end{equation}
and 
\begin{equation}
y_t = H(x_t)^\top e^{(T-t)G} \vec{p} + v_t, v_t \sim N(\textbf{0}, \Sigma_v). 
\label{eq:qua_yt}
\end{equation}

## Filtering methods

In this section, we use the notation 
\begin{align}
a_{t|t-1} &:= \mathbb{E}(x_t | \mathcal{F}_{t-1}),& P_{t|t-1} &:= Cov(x_t | \mathcal{F}_{t-1}), \nonumber \\
a_t &:= \mathbb{E}(x_t | \mathcal{F}_t),& P_t &:= Cov(x_t | \mathcal{F}_t). \nonumber
\end{align} 

![Flowcharts of EKF\label{fig:EKF}](paper-figures/EKF.jpg){ width=70% }

![Flowcharts of UKF\label{fig:UKF}](paper-figures/UKF.jpg){ width=60% }

The Kalman Filter (KF) [@harvey:1990] is a commonly used filtering method in estimating hidden state variables. However, KF can only deal with the linear Gaussian state model. To capture the non-linear dynamics in the PD model, we use Extended Kalman Filter (EKF) [@julier:1997] and Unscented Kalman Filter (UKF) [@julier:2004; @wan:2000]. Suppose we have a non-linear state-space model 
$$x_t = f(x_{t-1}) + w_t, w_t \sim N(\textbf{0}, \Sigma_w), $$
$$y_t = h(x_t) + v_t, v_t \sim N(\textbf{0}, \Sigma_v). $$
The EKF linearises the state and measurement equations through the first-order Taylor series. To run KF, we replace $J_f$ and $J_h$ with $E$ and $F_t$ respectively, where $J_f$ and $J_h$ are the Jacobians of $f(\cdot)$ and $h(\cdot)$. In contrast, the UKF uses a set of carefully chosen points, called sigma points, to represent the true distributions of state variables. Then, these sigma points are propagated through the state equation. The flowcharts of EKF and UKF are given in \autoref{fig:EKF} and \autoref{fig:UKF}. In this application, we use KF for the Schwartz-Smith model, and EKF/UKF for the polynomial diffusion model. 

# Statement of need

This application is aimed at researchers who are pricing commodity futures by Schwartz-Smith model or PD model. It has been designed with the following goals: 

\begin{itemize}
\item[1. ] To provide a simulation tool for the polynomial diffusion model. Users can declare all model specifications and parameters. The generated data is downloadable. 
\item[2. ] To provide two filtering methods, EKF and UKF, to estimate the futures prices and hidden state variables. Currently, there is no filtering toolbox for the polynomial diffusion model. 
\item[3. ] To provide well-designed visualisations. That includes the futures prices, the state variables, the estimates of futures prices and state variables, and some downloadable tables. Moreover, all these plots are interactive. Users can zoom in/out, highlight a specific curve, download these plots, etc. 
\item[4. ] To provide the estimation errors including root mean squared error (RMSE), mean absolute error (MSE) and mean relative error (MRE). These measures are presented in tables and plots. 
\item[5. ] To provide all functions listed above for the Schwartz-Smith model as a comparison. 
\end{itemize}

# Comparison with existing libraries

The R package "NFCP" [@aspinall:2022] was developed for multi-factor pricing of commodity futures, which is a generalisation of the Schwartz-Smith model. However, this package doesn't accommodate the polynomial diffusion model. There are no R packages available for PD models currently. 

There are many packages in R for KF, for example, "dse", "FKF", "sspir", "dlm", "KFAS": "dse" can only take time-invariant state and measurement transition matrices; "FKF" emphasizes computation speed but cannot run smoother; "sspir", "dlm" and "KFAS" have no deterministic inputs in state and measurement equations. For the non-linear state-space model, the functions "ukf" and "ekf" in package "bssm" run the EKF and UKF respectively. However, this package was designed for Bayesian inference where a prior distribution of unknown parameters is required. To achieve the best collaboration of filters and models, we developed functions of KF, EKF and UKF within this code. 

# Acknowledgments

We are indebted to Mark Van de Vyver for his constructive suggestions and ongoing support in improving the paper since its submission. Thanks to Mark's suggestions, which included helpful references to the previous work, we were able to achieve code containerisation, and develop and provide the results of the two unit tests.

We would also like to thank Sam Forbes, and Blake Rayfield for their valuable feedback and suggestions, which helped us to improve the PDSim's code.

# References
