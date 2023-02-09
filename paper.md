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
    orcid: 0000−0002−7853−4550
    equal-contrib: true
    affiliation: 2 
  - name: Gareth W. Peters
    orcid: 0000-0003-2768-8979
    equal-contrib: true
    affiliation: 3
  - name: Pavel V. Shevchenko
    orcid: 0000−0001−8104−8716
    equal-contrib: true
    affiliation: 1
affiliations:
 - name: Department of Actuarial Studies and Business Analytics, Macquarie University, Australia
   index: 1
 - name: Department of Mathematics and Statistics, Macquarie University, Australia
   index: 2
 - name: Department of Statistics and Applied Probability, University of California Santa Barbara, USA
   index: 3
date: 9 February 2023
bibliography: paper.bib

---

# Summary

The Schwartz-Smith two-factor model [@schwartz:2000] was commonly used in the pricing of crude oil futures and some other futures in the last two decades. In 2016, [@filipovic:2016] introduced a new polynomial diffusion framework which allows a more complicated structure of spot price. This framework has been applied to electricity forwards in [@kleisinger:2020], in which the spot price is modelled in a quadratic form of two factors. This application aims to estimate futures prices as well as the latent state variables, and provides well-designed visualisations. The application is available at [https://github.com/peilun-he/polynomial-diffusion-model-simulation-and-estimation](https://github.com/peilun-he/polynomial-diffusion-model-simulation-and-estimation). 

## Schwartz-Smith Two-Factor Model

Under the Schwartz-Smith framework, the spot price $S_t$ is modelled as the sum of two factors $\chi_t$ and $\xi_t$, 
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
where $\kappa, \gamma \in \mathbb{R}^+$ are the speed of mean-reversion parameters, $\mu_{\xi} \in \mathbb{R}$ is the mean level of the long-term factor, $\sigma_{\chi}, \sigma_{\xi} \in \mathbb{R}^+$ are the volatilities, and $\lambda_{\chi}, \lambda_{\xi} \in \mathbb{R}$ are risk premiums. The processes $(W_t^{\chi})_{t \ge 0}$ and $(W_t^{\xi})_{t \ge 0}$ are correlated standard Brownian Motions with 
$$\mathbb{E} \left(dW_t^{\chi} dW_t^{\xi}\right) = \rho dt. $$ 
We set $\lambda_{\chi} = \lambda_{\xi} = 0$ in \eqref{eq:SS_rn_chi} and \eqref{eq:SS_rn_xi} to get the real processes. We use the risk-neutral processes for futures pricing, and real processes for modelling of state variables. 

In discrete time, $\chi_t$ and $\xi_t$ are jointly normally distributed. Therefore, the spot price, which is equal to the sum of $\chi_t$ and $\xi_t$, is log-normally distributed. Moreover, under the arbitrage-free assumption, the futures price ($F_{t,T}$) at current time $t$ must be equal to the expected value of spot price at maturity time $T$, 
$$\log{(F_{t,T})} = \log{\left(\mathbb{E}^*(S_T | \mathcal{F}_t)\right)}, $$
where $\mathcal{F}_t$ is the filtration and $\mathbb{E}^*(\cdot)$ is the expectation under the risk-neutral processes \eqref{eq:SS_rn_chi} and \eqref{eq:SS_rn_xi}. Then we can get the linear Gaussian state space model: 
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

## Polynomial Diffusion Model

In this section, we present a general framework of the polynomial diffusion model first, then we give the application in the two-factor model. The mathematical foundations and applications of polynomial diffusion model in finance are provided in [@filipovic:2016].

Consider the stochastic differential equation
\begin{equation}
dX_t = b(X_t)dt + \sigma(X_t)dW_t, 
\label{eq:sde}
\end{equation}
where $W_t$ is a $d$-dimensional standard Brownian motion and map $\sigma: \mathbb{R}^d \to \mathbb{R}^{d \times d}$ is continuous. Define $a := \sigma \sigma^\top$. For maps $a: \mathbb{R}^d \to \mathbb{S}^{d}$ and $b: \mathbb{R}^d \to \mathbb{R}^d$, suppose we have 
$a_{ij} \in Pol_2$ and $b_i \in Pol_1$. $\mathbb{S}^d$ is the set of all real symmetric $d \times d$ matrices and $Pol_n$ is the set of all polynomials of degree at most $n$. Then the solution of \eqref{eq:sde} is a polynomial diffusion. 

Moreover, we define the generator $\mathcal{G}$ associated to the polynomial diffusion $X_t$ as
\begin{equation}
\mathcal{G}f(x) = \frac{1}{2} Tr\left( a(x) \nabla^2 f(x)\right) + b(x)^\top \nabla f(x)
\label{eq:generator}
\end{equation}
for $x \in \mathbb{R}^d$ and any $C^2$ function $f$. Let $N$ be the dimension of $Pol_n$, and $H: \mathbb{R}^d \to \mathbb{R}^N$ be a function whose components form a basis of $Pol_n$. Then for any $p \in Pol_n$, there exists a unique vector $\vec{p} \in \mathbb{R}^N$ such that 
\begin{equation}
p(x) = H(x)^\top \vec{p}
\label{eq:vec_p}
\end{equation}
and $\vec{p}$ is the coordinate representation of $p(x)$. Moreover, there exists a unique matrix representation $G \in \mathbb{R}^{N \times N}$ of the generator $\mathcal{G}$, such that $G \vec{p}$ is the coordinate vector of $\mathcal{G} p$. So we have
\begin{equation}
\mathcal{G} p(x) = H(x)^\top G \vec{p}. 
\label{eq:G}
\end{equation}

Theorem 1: Let $p(x) \in Pol_n$ be a polynomial with coordinate representation $\vec{p} \in \mathbb{R}^N$ satisfying \eqref{eq:vec_p}, $G \in \mathbb{R}^{N \times N}$ be a matrix representation of generator $\mathcal{G}$ satisfying \eqref{eq:G}, and $X_t \in \mathbb{R}^d$ satisfies \eqref{eq:sde}. Then for $0 \le t \le T$, we have
	$$\mathbb{E} \left(p(X_T)|\mathcal{F}_t \right) = H(X_t)^\top e^{(T-t)G} \vec{p}, $$
	where $\mathcal{F}_t$ is a $\sigma$-algebra generated up tp time $t$. 
	\label{th:pd}

The proof of Theorem 1 is given in [@filipovic:2016]. 

Next, we apply this theorem to the two-factor model. Assume the spot price $S_t$ is modelled as
\begin{equation}
S_t = p_n(x_t), 
\label{eq:poly_st}
\end{equation}
where $x_t = (\chi_t, \xi_t)^\top$ is a vector of state variables and $p_n(\cdot)$ is a polynomial function with a degree at most $n$. $\chi_t$ and $\xi_t$ are the short-term and long-term factors defined in \eqref{eq:SS_rn_chi} and \eqref{eq:SS_rn_xi} for risk-neutral processes. $x_t$ satisfies the stochastic differential equation \eqref{eq:sde}, with 
$$b(x_t) = \left[ \begin{matrix} -\kappa \chi_t - \lambda_{\chi} \\ \mu_{\xi} - \gamma \xi_t - \lambda_{\xi} \end{matrix} \right], \sigma(x_t) = \left[ \begin{matrix} \sigma_{\chi} & 0 \\ 0 & \sigma_{\xi} \end{matrix} \right], a(x_t) = \sigma(x_t) \sigma(x_t)^\top = \left[ \begin{matrix} \sigma_{\chi}^2 & 0 \\ 0 & \sigma_{\xi}^2 \end{matrix} \right]. $$
For any basis $H_n(x_t)$, the polynomial $p_n(x_t)$ can be uniquely represented as
$$p_n(x_t) = H_n(x_t)^\top \vec{p}. $$
The generator $\mathcal{G}$ is given by 
$$\mathcal{G}f(x) = \frac{1}{2} Tr \left( \left[ \begin{matrix} \sigma_{\chi}^2 & 0 \\ 0 & \sigma_{\xi}^2 \end{matrix} \right] \nabla^2 f(x) \right) + \left[ \begin{matrix} -\kappa \chi_t - \lambda_{\chi} \\ \mu_{\xi} - \gamma \xi_t - \lambda_{\xi} \end{matrix} \right]^\top \nabla f(x). $$
By applying $\mathcal{G}$ to each element of $H_n(x_t)$, we get the matrix representation $G$. Then, by Theorem \ref{th:pd}, the futures price $F_{t,T}$ is given by
\begin{equation}
F_{t,T} = \mathbb{E}^*(S_T | \mathcal{F}_t) = H(x_t)^\top e^{(T-t)G} \vec{p}. 
\label{eq:qua_ftt}
\end{equation}
Therefore, we have the non-linear state-space model 
\begin{equation}
x_t = c + E x_{t-1} + w_t, w_t \sim N(\textbf{0}, W), 
\label{eq:qua_xt}
\end{equation}
and 
\begin{equation}
y_t = H_n(x_t)^\top e^{(T-t)G} \vec{p} + v_t, v_t \sim N(\textbf{0}, V). 
\label{eq:qua_yt}
\end{equation}

In this application, we assume the spot price is a polynomial of state variables with degree $n = 2$, 
$S_t = \alpha_1 + \alpha_2 \chi_t + \alpha_3 \xi_t + \alpha_4 \chi_t^2 + \alpha_5 \chi_t \xi_t + \alpha_6 \xi_t^2$, and the dimension of $Pol_2$ is $N = 6$. The coordinate representation $\vec{p} = (\alpha_1, \alpha_2, \alpha_3, \alpha_4, \alpha_5, \alpha_6)^\top$. $\alpha_i, i = 1, \dots, 6$ are parameters which can be chosen by users. 

# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
