draw_user_guide <- function() {
  style = "font-size:25px"
  tabPanel("User Guide", 
           withMathJax(), 
           titlePanel(h1("User Guide")), 
           titlePanel(h3("Hints to use PDSim")), 
           p("1. Once users enter all parameters, the data will be generated automatically. 
             Users do NOT need to click any buttons. However, if users wish to generate more
             realisations under the same set of parameters, please click the 'Generate new 
             data' button. ", style = style), 
           p("2. The seed to generate random numbers is fixed, i.e., for the same set of 
             parameters, users will get exactly the same data every time they use PDSim. ", 
             style = style), 
           p("3. Futures prices in all tables / plots are REAL prices (NOT the logarithm), 
             no matter which model is used. ", style = style), 
           p("4. The 95% confidence interval is shown as a grey ribbon on each plot. ", style = style), 
           p("5. Because of the limitation of filtering methods, the standard error of 
             the estimated futures price on the first day is extremely large. All plots
             of contracts estimation start from the second day. ", style = style),
           p(paste("6. This is Version", ver, "of PDSim. The newest version is available at: "), 
             a("https://github.com/peilun-he/PDSim", 
               href = "https://github.com/peilun-he/PDSim", 
               target = "_blank"), style = style), 
           hr(),
           titlePanel(h3("Schwartz-Smith model")), 
           p("The spot price \\( S_t \\) is modelled as", style = style), 
           p("\\( \\log(S_t) = \\chi_t + \\xi_t, \\)", class = "text-center", style = style), 
           p("where \\( \\chi_t \\) represents the short-term fluctuation and \\( \\xi_t \\) is 
             the long-term equilibrium price level. We assume both \\( \\chi_t \\) and \\( \\xi_t \\) 
             follow an Ornstein–Uhlenbeck process", style = style), 
           p("\\( d\\chi_t = - \\kappa \\chi_t dt + \\sigma_{\\chi} dW_t^{\\chi} \\)", 
             class = "text-center", style = style),
           p("\\( d\\xi_t = (\\mu_{\\xi} - \\gamma \\xi_t) dt + \\sigma_{\\xi} dW_t^{\\xi} \\)", 
             class = "text-center", style = style),
           p("for real processes and", style = style), 
           p("\\( d\\chi_t = (- \\kappa \\chi_t - \\lambda_{\\chi}) dt + \\sigma_{\\chi} dW_t^{\\chi*} \\)", 
             class = "text-center", style = style),
           p("\\( d\\xi_t = (\\mu_{\\xi} - \\gamma \\xi_t - \\lambda_{\\xi}) dt + \\sigma_{\\xi} dW_t^{\\xi*} \\)", 
             class = "text-center", style = style),
           p("for risk-neutral processes. \\( \\kappa, \\gamma \\in \\mathbb{R}^+ \\) are speed of 
             mean-reversion parameters. \\( \\mu_{\\xi} \\in \\mathbb{R} \\) is the mean level 
             of the long-term factor. \\( \\sigma_{\\chi}, \\sigma_{\\xi} \\in \\mathbb{R}^+ \\) 
             are volatilities. \\( \\lambda_{\\chi}, \\lambda_{\\xi} \\in \\mathbb{R} \\) are risk 
             premiums. \\( W_t^{\\chi*} \\) and \\( W_t^{\\xi*} \\) are 
             correlated standard Brownian Motions with correlation coefficient \\( \\rho \\). 
             In the original Schwartz-Smith model, \\( \\gamma = 0 \\), but here we extend this model 
             to allow non-zero speed of mean-reversion parameter for the long-term factor. 
             Under the arbitrage-free assumption, the futures price \\( T_{t,T} \\) at current 
             time \\( t \\) with maturity time \\( T \\) must be equal to the expected value of 
             spot price at maturity time, i.e., ", 
             style = style), 
           p("\\( \\log(F_{t,T}) = \\log(\\mathbb{E}^*(S_T | \\mathcal{F}_t)), \\)", class = "text-center", style = style), 
           p("where \\( F_t \\) is the filtration until time t. In discrete time, we have the 
           linear Gaussian state-space model: ", style = style), 
           p("\\( x_t = c + E x_{t-1} + w_t, \\)", class = "text-center", style = style), 
           p("\\( y_t = d_t + F_t x_t + v_t, \\)", class = "text-center", style = style), 
           p("where \\( x_t = \\left[ \\begin{matrix} \\chi_t \\\\ \\xi_t \\end{matrix} \\right] \\), 
             \\( c = \\left[ \\begin{matrix} 0 \\\\ \\frac{\\mu_{\\xi}}{\\gamma} (1 - e^{-\\gamma \\Delta t}) \\end{matrix} \\right] \\), 
             \\( E = \\left[ \\begin{matrix} e^{-\\kappa \\Delta t} & 0 \\\\ 0 & e^{-\\gamma \\Delta t} \\end{matrix} \\right] \\), 
             \\( y_t = \\left( \\log{(F_{t,T_1})}, \\dots, \\log{(F_{t,T_m})} \\right)^\\top \\), 
             \\( d_t = \\left( A(T_1 - t), \\dots, A(T_m - t) \\right)^\\top \\), 
             \\( F_t = \\left[ \\begin{matrix} e^{-\\kappa (T_1 - t)}, \\dots, e^{-\\kappa (T_m - t)} \\\\ 
             e^{-\\gamma (T_1 - t)}, \\dots, e^{-\\gamma (T_m - t)} \\end{matrix} \\right]^\\top \\)
             and \\( m \\) is the number of contracts. The function \\( A(\\cdot) \\) is given by", 
             style = style), 
           p("\\( A(t) = -\\frac{\\lambda_{\\chi}}{\\kappa}(1 - e^{-\\kappa t}) + 
             \\frac{\\mu_{\\xi} - \\lambda_{\\xi}}{\\gamma}(1 - e^{-\\gamma t}) + 
             \\frac{1}{2}\\left(\\frac{1 - e^{-2\\kappa t}}{2\\kappa}\\sigma_{\\chi}^2 + 
             \\frac{1 - e^{-2\\gamma t}}{2\\gamma}\\sigma_{\\xi}^2 + 
             2\\frac{1 - e^{-(\\kappa + \\gamma)t}}{\\kappa + \\gamma}\\sigma_{\\chi}\\sigma_{\\xi}\\rho \\right). \\)", 
             class = "text-center", style = style), 
           p("\\( w_t \\) and \\( v_t \\) are multivariate Gaussian noises with mean 0 
             and covariance matrix \\( \\Sigma_w \\) and \\( \\Sigma_v \\) respectively, 
             where", style = style), 
           p("\\( \\Sigma_w = \\left[ \\begin{matrix} 
             \\frac{1 - e^{-2\\kappa \\Delta t}}{2\\kappa}\\sigma_{\\chi}^2 & 
             \\frac{1 - e^{-(\\kappa + \\gamma) \\Delta t}}{\\kappa + \\gamma} \\sigma_{\\chi}\\sigma_{\\xi}\\rho \\\\ 
             \\frac{1 - e^{-(\\kappa + \\gamma) \\Delta t}}{\\kappa + \\gamma} \\sigma_{\\chi}\\sigma_{\\xi}\\rho & 
             \\frac{1 - e^{-2\\gamma \\Delta t}}{2\\gamma}\\sigma_{\\xi}^2 \\end{matrix} \\right], \\)", 
             class = "text-center", style = style),
           p("and we assume \\( \\Sigma_v \\) is diagonal, i.e., ", style = style),
           p("\\( \\Sigma_v = \\left[ \\begin{matrix} 
             \\sigma_1^2 & 0 & \\dots & 0 \\\\ 
             0 & \\sigma_2^2 & \\dots & 0 \\\\ 
             \\vdots & \\vdots & \\ddots & \\vdots \\\\ 
             0 & 0 & \\dots & \\sigma_m^2 
             \\end{matrix} \\right]. \\)", class = "text-center", style = style), 
           p("Moreover, we assume \\( \\sigma_1, \\dots \\sigma_m \\) are evenly spaced, i.e., 
             \\( \\sigma_1 - \\sigma_2 = \\sigma_2 - \\sigma_3 = \\dots = \\sigma_{m-1} - \\sigma_m \\). ", 
             style = style), 
           hr(), 
           titlePanel(h3("Polynomial diffusion model")), 
           p("Consider the stochastic differential equation", style = style), 
           p("\\( dX_t = b(X_t)dt + \\sigma(X_t)dW_t, \\)", class = "text-center", style = style),
           p("where \\( W_t \\) is a \\( d \\)-dimensional standard Brownian Motion and map
             \\( \\sigma: \\mathbb{R}^d \\to \\mathbb{R}^{d \\times d} \\) is continuous. 
             define \\( a := \\sigma \\sigma^\\top \\). For maps \\( a: \\mathbb{R}^d \\to 
             \\mathbb{S}^{d} \\) and \\( b: \\mathbb{R}^d \\to \\mathbb{R}^d \\), suppose 
             we have \\( a_{ij} \\in Pol_2 \\) and \\( b_i \\in Pol_1 \\). \\( \\mathbb{S}^d \\) 
             is the set of all real symmetric \\( d \\times d \\) matrices and \\( Pol_n \\) is 
             the set of all polynomials of degree at most \\(n\\). Then the solution of the SDE is 
             a polynomial diffusion. Moreover, we define the generator \\( \\mathcal{G} \\) associated 
             to the polynomial diffusion \\( X_t \\) as ", style = style),
           p("\\( \\mathcal{G}f(x) = \\frac{1}{2} Tr\\left( a(x) \\nabla^2 f(x)\\right) + b(x)^\\top \\nabla f(x) \\)", 
             class = "text-center", style = style),
           p("for \\( x \\in \\mathbb{R}^d \\) and any \\( C^2 \\) function \\(f\\). Let 
             \\(N\\) be the dimension of \\(Pol_n\\), and \\(H: \\mathbb{R}^d \\to \\mathbb{R}^N\\) 
             be a function whose components form a basis of \\(Pol_n\\). Then for any \\( p \\in 
             Pol_n\\), there exists a unique vector \\( \\vec{p} \\in \\mathbb{R}^N\\) such that", 
             style = style),
           p("\\( p(x) = H(x)^\\top \\vec{p} \\)", class = "text-center", style = style),
           p("and \\( \\vec{p} \\) is the coordinate representation of \\( p(x) \\). 
             Moreover, there exists a unique matrix representation \\( G \\in \\mathbb{R}^{N \\times N} \\) 
             of the generator \\( \\mathcal{G} \\), such that \\( G \\vec{p} \\) is the 
             coordinate vector of \\( \\mathcal{G} p \\). So we have", style = style), 
           p("\\( \\mathcal{G} p(x) = H(x)^\\top G \\vec{p}. \\)", class = "text-center", style = style), 
           br(),
           p("Theorem 1: Let \\( p(x) \\in Pol_n \\) be a polynomial with coordinate 
             representation \\( \\vec{p} \\in \\mathbb{R}^N \\), \\( G \\in \\mathbb{R}^{N \\times N} \\) 
             be a matrix representation of generator \\( \\mathcal{G} \\), and \\( X_t \\in \\mathbb{R}^d \\)
             satisfies the SDE. Then for \\( 0 \\le t \\le T \\), we have", style = style),
           p("\\( \\mathbb{E} \\left[p(X_T)|\\mathcal{F}_t \\right] = H(X_t)^\\top e^{(T-t)G} \\vec{p}, \\)", 
             class = "text-center", style = style),
           p("where \\( \\mathcal{F}_t \\) represents all information available until time \\(t\\). ", style = style), 
           br(), 
           br(), 
           br(), 
           p("Next, we apply this theorem to the two-factor model. Assume the spot price \\(S_t\\) is modelled as", style = style), 
           p("\\( S_t = p_n(x_t) = \\alpha_1 + \\alpha_2 \\chi_t + \\alpha_3 \\xi_t + 
             \\alpha_4 \\chi_t^2 + \\alpha_5 \\chi_t \\xi_t + \\alpha_6 \\xi_t^2\\), ", class = "text-center", style = style), 
           p("which is a polynomial with order \\( n = 2 \\).  \\( x_t = (\\chi_t, \\xi_t)^\\top \\)
             is a vector of state variables. Obviously, \\( x_t \\) satisfies the SDE with", style = style), 
           p("\\( b(x_t) = \\left[ \\begin{matrix} -\\kappa \\chi_t - \\lambda_{\\chi} \\\\ 
             \\mu_{\\xi} - \\gamma \\xi_t - \\lambda_{\\xi} \\end{matrix} \\right], 
             \\sigma(x_t) = \\left[ \\begin{matrix} \\sigma_{\\chi} & 0 \\\\ 
             0 & \\sigma_{\\xi} \\end{matrix} \\right], 
             a(x_t) = \\sigma(x_t) \\sigma(x_t)^\\top = \\left[ \\begin{matrix} \\sigma_{\\chi}^2 & 0 \\\\ 
             0 & \\sigma_{\\xi}^2 \\end{matrix} \\right]. \\)", class = "text-center", style = style), 
           p("The basis", style = style), 
           p("\\( H(x_t) = (1, \\chi_t, \\xi_t, \\chi_t^2, \\chi_t \\xi_t, \\xi_t^2)^\\top \\)", class = "text-center", style = style), 
           p("has a dimension \\( N = 6 \\). The coordinate representation is ", style = style), 
           p("\\( \\vec{p} = (\\alpha_1, \\alpha_2, \\alpha_3, \\alpha_4, \\alpha_5, \\alpha_6)^\\top. \\)", class = "text-center", style = style), 
           p("By applying \\( \\mathcal{G} \\) to each element of \\( H(x_t) \\), 
             we get the matrix representation", style = style), 
           p("\\( G = \\left[ \\begin{matrix} 
             0 & -\\lambda_{\\chi} & \\mu_{\\xi} - \\lambda_{\\xi} & \\sigma_{\\chi}^2 & 0 & \\sigma_{\\xi}^2 \\\\
             0 & -\\kappa & 0 & -2 \\lambda_{\\chi} & \\mu_{\\xi} - \\lambda_{\\xi} & 0 \\\\
             0 & 0 & -\\gamma & 0 & -\\lambda_{\\chi} & 2\\mu_{\\xi} - 2\\lambda_{\\xi} \\\\
             0 & 0 & 0 & -2\\kappa & 0 & 0 \\\\ 
             0 & 0 & 0 & 0 & -\\kappa - \\gamma & 0 \\\\ 
             0 & 0 & 0 & 0 & 0 & -2\\gamma
             \\end{matrix} \\right]. \\)", class = "text-center", style = style),
           p("Then, by Theorem 1, the futures price \\(F_{t,T}\\) is given by", style = style), 
           p("\\( F_{t,T} = \\mathbb{E}^*(S_T | \\mathcal{F}_t) = H(x_t)^\\top e^{(T-t)G} \\vec{p}. \\)", class = "text-center", style = style),
           p("Therefore, we have the non-linear state-space model ", style = style), 
           p("\\( x_t = c + E x_{t-1} + w_t, \\)", class = "text-center", style = style), 
           p("\\( y_t = H(x_t)^\\top e^{(T-t)G} \\vec{p} + v_t. \\)", class = "text-center", style = style), 
           br()
           )
}
