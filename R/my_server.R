my_server <- function(input, output, session) {
  input <<- input # make input global
  output <<- output # make output global 
  session <<- session # make session global
  seed <<- 1234 # seed for generating random number 
  
  # Generate new data button
  observeEvent(input$new_data, {
    seed <<- seed + 1
  })
  
  # Download futures data button
  output$download_data<- downloadHandler(filename = "data.csv", content = function(file){
    dat <- sim_dat()
    write.csv(dat$yt, file)
  })
  
  # Download maturity data button
  output$download_mat <- downloadHandler(filename = "maturity.csv", content = function(file){
    dat <- sim_dat()
    write.csv(dat$mats, file)
  })
  
  # Simulate data
  sim_dat <<- reactive({
    a <- input$new_data # This value has no meaning. Just make sure every time the button is clicked, this function will be run again  
    par_org <- c(input$kappa, input$gamma, input$mu, input$sc, input$sx, input$rho, input$lc, input$lx, 
                 seq(from = input$s1, to = input$sn, length.out = input$n_contract)) # all required parameters
    x0 <- c(0, input$mu / input$gamma) # initialisation 
    noise <- "Gaussian"
    dt <- 1/360
    if (input$model == "SS2000") {
      n_coe <- 0 # number of model coefficients
      func_f <- function(xt, par) state_linear(xt, par, dt) # state equation
      func_g <- function(xt, par, mats) measurement_linear(xt, par, mats) # measurement equation
    } else if (input$model == "PD") {
      n_coe <- 6 # number of model coefficients
      par_coe <- c(input$constant, input$chi, input$xi, input$chi2, input$chi_xi, input$xi2) 
      par_org <- c(par_org, par_coe) 
      func_f <- function(xt, par) state_linear(xt, par, dt) # state equation
      func_g <- function(xt, par, mats) measurement_polynomial(xt, par, mats, 2, n_coe) # measurement equation
    }
    dat <- simulate_data(par_org, x0, input$n_obs, input$n_contract, func_f, func_g, n_coe, noise, seed)
    if (input$model == "SS2000") {
      dat$yt <- exp(dat$yt) # true futures prices 
    } 
    dat$yt <- data.frame(dat$yt)
    colnames(dat$yt) <- paste("Contract", 1: input$n_contract, sep = "")
    colnames(dat$mats) <- paste("Contract", 1: input$n_contract, sep = "")
    dat$xt <- data.frame(dat$xt)
    colnames(dat$xt) <- c("chi", "xi")
    dat$yt_long <- pivot_longer(dat$yt, 1: input$n_contract)
    dat$yt_long$observations <- rep(1: input$n_obs, each = input$n_contract)
    dat$xt_long <- pivot_longer(dat$xt, 1: 2)
    dat$xt_long$observations <- rep(1: input$n_obs, each = 2)
    return(dat)
  })
  
  # Estimate model
  estimation <<- reactive({
    dat <- sim_dat()
    yt <- as.matrix(dat$yt)
    if (input$model == "SS2000") {
      yt <- log(yt)
    }
    mats <- as.matrix(dat$mats)
    par_org <- c(input$kappa, input$gamma, input$mu, input$sc, input$sx, input$rho, input$lc, input$lx, 
                 seq(from = input$s1, to = input$sn, length.out = input$n_contract))
    x0 <- c(0, input$mu / input$gamma)
    noise <- "Gaussian"
    dt <- 1/360
    
    if (input$model == "SS2000") {
      func_f <- function(xt, par) state_linear(xt, par, dt)
      func_g <- function(xt, par, mats) measurement_linear(xt, par, mats)
      est <- KF(c(par_org, x0), yt, mats, 0, dt, FALSE, "None")
      est$yt_hat <- data.frame(exp(func_g(t(est$xt_filter), par_org, mats)$y))
    } else if (input$model == "PD") {
      n_coe <- 6
      par_coe <- c(input$constant, input$chi, input$xi, input$chi2, input$chi_xi, input$xi2)
      par_org <- c(par_org, par_coe)
      func_f <- function(xt, par) state_linear(xt, par, dt)
      func_g <- function(xt, par, mats) measurement_polynomial(xt, par, mats, 2, n_coe)
      if (input$filter == "Extended Kalman Filter") {
        est <- EKF(c(par_org, x0), yt, mats, func_f, func_g, dt, n_coe, noise)
      } else if (input$filter == "Unscented Kalman Filter") {
        est <- UKF(c(par_org, x0), yt, mats, func_f, func_g, dt, n_coe, noise)
      }
      est$yt_hat <- data.frame(func_g(t(est$xt_filter), par_org, mats)$y)
    }
    
    colnames(est$yt_hat) <- paste("Contract", 1: input$n_contract, sep = "")
    est$yt_hat_long <- pivot_longer(est$yt_hat, 1: input$n_contract)
    est$yt_hat_long$observations <- rep(1: input$n_obs, each = input$n_contract)
    est$xt_filter <- data.frame(est$xt_filter)
    colnames(est$xt_filter) <- c("chi", "xi")
    est$xt_filter_long <- pivot_longer(est$xt_filter, 1: 2)
    est$xt_filter_long$observations <- rep(1: input$n_obs, each = 2)
    return(est)
  })
  
  # Estimation errors 
  errors <<- reactive({
    yt <- sim_dat()$yt
    yt_hat <- estimation()$yt_hat
    
    err <- list( se = (yt - yt_hat)^2, 
                 ae = abs(yt - yt_hat), 
                 re = abs((yt - yt_hat)/yt) ) # square error, absolute error, relative error
    err$rmse <- sqrt( apply(err$se, 2, mean) )
    err$mae <- apply(err$ae, 2, mean)
    err$mre <- apply(err$re, 2, mean)
    return(err)
  })
  
  # Check if parameters are valid 
  check <<- reactive({
    if (input$n_obs <= 0 || !is.integer(input$n_obs)) {
      showModal(modalDialog("The number of trading days must be a positive integer. "))
    } else if (input$n_contract <= 0 || !is.integer(input$n_contract)) {
      showModal(modalDialog("The number of contracts must be a positive integer. "))
    } else if (input$contract <= 0 || !is.integer(input$contract) || input$contract > input$n_contract) {
      showModal(modalDialog("The contract must be a positive integer and no greater than the number of contracts. "))
    } else if (input$day <= 0 || !is.integer(input$day) || input$day > input$n_obs) {
      showModal(modalDialog("The day must be a positive integer and no greater than the number of trading days. "))
    } else if (input$kappa <= 0) {
      showModal(modalDialog(withMathJax("\\( \\kappa \\) must be positive. ")))
    } else if (input$gamma <= 0) {
      showModal(modalDialog(withMathJax("\\( \\gamma \\) must be positive. ")))
    } else if (input$sc <= 0) {
      showModal(modalDialog(withMathJax("\\( \\sigma_{\\chi} \\) must be positive. ")))
    } else if (input$sx <= 0) {
      showModal(modalDialog(withMathJax("\\( \\sigma_{\\xi} \\) must be positive. ")))
    } else if (input$rho < -1 || input$rho > 1) {
      showModal(modalDialog(withMathJax("\\( \\rho \\) must be between -1 and 1. ")))
    } else if (input$s1 <= 0) {
      showModal(modalDialog("Standard error must be positive. "))
    }else if (input$sn <= 0) {
      showModal(modalDialog("Standard error must be positive. "))
    }
  })
  
  server_welcome() # server for welcome page
  server_app() # server for app page
  server_members() # server for team members page
}