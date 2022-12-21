library(shiny)
library(ggplot2)
library(DT)
library(plotly)
library(shinythemes)
library(tidyr)
library(scales)

source("./Functions/decomposition_eigen.R")
source("./Functions/G_matrix.R")
source("./Functions/AofT.R")
source("./Functions/measurement_polynomial.R")
source("./Functions/measurement_linear.R")
source("./Functions/simulate_data.R")
source("./Functions/state_linear.R")
source("./Functions/taylor_coe.R")
source("./Functions/EKF.R")
source("./Functions/UKF.R")
source("./Functions/KF.R")

# Define UI for application
ui <- fluidPage(
  titlePanel(h1("State Variables and contracts estimations for the polynomial diffusion model")), 
  theme = shinytheme("cyborg"), 
  #themeSelector(), 
  sidebarLayout(
    sidebarPanel(
      titlePanel(h3("Model Specification")), 
      fluidRow(
        column(5, 
          numericInput("n_obs", "No. of trading dates", 100, step = 100)
        ), 
        column(5, offset = 2,
          numericInput("n_contract", "No. of contracts", 10, step = 1)
        )
      ), 
      numericInput("seed", "Enter a random seed", 1234, step = 1), 
      selectInput("model", "Select a model", selected = "PD", 
                  choices = c("Schwartz-Smith two-factor model" = "SS2000", "Polynomial diffusion model" = "PD")),
      #textOutput("selected_model"), 
      conditionalPanel(condition = "input.model == 'PD'", 
                       selectInput("filter", "Select a filtering method", selected = "Extended Kalman Filter", 
                                   choices = c("Extended Kalman Filter", "Unscented Kalman Filter")), 
                       fluidRow(
                         column(5, 
                                numericInput("constant", "Constant", 1), 
                                numericInput("chi", "chi", 1), 
                                numericInput("xi", "xi", 1)
                         ), 
                         column(5, offset = 2, 
                                numericInput("chi2", "chi^2", 1), 
                                numericInput("chi_xi", "chi * xi", 1),
                                numericInput("xi2", "xi^2", 1)
                         )
                       )
      ),
      actionButton("model_instructions", "Model instructions"), 
      hr(),
      
      titlePanel(h3("Estimations")), 
      numericInput("contract", "Enter the number of contracts you want to show estimations", 
                   value = 1, min = 1, step = 1),
      numericInput("day", "Enter a day you want to show estimations", 
                   value = 1, min = 1, step = 1),
      hr(), 
      
      titlePanel(h3("Model parameters")), 
      fluidRow(
        column(5, 
          numericInput("kappa", "kappa", 0.5, step = 0.1), 
          numericInput("mu", "mu_xi", 1, step = 0.1), 
          numericInput("sc", "sigma_chi", 1.5, step = 0.1), 
          numericInput("lc", "lambda_chi", 0.5, step = 0.1), 
          numericInput("s1", "Measurement error of the first available contract", 0.1, step = 0.01)
        ), 
        column(5, offset = 2, 
          numericInput("gamma", "gamma", 0.3, step = 0.1), 
          numericInput("rho", "rho", -0.3, step = 0.1), 
          numericInput("sx", "sigma_xi", 1.3, step = 0.1), 
          numericInput("lx", "lambda_xi", 0.3, step = 0.1), 
          numericInput("sn", "Measurement error of the last available contract ", 0.01, step = 0.01)
        )
      ), 
      #actionButton("run", "Run"), 
      actionButton("help", "Help")
    ), 
    mainPanel(
      tabsetPanel(
        tabPanel("Simulated futures contracts", plotlyOutput("plot_yt")), 
        tabPanel("Simulated state variables", plotlyOutput("plot_xt")), 
        tabPanel("Simulated futures", DTOutput("table_yt")), 
        tabPanel("Simulated maturities (in year)", DTOutput("table_mat"))
      ), 
      tabsetPanel(
        tabPanel("Simulated chi vs estimated chi", plotlyOutput("plot_chi")), 
        tabPanel("Simulated xi vs estimated xi", plotlyOutput("plot_xi")), 
        tabPanel("Simulated vs estimated contract", plotlyOutput("plot_yt_hat")), 
        tabPanel("Prices vs maturities on the given day", plotlyOutput("plot_price_mat"))
      ), 
      tabsetPanel(
        tabPanel("Estimation errors for each contract", DTOutput("table_errors")), 
        tabPanel("Mean relative errors on each trading date", plotlyOutput("plot_re")),
        tabPanel("Mean absolute errors on each trading date", plotlyOutput("plot_ae")), 
        tabPanel("Distribution of relative errors", plotlyOutput("plot_boxplot_re")), 
        tabPanel("Distribution of absolute errors", plotlyOutput("plot_boxplot_ae")),
      )
    )
  )
)

# Define server logic required 
server <- function(input, output) {
  # Model instructions button
  observeEvent(input$model_instructions, {
    showModal(modalDialog("Model instructions. "))
  })
  
  # Help button
  observeEvent(input$help, {
    showModal(modalDialog("Introduction of PD model and explanation of some notations. "))
  })
  
  # Simulate data
  #sim_dat <- eventReactive(input$run, {
  sim_dat <- reactive({
    par_org <- c(input$kappa, input$gamma, input$mu, input$sc, input$sx, input$rho, input$lc, input$lx, 
                 seq(from = input$s1, to = input$sn, length.out = input$n_contract))
    x0 <- c(0, input$mu / input$gamma)
    noise <- "Gaussian"
    dt <- 1/360
    if (input$model == "SS2000") {
      n_coe <- 0
      func_f <- function(xt, par) state_linear(xt, par, dt)
      func_g <- function(xt, par, mats) measurement_linear(xt, par, mats)
    } else if (input$model == "PD") {
      n_coe <- 6
      par_coe <- c(input$constant, input$chi, input$xi, input$chi2, input$chi_xi, input$xi2)
      par_org <- c(par_org, par_coe)
      func_f <- function(xt, par) state_linear(xt, par, dt)
      func_g <- function(xt, par, mats) measurement_polynomial(xt, par, mats, 2, n_coe)
    }
    dat <- simulate_data(par_org, x0, input$n_obs, input$n_contract, func_f, func_g, n_coe, noise, input$seed)
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
  estimation <- reactive({
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
  errors <- reactive({
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
  
  # Text: selected_model
  #output$selected_model <- renderText({
  #  if (input$model == "SS2000") {
  #    msg <- "You select the Schwartz-Smith two-factor model. The state variables are estimated through standard Kalman Filter. "
  #  } else if (input$model == "PD") {
  #    msg <- "You select the polynomial diffusion model. Please choose a filtering method and enter the values of coefficients. "
  #  }
  #})
  
  # Plot: simulated futures contracts
  output$plot_yt <- renderPlotly({
    dat <- sim_dat()
    ggplot(dat$yt_long, aes(x = observations, y = value, col = name)) + 
      geom_line() + 
      xlab("Dates") + 
      ylab("Futures prices")
  })
  
  # Plot: simulated state variables
  output$plot_xt <- renderPlotly({
    dat <- sim_dat()
    ggplot(dat$xt_long, aes(x = observations, y = value, col = name)) + 
      geom_line() + 
      xlab("Dates") + 
      ylab("Values")
  })
  
  # Table: simulated futures contracts
  output$table_yt <- renderDT({
    datatable(
      round(sim_dat()$yt, 4), 
      options = list(
        columnDefs = list( list(className = 'dt-center', targets = "_all") )
      ) )
  })
  
  # Table: simulated maturities
  output$table_mat <- renderDT({
    datatable(
      round(sim_dat()$mats, 4), 
      options = list(
        columnDefs = list( list(className = 'dt-center', targets = "_all") )
      ) )
  })
  
  # Plot: simulated chi vs estimated chi
  output$plot_chi <- renderPlotly({
    xt <- sim_dat()$xt
    xt_filter <- estimation()$xt_filter
    dat <- data.frame(observations = c( 1: input$n_obs, 1: input$n_obs ), 
                      chi = c(xt$chi, xt_filter$chi), 
                      source = c(rep("Simulated", input$n_obs), rep("Estimated", input$n_obs)))

    ggplot(dat, aes(x = observations, y = chi, col = source)) + 
      geom_line() + 
      scale_color_manual(values = c("red", "black")) + 
      xlab("Dates") + 
      ylab("Values")
  })
  
  # Plot: simulated xi vs estimated xi
  output$plot_xi <- renderPlotly({
    xt <- sim_dat()$xt
    xt_filter <- estimation()$xt_filter
    dat <- data.frame(observations = c( 1: input$n_obs, 1: input$n_obs ), 
                      xi = c(xt$xi, xt_filter$xi), 
                      source = c(rep("Simulated", input$n_obs), rep("Estimated", input$n_obs)))
    
    ggplot(dat, aes(x = observations, y = xi, col = source)) + 
      geom_line() + 
      scale_color_manual(values = c("red", "black")) + 
      xlab("Dates") + 
      ylab("Values")
  })
  
  # Plot: first available contract
  output$plot_yt_hat <- renderPlotly({
    yt <- sim_dat()$yt
    yt_hat <- estimation()$yt_hat
    dat <- data.frame(observations = c( 1: input$n_obs, 1: input$n_obs ), 
                      y = c(yt[, input$contract], yt_hat[, input$contract]), 
                      source = c(rep("Simulated", input$n_obs), rep("Estimated", input$n_obs)))
    
    ggplot(dat, aes(x = observations, y = y, col = source)) + 
      geom_line() + 
      scale_color_manual(values = c("red", "black")) + 
      xlab("Dates") + 
      ylab("Futures prices")
  })
  
  # Plot: prices on a given day
  output$plot_price_mat <- renderPlotly({
    yt <- sim_dat()$yt
    yt_hat <- estimation()$yt_hat
    mats <- sim_dat()$mats
    dat <- data.frame(y = c(unlist(yt[input$day, ]), unlist(yt_hat[input$day, ])), 
                      mats = c(mats[input$day, ], mats[input$day, ]), 
                      source = c(rep("Simulated", input$n_contract), rep("Estimated", input$n_contract)))
    print(dat)
    ggplot(dat, aes(x = mats, y = y, col = source)) + 
      geom_line() + 
      scale_color_manual(values = c("red", "black")) + 
      xlab("Time to maturities (in year)") + 
      ylab("Futures prices")
  })

  # Table: errors 
  output$table_errors <- renderDT({
    err <- errors()
    datatable(
      data.frame( `Root mean square error` = round(err$rmse, 4), 
                  `Mean absolute error` = round(err$mae, 4), 
                  `Mean relative error` = label_percent()(round(err$mre, 4)) ), 
      options = list(
        columnDefs = list( list(className = 'dt-center', targets = "_all") )
      ) )
  })
  
  # Plot: mean relative errors on each trading date
  output$plot_re <- renderPlotly({
    err <- errors()
    dat <- data.frame(observations = 1: input$n_obs, re = apply(err$re, 1, mean))
    ggplot(dat, aes(x = observations, y = re)) + 
      geom_line() + 
      xlab("Dates") + 
      ylab("Relative errors")
  })
  
  # Plot: mean absolute errors on each trading date
  output$plot_ae <- renderPlotly({
    err <- errors()
    dat <- data.frame(observations = 1: input$n_obs, ae = apply(err$ae, 1, mean))
    ggplot(dat, aes(x = observations, y = ae)) + 
      geom_line() + 
      xlab("Dates") + 
      ylab("Absolute errors")
  })
  
  # Plot: boxplot of relative errors
  output$plot_boxplot_re <- renderPlotly({
    err <- errors()
    dat <- pivot_longer(err$re, 1: input$n_contract)
    ggplot(dat, aes(x = name, y = value, col = name)) + 
      geom_boxplot() + 
      xlab("Contracts") + 
      ylab("Relative errors")
  })
  
  # Plot: boxplot of absolute errors
  output$plot_boxplot_ae <- renderPlotly({
    err <- errors()
    dat <- pivot_longer(err$ae, 1: input$n_contract)
    ggplot(dat, aes(x = name, y = value, col = name)) + 
      geom_boxplot() + 
      xlab("Contracts") + 
      ylab("Absolute errors")
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
