# Server logic for the app page
server_app <- function() {
  observe(updateNumericInput(session, "contract", max = input$n_contract)) # change the maximum value of input 
  observe(updateNumericInput(session, "day", max = input$n_obs)) # change the maximum value of input 
  
  # Plot: simulated futures contracts
  output$plot_yt <- renderPlotly({
    check() # check if all inputs are valid
    dat <- sim_dat()
    date <- dat$yt_long$observations
    price <- dat$yt_long$value
    contract <- factor(dat$yt_long$name, levels = paste("Contract", 1: input$n_contract, sep = ""))
    
    ggplot(mapping = aes(x = date, y = price, colour = contract)) + 
      geom_line() + 
      xlab("Dates") + 
      ylab("Futures prices") + 
      labs(col = "Contracts")
  })
  
  # Plot: 3D surface of futures 
  output$plot_3d_yt <- renderPlotly({
    dat <- sim_dat()
    fig <- plot_ly(x = 1: input$n_contract, y = 1: input$n_obs, z = as.matrix(dat$yt))
    fig <- add_surface(fig)
    fig <- layout(fig, scene = list(xaxis = list(title = "Maturities"), yaxis = list(title = "Dates"), zaxis = list(title = "Futures prices")))
    fig
  })
  
  # Plot: simulated state variables
  output$plot_xt <- renderPlotly({
    dat <- sim_dat()
    date <- dat$xt_long$observations
    value <- dat$xt_long$value
    state <- dat$xt_long$name
    
    ggplot(mapping = aes(x = date, y = value, col = state)) + 
      geom_line() + 
      scale_color_manual(name = "State variables", values = c("red", "black")) +
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
    Pt_filter <- estimation()$Pt_filter
    CI_lower <- xt_filter$chi - 1.96 * sqrt(Pt_filter[1,1,])
    CI_upper <- xt_filter$chi + 1.96 * sqrt(Pt_filter[1,1,])
    date <- 1: input$n_obs
    sim <- xt$chi
    est <- xt_filter$chi
    colors <- c("Simulated" = "black", "Estimated" = "red")
    
    ggplot(mapping = aes(x = date)) + 
      geom_line(aes(y = sim, color = "Simulated")) + 
      geom_line(aes(y = est, color = "Estimated")) +
      geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) + 
      labs(x = "Dates", y = "Values", color = "") +
      scale_color_manual(values = colors) 
  })
  
  # Plot: simulated xi vs estimated xi
  output$plot_xi <- renderPlotly({
    xt <- sim_dat()$xt
    xt_filter <- estimation()$xt_filter
    Pt_filter <- estimation()$Pt_filter
    CI_lower <- xt_filter$xi - 1.96 * sqrt(Pt_filter[2,2,])
    CI_upper <- xt_filter$xi + 1.96 * sqrt(Pt_filter[2,2,])
    date <- 1: input$n_obs
    sim <- xt$xi
    est <- xt_filter$xi
    colors <- c("Simulated" = "black", "Estimated" = "red")
    
    ggplot(mapping = aes(x = date)) + 
      geom_line(aes(y = sim, color = "Simulated")) + 
      geom_line(aes(y = est, color = "Estimated")) +
      geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) + 
      labs(x = "Dates", y = "Values", color = "") +
      scale_color_manual(values = colors) 
  })
  
  # Plot: estimation of one contract
  output$plot_yt_hat <- renderPlotly({
    yt <- sim_dat()$yt
    yt_hat <- estimation()$yt_hat
    cov_y <- estimation()$cov_y
    
    if (input$model == "SS2000") {
      CI_lower <- qlnorm(0.025, meanlog = log(yt_hat[, input$contract]), sdlog = sqrt(cov_y[input$contract,input$contract,]))
      CI_upper <- qlnorm(0.975, meanlog = log(yt_hat[, input$contract]), sdlog = sqrt(cov_y[input$contract,input$contract,]))
    } else if (input$model == "PD") {
      CI_lower <- yt_hat[, input$contract] - 1.96 * sqrt(cov_y[input$contract,input$contract,])
      CI_upper <- yt_hat[, input$contract] + 1.96 * sqrt(cov_y[input$contract,input$contract,])
    }
    
    date <- 2: input$n_obs
    sim <- yt[2: input$n_obs, input$contract]
    est <- yt_hat[2: input$n_obs, input$contract]
    CI_lower <- CI_lower[2: input$n_obs]
    CI_upper <- CI_upper[2: input$n_obs]
    colors <- c("Simulated" = "black", "Estimated" = "red")
    
    ggplot(mapping = aes(x = date)) + 
      geom_line(aes(y = sim, color = "Simulated")) + 
      geom_line(aes(y = est, color = "Estimated")) +
      geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) + 
      labs(x = "Dates", y = "Futures prices", color = "") +
      scale_color_manual(values = colors) 
  })
  
  # Plot: prices on a given day
  output$plot_price_mat <- renderPlotly({
    yt <- sim_dat()$yt
    yt_hat <- estimation()$yt_hat
    mats <- sim_dat()$mats
    cov_y <- estimation()$cov_y
    
    if (input$model == "SS2000") {
      CI_lower <- qlnorm(0.025, meanlog = log(unlist(yt_hat[input$day, ])), sdlog = sqrt(diag(cov_y[, , input$day])))
      CI_upper <- qlnorm(0.975, meanlog = log(unlist(yt_hat[input$day, ])), sdlog = sqrt(diag(cov_y[, , input$day])))
    } else if (input$model == "PD") {
      CI_lower <- unlist(yt_hat[input$day, ]) - 1.96 * sqrt(diag(cov_y[, , input$day]))
      CI_upper <- unlist(yt_hat[input$day, ]) + 1.96 * sqrt(diag(cov_y[, , input$day]))
    }
    
    date <- mats[input$day, ]
    sim <- unlist(yt[input$day, ])
    est <- unlist(yt_hat[input$day, ])
    colors <- c("Simulated" = "black", "Estimated" = "red")
    
    if (input$n_contract == 1) {
      ggplot() 
    } else {
      ggplot(mapping = aes(x = date)) + 
        geom_line(aes(y = sim, color = "Simulated")) + 
        geom_line(aes(y = est, color = "Estimated")) + 
        geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2) + 
        labs(x = "Time to maturities (in year)", y = "Futures prices", color = "") +
        scale_color_manual(values = colors) 
    }
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
    date <- dat$observations
    re <- dat$re
    
    ggplot(mapping = aes(x = date, y = re)) + 
      geom_line() + 
      xlab("Dates") + 
      ylab("Relative errors")
  })
  
  # Plot: mean absolute errors on each trading date
  output$plot_ae <- renderPlotly({
    err <- errors()
    dat <- data.frame(observations = 1: input$n_obs, ae = apply(err$ae, 1, mean))
    date <- dat$observations
    ae <- dat$ae
    
    ggplot(mapping = aes(x = date, y = ae)) + 
      geom_line() + 
      xlab("Dates") + 
      ylab("Absolute errors")
  })
  
  # Plot: boxplot of relative errors
  output$plot_boxplot_re <- renderPlotly({
    err <- errors()
    dat <- pivot_longer(err$re, 1: input$n_contract)
    dat$name <- factor(dat$name, levels = paste("Contract", 1: input$n_contract, sep = ""))
    ggplot(dat, aes(x = name, y = value, col = name)) + 
      geom_boxplot() + 
      xlab("Contracts") + 
      ylab("Relative errors") + 
      labs(col = "Contracts")
  })
  
  # Plot: boxplot of absolute errors
  output$plot_boxplot_ae <- renderPlotly({
    err <- errors()
    dat <- pivot_longer(err$ae, 1: input$n_contract)
    dat$name <- factor(dat$name, levels = paste("Contract", 1: input$n_contract, sep = ""))
    
    ggplot(dat, aes(x = name, y = value, col = name)) + 
      geom_boxplot() + 
      xlab("Contracts") + 
      ylab("Absolute errors") + 
      labs(col = "Contracts")
  })
}



