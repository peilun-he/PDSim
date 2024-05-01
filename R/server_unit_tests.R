# Server logic for the unit tests page
server_unit_tests <- function() {
  # Massege for unit test
  output$unit_tests_msg <- renderText({
    coverage_rate <- unit_test()$coverage_rate$Coverage.rate
    mean_coverage_rate <- mean(coverage_rate)
    
    if (mean_coverage_rate > 0.9) {
      msg <- paste("Congratulations! The mean coverage rate for this set of parameters is ", 
                   round(mean_coverage_rate*100, 2), 
                   "%, unit test is successful!", sep = "")
    } else {
      msg <- paste("The mean coverage is ", 
                   round(mean_coverage_rate*100, 2), 
                   ", the unit test is failed. Please try a different set of parameters.", sep = "")
    }
    return(msg)
  })
  
  # Plot for the best trajectory
  output$plot_best_unit_test <- renderPlotly({
    best <- unit_test()$best
    
    colors <- c("Simulated" = "black", "Estimated" = "red")
    
    ggplot(mapping = aes(x = 2: input$n_obs)) + 
      geom_line(aes(y = best$yt[2: input$n_obs], color = "Simulated")) + 
      geom_line(aes(y = best$yt_hat[2: input$n_obs], color = "Estimated")) +
      geom_ribbon(aes(ymin = best$CI_lower[2: input$n_obs], ymax = best$CI_upper[2: input$n_obs]), alpha = 0.2) + 
      labs(x = "Dates", y = "Futures prices", color = "") +
      scale_color_manual(values = colors) 
  })
  
  # Plot for the worst trajectory
  output$plot_worst_unit_test <- renderPlotly({
    worst <- unit_test()$worst
    
    colors <- c("Simulated" = "black", "Estimated" = "red")
    
    ggplot(mapping = aes(x = 2: input$n_obs)) + 
      geom_line(aes(y = worst$yt[2: input$n_obs], color = "Simulated")) + 
      geom_line(aes(y = worst$yt_hat[2: input$n_obs], color = "Estimated")) +
      geom_ribbon(aes(ymin = worst$CI_lower[2: input$n_obs], ymax = worst$CI_upper[2: input$n_obs]), alpha = 0.2) + 
      labs(x = "Dates", y = "Futures prices", color = "") +
      scale_color_manual(values = colors) 
  })
  
  # Table for the coverage rate of each trajectory
  output$table_coverage <- renderDT({
    datatable(
      round(unit_test()$coverage_rate, 4), 
      options = list(
        columnDefs = list( list(className = 'dt-center', targets = "_all") )
      ) )
  })
}