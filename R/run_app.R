#' Starts the PDSim app
#' 
#' Web application for the polynomial diffusion model. 
#' This app generates a futures price data by providing all parameters. 
#' Also, it gives state variables and contracts estimations through 
#' Extended Kalman Filter (EKF) or Unscented Kalman Filter (UKF). 
#' The Schwartz and Smith's two-factor model is also provided for comparison. 
#' 
#' @import shiny lubridate ggplot2 DT plotly shinythemes tidyr scales MASS
#' @export
#' @examples 
#' \dontrun{
#' library(PDSim)
#' PDSim::run_app()
#' }

run_app <- function() {
  shinyApp(ui = my_ui, server = my_server)
}
