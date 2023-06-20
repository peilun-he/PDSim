draw_app <- function() {
  tabPanel("App", 
           withMathJax(),
           titlePanel(h1("Polynomial Diffusion Model Simulation and Estimation")), 
           #themeSelector(), 
           sidebarLayout(
             sidebarPanel(
               titlePanel(h3("Model Specification")), 
               fluidRow(
                 column(5, 
                        numericInput("n_obs", "No. of trading days", 100, step = 100, min = 100)
                 ), 
                 column(5, offset = 2,
                        numericInput("n_contract", "No. of contracts", 10, step = 1, min = 1)
                 )
               ), 
               selectInput("model", "Select a model", selected = "PD", 
                           choices = c("Schwartz-Smith two-factor model" = "SS2000", "Polynomial diffusion model" = "PD")),
               conditionalPanel(condition = "input.model == 'PD'", 
                                selectInput("filter", "Select a filtering method", selected = "Extended Kalman Filter", 
                                            choices = c("Extended Kalman Filter", "Unscented Kalman Filter")), 
                                fluidRow(
                                  column(5, 
                                         numericInput("constant", "\\( \\alpha_1 \\)", 1, step = 0.1), 
                                         numericInput("xi", "\\( \\alpha_3 \\)", 1, step = 0.1), 
                                         numericInput("chi_xi", "\\( \\alpha_5 \\)", 1, step = 0.1)
                                  ), 
                                  column(5, offset = 2, 
                                         numericInput("chi", "\\( \\alpha_2 \\)", 1, step = 0.1),
                                         numericInput("chi2", "\\( \\alpha_4 \\)", 1, step = 0.1), 
                                         numericInput("xi2", "\\( \\alpha_6 \\)", 1, step = 0.1)
                                  )
                                )
               ),
               p(actionButton("new_data", "Generate new data", style = "width:10cm"), class = "text-center"), 
               p(downloadButton("download_data", "Download prices", style = "width:5cm"), 
                 downloadButton("download_mat", "Download maturities", style = "width:5cm"), class = "text-center"), 
               hr(),
               
               titlePanel(h3("Estimations")), 
               numericInput("contract", "Enter the contract number you want to show estimations", 
                            value = 1, min = 1, step = 1),
               numericInput("day", "Enter a day you want to show estimations", 
                            value = 2, min = 2, step = 1),
               hr(), 
               
               titlePanel(h3("Model parameters")), 
               fluidRow(
                 column(5, 
                        numericInput("kappa", "\\( \\kappa \\)", 0.5, min = 0.1, step = 0.1), 
                        numericInput("mu", "\\( \\mu_{\\xi} \\)", 1, step = 0.1), 
                        numericInput("sc", "\\( \\sigma_{\\chi} \\)", 1.5, min = 0.1, step = 0.1), 
                        numericInput("lc", "\\( \\lambda_{\\chi} \\)", 0.5, step = 0.1), 
                        numericInput("s1", "\\( \\sigma_1 \\)", 0.1, min = 0.001, step = 0.01)
                 ), 
                 column(5, offset = 2, 
                        numericInput("gamma", "\\( \\gamma \\)", 0.3, min = 0.1, step = 0.1), 
                        numericInput("rho", "\\( \\rho \\)", -0.3, min = -0.99, max = 0.99, step = 0.01), 
                        numericInput("sx", "\\( \\sigma_{\\xi} \\)", 1.3, min = 0.1, step = 0.1), 
                        numericInput("lx", "\\( \\lambda_{\\xi} \\)", 0.3, step = 0.1), 
                        numericInput("sn", "\\( \\sigma_m \\)", 0.01, min = 0.001, step = 0.01)
                 )
               )
             ), 
             mainPanel(
               tabsetPanel(
                 tabPanel("Simulated futures contracts", plotlyOutput("plot_yt")), 
                 tabPanel("3D surface of futures prices", plotlyOutput("plot_3d_yt")), 
                 tabPanel("Simulated state variables", plotlyOutput("plot_xt")), 
                 tabPanel("Simulated futures prices", DTOutput("table_yt")), 
                 tabPanel("Simulated time to maturities (in year)", DTOutput("table_mat"))
               ), 
               tabsetPanel(
                 tabPanel("Simulated \\( \\chi \\) vs estimated \\( \\chi \\)", plotlyOutput("plot_chi")), 
                 tabPanel("Simulated \\( \\xi \\) vs estimated \\( \\xi \\)", plotlyOutput("plot_xi")), 
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
}
