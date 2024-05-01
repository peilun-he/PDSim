draw_unit_tests <- function() {
  style = "font-size:25px"
  tabPanel("Unit Tests", 
           withMathJax(),
           titlePanel(h1("Unit Tests")),
           p("Users can undergo a unit test here to ensure that all functionalities 
             of PDSim are operating correctly. This test sequence entails several 
             key steps: initially, users define the desired number of trajectories and 
             relevant parameters. Subsequently, PDSim executes simulations based on 
             these specifications, generating simulated trajectories. Upon simulation 
             completion, we employ KF/EKF/UKF methodologies to estimate trajectories 
             alongside their 95% confidence intervals. The coverage rate, indicating 
             the proportion of trajectories where over 95% of points fall within the 
             confidence interval, is then computed. Given our knowledge of the true 
             parameter values, a high coverage rate, ideally close to 100%, is expected.", 
             style = style), 
           p("Users receive detailed feedback under the 'Results' tab panel. Moreover, 
             PDSim generates two plots: one illustrating the trajectory with the highest 
             coverage rate and another depicting the trajectory with the lowest coverage 
             rate. Additionally, a table presents the coverage rate for each trajectory.",
             style = style),
           p("It's important to note a few considerations: firstly, for simplicity, only 
             a single contract is simulated, regardless of the number specified by the user. 
             Secondly, if the coverage rate falls below 99%, users are advised to either 
             increase the number of trajectories or adjust parameters. Lastly, users are 
             informed that extensive simulations may lead to longer processing times; 
             for instance, generating results for 100 trajectories typically requires around 
             15 seconds on a standard laptop.",
             style = style), 
           titlePanel(h3("Please specify your model in the tab panel 'App', 
                         specify the number of trajectories below, and then press the 
                         'Start unit test' button to start a unit test.")),
           sidebarLayout(
             sidebarPanel(
               numericInput("n_tra", "No. of trajectories", 1, step = 10, min = 1), 
               p(actionButton("start", "Start unit test", style = "width:10cm"), class = "text-center")
             ), 
             mainPanel(
               tabsetPanel(
                 tabPanel("Results", textOutput("unit_tests_msg")), 
                 tabPanel("The best trajectory", plotlyOutput("plot_best_unit_test")),
                 tabPanel("The worst trajectory", plotlyOutput("plot_worst_unit_test")),
                 tabPanel("Coverage rate for each trajectory", DTOutput("table_coverage"))
               )
             )
          )
  )
}