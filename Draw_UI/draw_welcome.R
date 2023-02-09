draw_welcome <- function() {
  tabPanel("Welcome!", icon = icon("home", lib = "glyphicon"), 
           imageOutput("home_img")
           )
}