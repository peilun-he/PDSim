# Server logic for the welcome page
# Image: home image
server_welcome <- function() {
  output$home_img <- renderImage({
    list(src = system.file("extdata/Home.jpg", package = "PDSim"), width = "100%", height = "200%")
  }, deleteFile = FALSE)
}