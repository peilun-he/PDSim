# Server logic for the welcome page
# Image: home image
server_welcome <- function() {
  output$home_img <- renderImage({
    list(src = "Images/Home.jpg", width = "100%", height = "200%")
  }, deleteFile = FALSE)
}