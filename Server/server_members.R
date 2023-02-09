# Server logic for the team members page
server_members <- function() {
  # Image: Peilun He
  output$ph <- renderImage({
    list(src = "Images/peilun.jpg", width = 400, height = 400)
  }, deleteFile = FALSE)
  
  # Image: Nino Kordzakhia
  output$nk <- renderImage({
    list(src = "Images/nino.jpg", width = 300, height = 400)
  }, deleteFile = FALSE)
  
  # Image: Gareth Peters
  output$gp <- renderImage({
    list(src = "Images/gareth.jpg", width = 400, height = 400)
  }, deleteFile = FALSE)
  
  # Image: Pavel Shevchenko
  output$ps <- renderImage({
    list(src = "Images/pavel.jpg", width = 300, height = 400)
  }, deleteFile = FALSE)
}