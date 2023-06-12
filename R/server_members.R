# Server logic for the team members page
server_members <- function() {
  # Image: Peilun He
  output$ph <- renderImage({
    list(src = system.file("extdata/peilun.jpg", package = "PDSim"), width = 400, height = 400)
  }, deleteFile = FALSE)
  
  # Image: Nino Kordzakhia
  output$nk <- renderImage({
    list(src = system.file("extdata/nino.jpg", package = "PDSim"), width = 300, height = 400)
  }, deleteFile = FALSE)
  
  # Image: Gareth Peters
  output$gp <- renderImage({
    list(src = system.file("extdata/gareth.jpg", package = "PDSim"), width = 400, height = 400)
  }, deleteFile = FALSE)
  
  # Image: Pavel Shevchenko
  output$ps <- renderImage({
    list(src = system.file("extdata/pavel.jpg", package = "PDSim"), width = 300, height = 400)
  }, deleteFile = FALSE)
}