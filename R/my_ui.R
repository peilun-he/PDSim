ver <<- "2.1.2"

my_ui <- fluidPage(
  theme = shinytheme("cyborg"), 
  navbarPage(title = paste("PDSim (v", ver, ")", sep = ""), 
             draw_welcome(), # draw welcome page
             draw_app(), # draw app page
             draw_user_guide(), # draw user guide page
             draw_members() # draw team members page
  )
)