draw_members <- function() {
  tabPanel("Team Members", 
           titlePanel(h1("Team Members")), 
           fluidRow(
             column(5, 
                    imageOutput("ph"), 
                    p(a("Peilun He", href = "https://scholar.google.com/citations?authuser=2&user=ixbORcQAAAAJ", target = "_blank"), style = "font-size:25px"), 
                    p("e-mail: peilun.he93@gmail.com", style = "font-size:20px"),
                    br(), 
                    br(),
                    imageOutput("gp"), 
                    p(a("Prof. Gareth W. Peters", href = "https://www.qrslab.com/", targets = "_blank"), style = "font-size:25px"), 
                    p("e-mail: garethpeters@ucsb.edu", style = "font-size:20px")
                    ), 
             column(5, offset = 2, 
                    imageOutput("nk"), 
                    p(a("Dr. Nino Kordzakhia", href = "https://researchers.mq.edu.au/en/persons/nino-kordzakhia", target = "_blank"), style = "font-size:25px"), 
                    p("email: nino.kordzakhia@mq.edu.au", style = "font-size:20px"), 
                    br(),
                    br(),
                    imageOutput("ps"), 
                    p(a("Prof. Pavel V. Shevchenko", href = "https://www.mq.edu.au/research/research-centres-groups-and-facilities/prosperous-economies/centres/centre-for-risk-analytics/our-people/pavel-shevchenko", target = "_blank"), style = "font-size:25px"),
                    p("e-mail: pavel.shevchenko@mq.edu.au", style = "font-size:20px"))
                    )
           )
}