# This file is used for building Docker image only. 

install.packages(c('DT', 
                   'ggplot2', 
                   'lubridate', 
                   'plotly', 
                   'scales', 
                   'shiny', 
                   'shinythemes', 
                   'tidyr'))

install.packages('/home/PDSim', repos = NULL, type="source")
