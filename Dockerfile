FROM rocker/shiny-verse

FROM rocker/rstudio:4.3.2

RUN mkdir home/PDSim

COPY . home/PDSim

RUN Rscript /home/PDSim/docker_installation.R