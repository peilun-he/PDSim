# Polynomial diffusion model simulation and estimation (V2.0)
Web application for the polynomial diffusion model. This app generates a futures price data by providing all parameters. Also, it gives state variables and contracts estimations through Extended Kalman Filter (EKF) or Unscented Kalman Filter (UKF). The Schwartz and Smith's two-factor model is also provided for comparison. 

Future plan: 
- Add decomposition of data through "seasonal" package 
- Add forecasting and smoothing 

## Installation
PDSim can be accessed in two ways. 

1. You can access PDSim at: https://peilunhe.shinyapps.io/pdsim/. 

2. Additionally, you can download and run PDSim locally, by running the following R code: 

```r
# install.packages("devtools") # uncomment if you do not have devtools installed
devtools::install_github("peilun-he/PDSim")
PDSim::run_app()
```

## Version history 
**Version 2.0** (current version): 
- Add navigation bar: welcome page, app, user guide, team members
- Descriptions of models and some hints are added to the user guide page
- Allow users to download simulated data as csv files
- Add 95% confidence intervals to all estimations 
- Add a 3D surface of data
- Allow users to generate new realisations of data using same set of parameters
- Bugs fixed 

**Version 1.0**: basic functions and UI
