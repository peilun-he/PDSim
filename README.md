# Polynomial Diffusion Model Simulation and Estimation (V2.1)
## Introduction
Web application for the polynomial diffusion model. This app generates a futures price data by providing all parameters. Also, it gives state variables and contracts estimations through Extended Kalman Filter (EKF) or Unscented Kalman Filter (UKF). The Schwartz and Smith's two-factor model is also provided for comparison. 

Future plan: 
- Add decomposition of data through "seasonal" package 
- Add forecasting and smoothing 

## Installation
PDSim can be accessed in two ways: 

1. You can use PDSim on the Shiny server. This way, you don't need to have R installed on your computer. Just go to https://peilunhe.shinyapps.io/pdsim/ and use it there. 

2. Additionally, you can download and run PDSim locally, by running the following R code: 

```r
# install.packages("devtools") # uncomment if you do not have devtools installed
devtools::install_github("peilun-he/PDSim")
PDSim::run_app()
```

## User Guide

### Schwartz-Smith Model

### Polynomial Diffusion Model

### Some Other Hints
1. Once users enter all parameters, the data will be generated automatically. Users do NOT need to click any buttons. However, if users wish to generate more realisations under the same set of parameters, please click the 'Generate new data' button.

2. The seed to generate random numbers is fixed, i.e., for the same set of parameters, users will get exactly the same data every time they use PDSim.

3. Futures prices in all tables / plots are REAL prices (NOT the logarithm), no matter which model is used.

4. The 95% confidence interval is shown as a grey ribbon on each plot.

5. Because of the limitation of filtering methods, the standard error of the estimated futures price on the first day is extremely large. All plots of contracts estimation start from the second day.

## Version History 
**Version 2.1** (current version): 
- PDSim is packaged into an R package. Some structures is changed to achieve this. 
- A exported function "run_app" is added to run PDSim. 
- Add some documentation.  

**Version 2.0**: 
- Add navigation bar: welcome page, app, user guide, team members. 
- Descriptions of models and some hints are added to the user guide page. 
- Allow users to download simulated data as csv files. 
- Add 95% confidence intervals to all estimations.  
- Add a 3D surface of data. 
- Allow users to generate new realisations of data using same set of parameters. 
- Bugs fixed. 

**Version 1.0**: basic functions and UI
