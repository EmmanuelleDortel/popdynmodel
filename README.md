# popdynmodel

<!-- badges: start -->
<!-- badges: end -->
*popdynmodel* provides tools and functions to implement population dynamics models in a hierarchical Bayesian framework. The models are intended for occupancy, abundance and/or biomass time series. They allow the estimation of occupancy change rates and population growth rates, for one or more taxa, at various spatial and temporal scales. This package incudes different types of functions defined by a specific prefixe:
- get : getting information about package modelling functions
- mef : dataframes formatting (useful for meeting modelling functions requirements)
- mod : simply statistical relationship-based modelling
- modenv : environmental relationship-based modelling
- wri : writing model and model data


*popdynmodel* was initially developed to analyse the river fish monitoring data from the ASPE database managed by the French Biodiversity Agency (OFB), but it can be used for a wide variety of taxa.

# Installation

*popdynmodel* requires the following packages to be installed: *dplyr*, *coda*, *magrittr*, *MCMCvis*, *nimble*, *rlang*, *stringr*; *tidyr* and *tidyverse*. You can install the missing packages from CRAN using the following command:

``` r 
install.packages(c("dplyr", "coda", "magrittr", "MCMCvis", "nimble", "rlang", "stringr", "tidyr", "tidyverse"))
```

Then you can install the lasted development version of *popdynmodel* from Github using :

``` r
library(devtools)
install_github("mmanuelleDortel/popdynmodel")
```

