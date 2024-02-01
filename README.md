# popdynmodel

<!-- badges: start -->
<!-- badges: end -->
*popdynmodel* provides tools to implement population dynamics models in a hierarchical Bayesian framework. Models are intended for occupancy, abundance and/or biomass time series and allow for estimations of occupancy change rates and population growth rates at various spatial and temporal scales. 

*popdynmodel* was initially developed to analyse the river fish monitoring data from the ASPE database managed by the French Biodiversity Agency (OFB), but it can be used for a wide variety of taxa.

## Installation

*popdynmodel* requires the following packages to be installed: *dplyr*, *magrittr*, *MCMCvis*, *nimble*, *rlang*, *stringr*; *tidyr* and *tidyverse*. You can install the missing packages from CRAN using the following command:

``` r 
install.packages(c("coda","rjags","R2jags","stringr","magrittr","dplyr"))
```

Then you can install the lasted development version of *popdynmodel* from Github using :

``` r
library(devtools)
install_github("manue6/popdynmodel")
```

You can also install from the downloadable tar archive [*popdynmodel_1.0.tar.gz*](https://github.com/manue6/popdynmodel_tuto) using :

``` r
install.packages("popdynmodel_1.0.tar.gz", repos = NULL)
```
