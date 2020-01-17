
<!-- README.md is generated from README.Rmd. Please edit that file -->
uroc
====

<!-- badges: start -->
<!-- badges: end -->
The uroc package provides the functionality of creating an animated ROC movie (ROCM), a universal ROC (UROC) curve and to compute the coefficient of predictive ability (CPA). These tools generalize the classical ROC curve and AUC and can be applied to assess the predictive abilities of features, markers and tests for not only binary classification problems but for just any ordinal or real-valued outcome.
For more information see: https://arxiv.org/abs/1912.01956

Installation
------------

You can install the latest version of uroc from Github with:

``` r
library(devtools)
install_github("evwalz/uroc")
```

Example
-------

The following basic example shows how to create a UROC curve:

``` r
library(uroc)

# Use build-in data set in R: Longley's Economic Regression Data 
data(longley)

# define response and predictor variable: Use gross national product (GNP) as feaure marker to predict the number of employed people (Employed)
response <- longley$Employed
predictor <- longley$GNP

# compute coefficient of predictive ability
cpa(response, predictor)

# create uROC curve
uroc(response, predictor)
```

