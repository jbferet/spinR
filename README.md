---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- ```{r include = FALSE} -->
<!-- knitr::opts_chunk$set( -->
<!--   collapse = TRUE, -->
<!--   comment = "#>", -->
<!--   fig.path = "man/figures/README-", -->
<!--   out.width = "100%" -->
<!-- ) -->
<!-- ``` -->

# prospect

<!-- badges: start -->
<!-- badges: end -->

The goal of prospect is to ...

## Installation

You can install the released version of prospect from [gitlab](https://gitlab.com/jbferet/prospect) with:

``` r
devtools::install_gitlab('jbferet/prospect')
```

## Example

This is a basic example which shows you how to simulate leaf optical properties using __prospect__:

``` r
library(prospect)
#############################################################
###         run PROSPECT using default parameters         ###
#############################################################
LRT   = PROSPECT(SpecPROSPECT)

#############################################################
###   run PROSPECT using user defind set of parameters    ###
#############################################################
LRT2   = PROSPECT(SpecPROSPECT,N = 1.4,CHL = 30,CAR = 10,EWT = 0.02,LMA = 0.01)
```