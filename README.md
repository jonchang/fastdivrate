
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastdivrate

<!-- badges: start -->

[![R-CMD-check](https://github.com/jonchang/fastdivrate/workflows/R-CMD-check/badge.svg)](https://github.com/jonchang/fastdivrate/actions)
<!-- badges: end -->

The goal of fastdivrate is to quickly compute tip-specific
diversification rates using a variety of statistics.

## Installation

Install the development version with:

``` r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("jonchang/fastdivrate")
```

Once this is accepted on CRAN, you could install the released version of
fastdivrate from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("fastdivrate")
```

(Note: Not on CRAN just yet.)

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(fastdivrate)
library(ape)
# Simulate a big tree
tree <- rcoal(10000)
# Calculate its diversification rates, quickly
rates <- DR_statistic_C(tree)
stem(rates)
#> 
#>   The decimal point is 3 digit(s) to the right of the |
#> 
#>    0 | 00000000000000000000000000000000000000000000000000000000000000000000+6352
#>    2 | 00000000000000000000000000000000000000000000000000000000000000000000+2160
#>    4 | 00000000000000000000000000000000000000000000001111111111111111111111+641
#>    6 | 00000000000000000000000000001111111111111111111122222222222222222222+265
#>    8 | 00001111111111111122333333333444455666666666667777777888888999999990+44
#>   10 | 0000002222223344446666666666778899990000002233446666888899
#>   12 | 11112222222234488889911334444555566
#>   14 | 22336688884455
#>   16 | 44555566880
#>   18 | 00334455779977
#>   20 | 88
#>   22 | 77
#>   24 | 
#>   26 | 
#>   28 | 44
```
