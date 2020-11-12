
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fastdivrate

<!-- badges: start -->
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
#>    0 | 00000000000000000000000000000000000000000000000000000000000000000000+6532
#>    2 | 00000000000000000000000000000000000000000000000000000000000000000000+1952
#>    4 | 00000000000000000000000000000000000000000001111111111111111111111111+625
#>    6 | 00000000000000000000000011111111111111111112222222222222222222222222+246
#>    8 | 00000000011111111111112222222233344444444444555555666667777777888889+52
#>   10 | 00001122334455555566666666688888888888800000011112233333344446666667
#>   12 | 0000111111223333444555566666778899990002222334666688
#>   14 | 0011111222266788000778899
#>   16 | 1122447733888899
#>   18 | 2233556666
#>   20 | 177
#>   22 | 33
#>   24 | 4444
#>   26 | 
#>   28 | 7700
#>   30 | 88
```
