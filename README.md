
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mtsp

Genetic algorithm to solve multiple travelling salespersons. Port of
matlab code found at
<https://uk.mathworks.com/matlabcentral/fileexchange/31814-mdmtspv_ga-multiple-depot-multiple-traveling-salesmen-problem-solved-by-genetic-algorithm>
and
<https://uk.mathworks.com/matlabcentral/fileexchange/19049-multiple-traveling-salesmen-problem-genetic-algorithm>.

Install package dependencies:

``` r
pkgs = c("Rcpp", "RcppArmadillo", "inline", "remotes" )
install.packages(pkgs)
```

Then install `mtsp` package from github:

``` r
remotes::install_github("daffp/mtsp")
```

Run some examples to see that it is working

``` r
library(mtsp)
set.seed(1)

# Matrix of node positions to visit
n = 50
xy = matrix(rnorm(n*2), ncol=2)

# Run search
run = mtsp(xy, nSalesmen=5, CostType=2, popSize=80, numIter=10, algorithm="mtsp_ga", return_all=TRUE)
run
#> Total distance travelled for all salespersons = 67.06
#> Maximum distance travelled by a single salesperson = 14.15

summary(run)
#> Total distance travelled: 67.06
#> ================================
#> Salesperson 1 : distance travelled = 11.91
#> Tour: 45 -> 20 -> 14 -> 35 -> 25 -> 8 -> 26 -> 45
#> 
#> Salesperson 2 : distance travelled = 14.15
#> Tour: 24 -> 30 -> 34 -> 41 -> 36 -> 46 -> 32 -> 23 -> 5 -> 1 -> 17 -> 24
#> 
#> Salesperson 3 : distance travelled = 13.92
#> Tour: 28 -> 10 -> 3 -> 31 -> 40 -> 12 -> 29 -> 50 -> 4 -> 44 -> 27 -> 28
#> 
#> Salesperson 4 : distance travelled = 13.41
#> Tour: 19 -> 15 -> 13 -> 37 -> 42 -> 7 -> 39 -> 6 -> 9 -> 18 -> 43 -> 19
#> 
#> Salesperson 5 : distance travelled = 13.67
#> Tour: 16 -> 49 -> 48 -> 2 -> 38 -> 11 -> 21 -> 33 -> 22 -> 47 -> 16
```
