# mtsp
Genetic algorithm to solve multiple travelling salespersons. Port of matlab code found at https://uk.mathworks.com/matlabcentral/fileexchange/31814-mdmtspv_ga-multiple-depot-multiple-traveling-salesmen-problem-solved-by-genetic-algorithm and https://uk.mathworks.com/matlabcentral/fileexchange/19049-multiple-traveling-salesmen-problem-genetic-algorithm.


Install package dependencies:

```
pkgs = c("Rcpp",
            "RcppArmadillo",
            "inline",
            "remotes"
            )
for(PKG in pkgs) {
  if(!suppressWarnings(suppressPackageStartupMessages(
        require(PKG, character.only=TRUE, quietly=TRUE)))) install.packages(PKG)
        suppressPackageStartupMessages(library(PKG, character.only=TRUE))
}
```

then install `mtsp` package

```
remotes::install_github("daffp/mtsp)
```

Run some examples to see that it is working

```
library(mtsp)
set.seed(1)
n = 50
xy = matrix(rnorm(n*2), ncol=2)
run = mtsp(xy, nSalesmen=5, CostType=2, popSize=80, numIter=10,
algorithm="mtsp_ga", return_all=TRUE)
run
summary(run)
```