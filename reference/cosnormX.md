# Cosine normalization of a matrix

Normalizes a numeric matrix dividing each row (if rows=TRUE) or column
(if rows=FALSE) by their L2 norm. Thus, each row (or column) has unit
norm.

## Usage

``` r
cosnormX(X, rows = TRUE)
```

## Arguments

- X:

  Numeric matrix or data.frame of any size.

- rows:

  If TRUE, the operation is done by row; otherwise, it is done by
  column. (Defaults: TRUE).

## Value

Cosine-normalized X.

## Examples

``` r
dat <- matrix(rnorm(50),ncol=5,nrow=10)
cosnormX(dat)
#>              [,1]        [,2]         [,3]       [,4]         [,5]
#>  [1,]  0.71721212  0.01907978  0.569199016 -0.3812222  0.126193741
#>  [2,]  0.51040596 -0.22596499 -0.537773478 -0.4847884  0.405222730
#>  [3,]  0.01968701 -0.78920435  0.439350374  0.2174163  0.369418882
#>  [4,] -0.52545035  0.61017754  0.121323522 -0.3607343  0.454683015
#>  [5,]  0.44029946 -0.40900940 -0.683972537 -0.4135165  0.005779533
#>  [6,] -0.48661733  0.20868852  0.792100883 -0.2766428  0.125290074
#>  [7,] -0.01082222  0.25254037 -0.004292825 -0.8967659  0.363178560
#>  [8,]  0.29817221  0.50761955 -0.612241575 -0.5277285  0.008863592
#>  [9,] -0.15562403  0.30856696  0.261952796  0.5782168  0.691095965
#> [10,]  0.69492315  0.15075159  0.523291311 -0.1909413 -0.429026090
```
