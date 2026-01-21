# Centering a squared matrix by row or column

It centers a numeric matrix with dimension *N x N* by row (rows=TRUE) or
column (rows=FALSE).

## Usage

``` r
centerX(X, rows = TRUE)
```

## Arguments

- X:

  Numeric matrix or data.frame of any size.

- rows:

  If TRUE, the operation is done by row; otherwise, it is done by
  column. (Defaults: TRUE).

## Value

Centered X (class "matrix").

## Examples

``` r
dat <- matrix(rnorm(25),ncol=5,nrow=5)
centerX(dat)
#>            [,1]        [,2]        [,3]        [,4]      [,5]
#> [1,]  0.4813067  0.30514811 -0.92704331 -0.03199894 0.1725874
#> [2,] -0.1955920 -0.03207894 -0.44163425 -0.89187582 1.5611810
#> [3,] -1.1213657 -0.41059101 -0.06197111  0.90113236 0.6927955
#> [4,] -0.8567685  0.55807998  0.26975244 -0.36846556 0.3974016
#> [5,] -0.6740542  1.08171252 -0.27616077 -0.59329251 0.4617950
```
