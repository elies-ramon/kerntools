# Gamma hyperparameter estimation (RBF kernel)

This function returns an estimation of the optimum value for the gamma
hyperparameter (required by the RBF kernel function) using different
heuristics:

- *D* criterion:

  It returns the inverse of the number of features in X.

- Scale criterion:

  It returns the inverse of the number of features, normalized by the
  total variance of X.

- Quantiles criterion:

  A range of values, computed with the function \`kernlab::sigest()\`.

## Usage

``` r
estimate_gamma(X)
```

## Arguments

- X:

  Matrix or data.frame that contains real numbers ("integer", "float" or
  "double").

## Value

A list with the gamma value estimation according to different criteria.

## Examples

``` r
data <- matrix(rnorm(150),ncol=50,nrow=30)
gamma <- estimate_gamma(data)
gamma
#> $d_criterion
#> [1] 0.02
#> 
#> $scale_criterion
#> [1] 0.02031485
#> 
#> $quantiles_criterion
#>         90%         50%         10% 
#> 0.006224311 0.010776361 0.016232140 
#> 
K <- RBF(data, g = gamma$scale_criterion)
K[1:5,1:5]
#>            [,1]      [,2]       [,3]       [,4]       [,5]
#> [1,] 1.00000000 0.0330113 0.01602149 0.49088496 0.02523561
#> [2,] 0.03301130 1.0000000 0.78524768 0.14693995 0.72818199
#> [3,] 0.01602149 0.7852477 1.00000000 0.09096069 0.61614457
#> [4,] 0.49088496 0.1469400 0.09096069 1.00000000 0.14702329
#> [5,] 0.02523561 0.7281820 0.61614457 0.14702329 1.00000000
```
