# Linear kernel

\`Linear()\` computes the inner product between all possible pairs of
rows of a matrix or data.frame with dimension *NxD*.

## Usage

``` r
Linear(X, cos.norm = FALSE, coeff = NULL)
```

## Arguments

- X:

  Matrix or data.frame that contains real numbers ("integer", "float" or
  "double").

- cos.norm:

  Should the resulting kernel matrix be cosine normalized? (Defaults:
  FALSE).

- coeff:

  (optional) A vector of length *D* that weights each one of the
  features (columns). When cos.norm=TRUE, \`Linear()\` first does the
  weighting and then the cosine-normalization.

## Value

Kernel matrix (dimension: *NxN*).

## Examples

``` r
dat <- matrix(rnorm(250),ncol=50,nrow=5)
Linear(dat)
#>              [,1]        [,2]        [,3]       [,4]      [,5]
#> [1,]  56.36893864  0.09496672 -15.2467464 -14.050418 11.868900
#> [2,]   0.09496672 58.56731315   0.7992432  -2.150232  4.124800
#> [3,] -15.24674635  0.79924315  58.9956748  10.422470 -5.152776
#> [4,] -14.05041796 -2.15023246  10.4224699  56.511142  4.236247
#> [5,]  11.86890032  4.12479992  -5.1527758   4.236247 43.785283
```
