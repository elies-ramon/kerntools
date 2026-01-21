# Gaussian RBF (Radial Basis Function) kernel

\`RBF()\` computes the RBF kernel between all possible pairs of rows of
a matrix or data.frame with dimension *NxD*.

## Usage

``` r
RBF(X, g = NULL)
```

## Arguments

- X:

  Matrix or data.frame that contains real numbers ("integer", "float" or
  "double").

- g:

  Gamma hyperparameter. If g=0 or NULL, \`RBF()\` returns the matrix of
  squared Euclidean distances instead of the RBF kernel matrix.
  (Defaults=NULL).

## Value

Kernel matrix (dimension: *NxN*).

## Details

Let \\x_i,x_j\\ be two real vectors. Then, the RBF kernel is defined as:
\$\$K\_{RBF}(x_i,x_j)=\exp(-\gamma \\x_i - x_j \\^2)\$\$

Sometimes the RBF kernel is given a hyperparameter called sigma. In that
case: \\\gamma = 1/\sigma^2\\.

## Examples

``` r
dat <- matrix(rnorm(250),ncol=50,nrow=5)
RBF(dat,g=0.1)
#>              [,1]         [,2]         [,3]         [,4]         [,5]
#> [1,] 1.000000e+00 2.601581e-05 7.283653e-04 0.0015929600 6.536862e-05
#> [2,] 2.601581e-05 1.000000e+00 6.283656e-05 0.0002441566 5.776624e-05
#> [3,] 7.283653e-04 6.283656e-05 1.000000e+00 0.0008188197 3.386474e-04
#> [4,] 1.592960e-03 2.441566e-04 8.188197e-04 1.0000000000 3.185288e-04
#> [5,] 6.536862e-05 5.776624e-05 3.386474e-04 0.0003185288 1.000000e+00
```
