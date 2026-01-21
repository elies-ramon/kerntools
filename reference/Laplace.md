# Laplacian kernel

\`Laplace()\` computes the laplacian kernel between all possible pairs
of rows of a matrix or data.frame with dimension *NxD*.

## Usage

``` r
Laplace(X, g = NULL)
```

## Arguments

- X:

  Matrix or data.frame that contains real numbers ("integer", "float" or
  "double").

- g:

  Gamma hyperparameter. If g=0 or NULL, \`Laplace()\` returns the
  Manhattan distance (L1 norm between two vectors). (Defaults=NULL)

## Value

Kernel matrix (dimension: *NxN*).

## Details

Let \\x_i,x_j\\ be two real vectors. Then, the laplacian kernel is
defined as: \$\$K\_{Lapl}(x_i,x_j)=\exp(-\gamma \\x_i - x_j \\\_1)\$\$

## Examples

``` r
dat <- matrix(rnorm(250),ncol=50,nrow=5)
Laplace(dat,g=0.1)
#>             1           2           3            4            5
#> 1 1.000000000 0.002934677 0.013209824 0.0055161209 0.0038626609
#> 2 0.002934677 1.000000000 0.004009897 0.0024804715 0.0026005123
#> 3 0.013209824 0.004009897 1.000000000 0.0042107338 0.0031150382
#> 4 0.005516121 0.002480471 0.004210734 1.0000000000 0.0009941164
#> 5 0.003862661 0.002600512 0.003115038 0.0009941164 1.0000000000
```
