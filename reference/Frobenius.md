# Frobenius kernel

\`Frobenius()\` computes the Frobenius kernel between numeric matrices.

## Usage

``` r
Frobenius(DATA, cos.norm = FALSE, feat_space = FALSE)
```

## Arguments

- DATA:

  A list of *M* matrices or data.frames containing only real numbers
  (class "integer", "float" or "double"). All matrices or data.frames
  should have the same number of rows and columns.

- cos.norm:

  Should the resulting kernel matrix be cosine normalized? (Defaults:
  FALSE).

- feat_space:

  If FALSE, only the kernel matrix is returned. Otherwise, the feature
  space is also returned. (Defaults: FALSE).

## Value

Kernel matrix (dimension:*NxN*), or a list with the kernel matrix and
the feature space.

## Details

The Frobenius kernel is the same than the Frobenius inner product
between matrices.

## Examples

``` r
data1 <- matrix(rnorm(250000),ncol=500,nrow=500)
data2 <- matrix(rnorm(250000),ncol=500,nrow=500)
data3 <- matrix(rnorm(250000),ncol=500,nrow=500)

Frobenius(list(data1,data2,data3))
#>               [,1]        [,2]          [,3]
#> [1,] 250153.701922   -197.6619      2.945187
#> [2,]   -197.661925 249707.9675  -1017.169660
#> [3,]      2.945187  -1017.1697 249701.128325
```
