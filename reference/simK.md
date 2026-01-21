# Kernel matrix similarity

\`simK()\` computes the similarity between kernel matrices.

## Usage

``` r
simK(Klist)
```

## Arguments

- Klist:

  A list of *M* kernel matrices with identical *NxN* dimension.

## Value

Kernel matrix (dimension: *MxM*).

## Details

It is a wrapper of \`Frobenius()\`.

## Examples

``` r
K1 <- Linear(matrix(rnorm(7500),ncol=150,nrow=50))
K2 <- Linear(matrix(rnorm(7500),ncol=150,nrow=50))
K3 <- Linear(matrix(rnorm(7500),ncol=150,nrow=50))

simK(list(K1,K2,K3))
#> Remember that Klist should contain only kernel matrices (i.e. squared, symmetric and PSD).
#>   This function does NOT verify the symmetry and PSD criteria.
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.7459548 0.7449204
#> [2,] 0.7459548 1.0000000 0.7564083
#> [3,] 0.7449204 0.7564083 1.0000000
```
