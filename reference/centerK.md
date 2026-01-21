# Centering a kernel matrix

It is equivalent to compute \`K\` over centered data (i.e. the mean of
each column is subtracted) in Feature Space.

## Usage

``` r
centerK(K)
```

## Arguments

- K:

  Kernel matrix (class "matrix").

## Value

Centered \`K\` (class "matrix").

## Examples

``` r
dat <- matrix(rnorm(250),ncol=50,nrow=5)
K <- Linear(dat)
centerK(K)
#>            [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,]  36.683262 -13.383489 -12.032447  -9.653916  -1.613410
#> [2,] -13.383489  45.009171 -12.265942 -11.884602  -7.475137
#> [3,] -12.032447 -12.265942  37.211863  -4.253565  -8.659909
#> [4,]  -9.653916 -11.884602  -4.253565  43.311665 -17.519581
#> [5,]  -1.613410  -7.475137  -8.659909 -17.519581  35.268038
```
