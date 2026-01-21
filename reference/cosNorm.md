# Cosine normalization of a kernel matrix

It is equivalent to compute K using the normalization
\`X/sqrt(sum(X^2))\` in Feature Space.

## Usage

``` r
cosNorm(K)
```

## Arguments

- K:

  Kernel matrix (class "matrix").

## Value

Cosine-normalized K (class "matrix").

## References

Ah-Pine, J. (2010). Normalized kernels as similarity indices. In
Advances in Knowledge Discovery and Data Mining: 14th Pacific-Asia
Conference, PAKDD 2010, Hyderabad, India, June 21-24, 2010. Proceedings.
Part II 14 (pp. 362-373). Springer Berlin Heidelberg.
[Link](https://hal.science/hal-01504523/document)

## Examples

``` r
dat <- matrix(rnorm(250),ncol=50,nrow=5)
K <- Linear(dat)
cosNorm(K)
#>             [,1]        [,2]        [,3]        [,4]        [,5]
#> [1,]  1.00000000 -0.07455420  0.04626991  0.14706929  0.10853172
#> [2,] -0.07455420  1.00000000 -0.09079827  0.06999559  0.19282377
#> [3,]  0.04626991 -0.09079827  1.00000000  0.18878683 -0.12107437
#> [4,]  0.14706929  0.06999559  0.18878683  1.00000000 -0.03161472
#> [5,]  0.10853172  0.19282377 -0.12107437 -0.03161472  1.00000000
```
