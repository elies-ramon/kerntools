# Kernel-target alignment

\`KTA()\` computes the alignment between a kernel matrix and a target
variable.

## Usage

``` r
KTA(K, y)
```

## Arguments

- K:

  A kernel matrix (class: "matrix").

- y:

  The target variable. A factor with two levels.

## Value

Alignment value.

## Examples

``` r
K1 <- RBF(iris[1:100,1:4],g=0.1)
y <- factor(iris[1:100,5])
KTA(K1,y)
#> Remember that Klist should contain only kernel matrices (i.e. squared, symmetric and PSD).
#>   This function does NOT verify the symmetry and PSD criteria.
#> [1] 0.4058431
```
