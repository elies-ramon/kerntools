# Von Neumann entropy

\`vonNeumann()\` computes the von Neumann entropy of a kernel matrix.
Entropy values close to 0 indicate that all its elements are very
similar, which may result in underfitting when training a prediction
model. Instead, values close to 1 indicate a high variability which may
produce overfitting.

## Usage

``` r
vonNeumann(K)
```

## Arguments

- K:

  Kernel matrix (class "matrix").

## Value

Von Neumann entropy (a single value).

## References

Belanche-Muñoz, L.A. and Wiejacha, M. (2023) Analysis of Kernel Matrices
via the von Neumann Entropy and Its Relation to RVM Performances.
Entropy, 25, 154. doi:10.3390/e25010154.
[Link](https://www.mdpi.com/1099-4300/25/1/154)

## Examples

``` r
data <- matrix(rnorm(150),ncol=50,nrow=30)
K <- Linear(data)
vonNeumann(K)
#> [1] 0.4328801
```
