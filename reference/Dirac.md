# Kernels for categorical variables

From a matrix or data.frame with dimension *NxD*, where *N*\>1, *D*\>0,
\`Dirac()\` computes the simplest kernel for categorical data. Samples
should be in the rows and features in the columns. When there is a
single feature, \`Dirac()\` returns 1 if the category (or class, or
level) is the same in two given samples, and 0 otherwise. Instead, when
*D*\>1, the results for the *D* features are combined doing a sum, a
mean, or a weighted mean.

## Usage

``` r
Dirac(X, comp = "mean", coeff = NULL, feat_space = FALSE)
```

## Arguments

- X:

  Matrix (class "character") or data.frame (class "character", or
  columns = "factor"). The elements in X are assumed to be categorical
  in nature.

- comp:

  When *D*\>1, this argument indicates how the variables of the dataset
  are combined. Options are: "mean", "sum" and "weighted". (Defaults:
  "mean")

  - "sum" gives the same importance to all variables, and returns an
    unnormalized kernel matrix.

  - "mean" gives the same importance to all variables, and returns a
    normalized kernel matrix (all its elements range between 0 and 1).

  - "weighted" weights each variable according to the \`coeff\`
    parameter, and returns a normalized kernel matrix.

- coeff:

  (optional) A vector of weights with length *D*.

- feat_space:

  If FALSE, only the kernel matrix is returned. Otherwise, the feature
  space is also returned. (Defaults: FALSE).

## Value

Kernel matrix (dimension: *NxN*), or a list with the kernel matrix and
the feature space.

## References

Belanche, L. A., and Villegas, M. A. (2013). Kernel functions for
categorical variables with application to problems in the life sciences.
Artificial Intelligence Research and Development (pp. 171-180). IOS
Press.
[Link](https://upcommons.upc.edu/bitstream/handle/2117/23347/KernelCATEG_CCIA2013.pdf)

## Examples

``` r
# Categorical data
summary(CO2)
#>      Plant             Type         Treatment       conc          uptake     
#>  Qn1    : 7   Quebec     :42   nonchilled:42   Min.   :  95   Min.   : 7.70  
#>  Qn2    : 7   Mississippi:42   chilled   :42   1st Qu.: 175   1st Qu.:17.90  
#>  Qn3    : 7                                    Median : 350   Median :28.30  
#>  Qc1    : 7                                    Mean   : 435   Mean   :27.21  
#>  Qc3    : 7                                    3rd Qu.: 675   3rd Qu.:37.12  
#>  Qc2    : 7                                    Max.   :1000   Max.   :45.50  
#>  (Other):42                                                                  
Kdirac <- Dirac(CO2[,1:3])
## Display a subset of the kernel matrix:
Kdirac[c(1,15,50,65),c(1,15,50,65)]
#>            1        15        50        65
#> 1  1.0000000 0.6666667 0.3333333 0.0000000
#> 15 0.6666667 1.0000000 0.3333333 0.0000000
#> 50 0.3333333 0.3333333 1.0000000 0.3333333
#> 65 0.0000000 0.0000000 0.3333333 1.0000000
```
