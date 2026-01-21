# Total Sum Scaling

This function transforms a dataset from absolute to relative frequencies
(by row or column).

## Usage

``` r
TSS(X, rows = TRUE)
```

## Arguments

- X:

  Numeric matrix or data.frame of any size containing absolute
  frequencies.

- rows:

  If TRUE, the operation is done by row; otherwise, it is done by
  column. (Defaults: TRUE).

## Value

A relative frequency matrix or data.frame with the same dimension than
X.

## Examples

``` r
dat <- matrix(rnorm(50),ncol=5,nrow=10)
TSS(dat) #It can be checked that, after scaling, the sum of each row is equal to 1.
#>               [,1]        [,2]        [,3]       [,4]        [,5]
#>  [1,]  0.644153574  0.27332252  0.07910102  0.1468468 -0.14342395
#>  [2,]  0.413748746  0.54177195 -0.20498894 -0.2033853  0.45285352
#>  [3,] -0.751742108  0.43041617  0.43017256  0.3982113  0.49294207
#>  [4,] -0.598075679  0.24131708  0.74357827  0.7466196 -0.13343927
#>  [5,]  0.586903709 -0.09728132 -0.01723159  0.1706674  0.35694180
#>  [6,] -0.062023578 -0.14948124  0.65351494  0.2434817  0.31450817
#>  [7,] -0.173923654  0.62158551  0.91337723 -0.7301574  0.36911831
#>  [8,] -0.007139825  0.37273957  1.75723180 -0.4488926 -0.67393899
#>  [9,] -0.079137230  0.52411641  0.22308516 -0.1439274  0.47586309
#> [10,]  0.723565156  0.42500305 -0.35864401  0.2547892 -0.04471336
```
