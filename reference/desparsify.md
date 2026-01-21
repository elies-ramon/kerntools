# This function deletes those columns and/or rows in a matrix/data.frame that only contain 0s.

This function deletes those columns and/or rows in a matrix/data.frame
that only contain 0s.

## Usage

``` r
desparsify(X, dim = 2)
```

## Arguments

- X:

  Numeric matrix or data.frame of any size.

- dim:

  A numeric vector. 1 indicates that the function should be applied to
  rows, 2 to columns, c(1, 2) indicates rows and columns. (Defaults: 2).

## Value

X with less rows or columns. (Class: the same than X).

## Examples

``` r
dat <- matrix(rnorm(150),ncol=50,nrow=30)
dat[c(2,6,12),] <- 0
dat[,c(30,40,50)] <- 0
dim(desparsify(dat))
#> [1] 30 47
dim(desparsify(dat,dim=c(1,2)))
#> [1] 27 47
```
