# Minmax normalization

Minmax normalization. Custom min/max values may be passed to the
function.

## Usage

``` r
minmax(X, rows = FALSE, values = NULL)
```

## Arguments

- X:

  Numeric matrix or data.frame of any size.

- rows:

  If TRUE, the minmax normalization is done by row; otherwise, it is
  done by column. (Defaults: FALSE)

- values:

  (optional) A list containing two elements, the "max" values and the
  "min" values. If no value is passed, the typical minmax normalization
  (which normalizes the dataset between 0 and 1) is computed with the
  observed maximum and minimum value in each column (or row) of X.

## Value

Minmax-normalized X.

## Examples

``` r
dat <- matrix(rnorm(100),ncol=10,nrow=10)
dat_minmax <- minmax(dat)
apply(dat_minmax,2,min) ## Min values = 0
#>  [1] 0 0 0 0 0 0 0 0 0 0
apply(dat_minmax,2,max) ## Max values = 1
#>  [1] 1 1 1 1 1 1 1 1 1 1
# We can also explicitly state the max and min values:
values <- list(min=apply(dat,2,min),max=apply(dat,2,max))
dat_minmax <- minmax(dat,values=values)
```
