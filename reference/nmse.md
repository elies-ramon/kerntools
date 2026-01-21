# NMSE (Normalized Mean Squared Error)

\`nmse()\` computes the Normalized Mean Squared Error between the output
of a regression model and the actual values of the target.

## Usage

``` r
nmse(target, pred)
```

## Arguments

- target:

  Numeric vector containing the actual values.

- pred:

  Numeric vector containing the predicted values. (The order should be
  the same than in the target)

## Value

The normalized mean squared error (a single value).

## Details

The Normalized Mean Squared error is defined as:

\$\$NMSE=MSE/((N-1)\*var(target))\$\$

where MSE is the Mean Squared Error.

## Examples

``` r
y <- 1:10
y_pred <- y+rnorm(10)
nmse(y,y_pred)
#> [1] 0.1255219
```
