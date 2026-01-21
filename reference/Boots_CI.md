# Confidence Interval using Bootstrap

\`Boots_CI()\` computes the Confidence Interval (CI) of a performance
measure (for instance, accuracy) via bootstrapping.

## Usage

``` r
Boots_CI(target, pred, index = "acc", nboots, confidence = 95, ...)
```

## Arguments

- target:

  Numeric vector containing the actual values.

- pred:

  Numeric vector containing the predicted values. (The order should be
  the same than the target's).

- index:

  Performance measure name, in lowercase. (Defaults: "acc").

- nboots:

  Number of bootstrapping replicas.

- confidence:

  Confidence level; for instance, 95% or 99%. (Defaults: 95).

- ...:

  Further arguments to be passed to the performance measures functions;
  notably, multi.class="macro" or multi.class="micro" for the macro or
  micro performance measures. (Defaults: "macro").

## Value

A vector containing the bootstrap estimate of the performance and its
CI.

## Examples

``` r
y <- c(rep("a",30),rep("b",20))
y_pred <- c(rep("a",20),rep("b",30))
# Computing Accuracy with their 95%CI
Boots_CI(target=y, pred=y_pred, index="acc", nboots=1000, confidence=95)
#>            2.5%   97.5% 
#> 0.79892 0.68000 0.90000 
```
