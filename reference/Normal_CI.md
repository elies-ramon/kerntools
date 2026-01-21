# Confidence Interval using Normal Approximation

\`Normal_CI()\` computes the Confidence Interval (CI) of a performance
measure (for instance, accuracy) using normal approximation. Thus, it is
advisable that the test has a size of, at least, 30 instances.

## Usage

``` r
Normal_CI(value, ntest, confidence = 95)
```

## Arguments

- value:

  Performance value (a single value).

- ntest:

  Test set size (a single value).

- confidence:

  Confidence level; for instance, 95% or 99%. (Defaults: 95).

## Value

A vector containing the CI.

## Examples

``` r
# Computing accuracy
y <- c(rep("a",30),rep("b",20))
y_pred <- c(rep("a",20),rep("b",30))
ct <- table(y,y_pred)
accuracy <- Acc(ct)
# Computing 95%CI
Normal_CI(accuracy, ntest=length(y), confidence=95)
#> [1] 0.6891277 0.9108723
```
