# Precision or PPV

\`Prec()\` computes the Precision of PPV (Positive Predictive Value)
between the output of a classification model and the actual values of
the target. The precision of each class can be aggregated.
Macro-precision is the average of the precision of each classes.
Micro-precision is the weighted average.

## Usage

``` r
Prec(ct, multi.class = "macro")
```

## Arguments

- ct:

  Confusion Matrix.

- multi.class:

  Should the results of each class be aggregated, and how? Options:
  "none", "macro", "micro". (Defaults: "macro").

## Value

PPV (a single value).

## Examples

``` r
y <- c(rep("a",3),rep("b",2))
y_pred <- c(rep("a",2),rep("b",3))
ct <- table(y,y_pred)
Prec(ct)
#> It is identical to weighted Accuracy
#> [1] 0.8333333
```
