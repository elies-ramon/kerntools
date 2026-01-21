# Specificity or TNR

\`Spe()\` computes the Specificity or TNR (True Negative Rate) between
the output of a classification prediction model and the actual values of
the target.

## Usage

``` r
Spe(ct, multi.class = "macro")
```

## Arguments

- ct:

  Confusion Matrix.

- multi.class:

  Should the results of each class be aggregated, and how? Options:
  "none", "macro", "micro". (Defaults: "macro").

## Value

TNR (a single value).

## Examples

``` r
y <- c(rep("a",3),rep("b",2))
y_pred <- c(rep("a",2),rep("b",3))
ct <- table(y,y_pred)
Spe(ct)
#> [1] 0.8333333
```
