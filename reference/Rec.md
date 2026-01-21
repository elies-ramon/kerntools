# Recall or Sensitivity or TPR

\`Rec()\` computes the Recall, also known as Sensitivity or TPR (True
Positive Rate), between the output of a classification model and the
actual values of the target.

## Usage

``` r
Rec(ct, multi.class = "macro")
```

## Arguments

- ct:

  Confusion Matrix.

- multi.class:

  Should the results of each class be aggregated, and how? Options:
  "none", "macro", "micro". (Defaults: "macro").

## Value

TPR (a single value).

## Examples

``` r
y <- c(rep("a",3),rep("b",2))
y_pred <- c(rep("a",2),rep("b",3))
ct <- table(y,y_pred)
Rec(ct)
#> [1] 0.8333333
```
