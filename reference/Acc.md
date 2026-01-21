# Accuracy

\`Acc()\` computes the accuracy between the output of a classification
model and the actual values of the target. It can also compute the
weighted accuracy, which is useful in imbalanced classification
problems. The weighting is applied according to the class frequencies in
the target. In balanced problems, weighted Acc = Acc.

## Usage

``` r
Acc(ct, weighted = FALSE)
```

## Arguments

- ct:

  Confusion Matrix.

- weighted:

  If TRUE, the weighted accuracy is returned. (Defaults: FALSE).

## Value

Accuracy of the model (a single value).

## Examples

``` r
y <- c(rep("a",3),rep("b",2))
y_pred <- c(rep("a",2),rep("b",3))
ct <- table(y,y_pred)
Acc(ct)
#> [1] 0.8
Acc(ct,weighted=TRUE)
#> [1] 0.8333333
```
