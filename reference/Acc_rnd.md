# Accuracy of a random model

\`Acc_rnd()\` computes the expected accuracy of a random classifier
based on the class frequencies of the target. This measure can be used
as a benchmark when contrasted to the accuracy (in test) of a given
prediction model.

## Usage

``` r
Acc_rnd(target, freq = FALSE)
```

## Arguments

- target:

  A character vector or a factor. Alternatively, a numeric vector (see
  below).

- freq:

  TRUE if \`target\` contains the frequencies of the classes (in this
  case, \`target\` should be numeric), FALSE otherwise. (Defaults:
  FALSE).

## Value

Expected accuracy of a random classification model (a single value).

## Examples

``` r
# Expected accuracy of a random model:
target <- c(rep("a",5),rep("b",2))
Acc_rnd(target)
#> [1] 0.5918367
# This is the same than:
freqs <- c(5/7,2/7)
Acc_rnd(freqs,freq=TRUE)
#> [1] 0.5918367
```
