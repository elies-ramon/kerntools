# Aggregate importances

\`aggregate_imp()\` sums the importances present in a matrix or
data.frame according to some user-specified grouping criterion.

## Usage

``` r
aggregate_imp(X, lev = NULL, samples = "rows")
```

## Arguments

- X:

  Matrix or data.frame containing the importances (in rows or in
  columns).

- lev:

  (optional) The grouping elements. \`lev\` should be as long as the
  dimension (cols or rows) that one wants to aggregate. If this
  parameter is absent, the colnames (if samples="rows") or rownames will
  be used to that effect. In that case, it is expected that the
  col/rownames follow this pattern: "V_Y", and the variables with the
  same "V" will be summed. (Check the colnames of a typical output of
  \`dummy_data()\` for more info).

- samples:

  Samples are in rows or in columns? (Defaults: "rows").

## Value

X, a matrix or data.frame containing the aggregated importances.

## Examples

``` r
importances <- matrix(rnorm(90),nrow=3,ncol=30)
rownames(importances) <- c("sample1","sample2","sample3")
colnames(importances) <- paste0("Feat",
rep(1:5,times=2*(1:5)), "_", unlist(lapply(2*(1:5),function(x)LETTERS[1:x])))

## The grouping criterion is:
groups <- paste0("Feat",1:5)
aggregate_imp(X=importances,samples="rows",lev=groups)
#>              Feat1       Feat2       Feat3    Feat4     Feat5
#> sample1  0.5234091 -0.01500287 -0.02722235 1.663339  2.098130
#> sample2 -0.1697910 -2.83341389 -0.68369776 1.039239  5.246396
#> sample3 -1.1724617 -1.50500564  1.84112488 2.194019 -1.930851
## We can also use the colnames:
colnames(importances)
#>  [1] "Feat1_A" "Feat1_B" "Feat2_A" "Feat2_B" "Feat2_C" "Feat2_D" "Feat3_A"
#>  [8] "Feat3_B" "Feat3_C" "Feat3_D" "Feat3_E" "Feat3_F" "Feat4_A" "Feat4_B"
#> [15] "Feat4_C" "Feat4_D" "Feat4_E" "Feat4_F" "Feat4_G" "Feat4_H" "Feat5_A"
#> [22] "Feat5_B" "Feat5_C" "Feat5_D" "Feat5_E" "Feat5_F" "Feat5_G" "Feat5_H"
#> [29] "Feat5_I" "Feat5_J"
aggregate_imp(X=importances,samples="rows")
#>              Feat1       Feat2       Feat3    Feat4     Feat5
#> sample1  0.5234091 -0.01500287 -0.02722235 1.663339  2.098130
#> sample2 -0.1697910 -2.83341389 -0.68369776 1.039239  5.246396
#> sample3 -1.1724617 -1.50500564  1.84112488 2.194019 -1.930851
```
