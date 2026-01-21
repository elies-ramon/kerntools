# Kendall's tau kernel

\`Kendall()\` computes the Kendall's tau, which happens to be a kernel
function for ordinal variables, ranks or permutations.

## Usage

``` r
Kendall(X, NA.as.0 = TRUE, samples.in.rows = FALSE, comp = "mean")
```

## Arguments

- X:

  When evaluating a single ordinal feature, X should be a numeric matrix
  or data.frame. If data is multivariate, X should be a list, and each
  ordinal/ranking feature should be placed in a different element of the
  list (see Examples).

- NA.as.0:

  Should NAs be converted to 0s? (Defaults: TRUE).

- samples.in.rows:

  If TRUE, the samples are considered to be in the rows. Otherwise, it
  is assumed that they are in the columns. (Defaults: FALSE).

- comp:

  If X is a list, this argument indicates how the ordinal/ranking
  variables are combined. Options are: "mean" and "sum". (Defaults:
  "mean").

## Value

Kernel matrix (dimension: *NxN*).

## References

Jiao, Y. and Vert, J.P. The Kendall and Mallows kernels for
permutations. International Conference on Machine Learning. PMLR, 2015.
[Link](https://proceedings.mlr.press/v37/jiao15.html)

## Examples

``` r
# 3 people are given a list of 10 colors. They rank them from most (1) to least
# (10) favorite
color_list <-  c("black","blue","green","grey","lightblue","orange","purple",
"red","white","yellow")
survey1 <- 1:10
survey2 <- 10:1
survey3 <- sample(10)
color <- cbind(survey1,survey2,survey3) # Samples in columns
rownames(color) <- color_list
Kendall(color)
#>             survey1     survey2     survey3
#> survey1  1.00000000 -1.00000000 -0.02222222
#> survey2 -1.00000000  1.00000000  0.02222222
#> survey3 -0.02222222  0.02222222  1.00000000

# The same 3 people are asked the number of times they ate 5 different kinds of
# food during the last month:
food <- matrix(c(10, 1,18, 25,30, 7, 5,20, 5, 12, 7,20, 20, 3,22),ncol=5,nrow=3)
rownames(food) <- colnames(color)
colnames(food) <- c("spinach", "chicken", "beef" , "salad","lentils")
# (we can observe that, for person 2, vegetables << meat, while for person 3
# is the other way around)
Kendall(food,samples.in.rows=TRUE)
#>         survey1 survey2 survey3
#> survey1     1.0     0.2     0.4
#> survey2     0.2     1.0    -0.4
#> survey3     0.4    -0.4     1.0

# We can combine this results:
dataset <- list(color=color,food=t(food)) #All samples in columns
Kendall(dataset)
#> Composition: Mean
#>            [,1]       [,2]       [,3]
#> [1,]  1.0000000 -0.4000000  0.1888889
#> [2,] -0.4000000  1.0000000 -0.1888889
#> [3,]  0.1888889 -0.1888889  1.0000000
```
