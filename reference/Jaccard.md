# Kernels for sets

\`Intersect()\` or \`Jaccard()\` compute the kernel functions of the
same name, which are useful for set data. Their input is a matrix or
data.frame with dimension *NxD*, where *N*\>1, *D*\>0. Samples should be
in the rows and features in the columns. When there is a single feature,
\`Jaccard()\` returns 1 if the elements of the set are exactly the same
in two given samples, and 0 if they are completely different (see
Details). Instead, in the multivariate case (*D*\>1), the results (for
both \`Intersect()\` and \`Jaccard()\`) of the *D* features are combined
with a sum, a mean, or a weighted mean.

## Usage

``` r
Jaccard(X, elements = LETTERS, comp = "sum", coeff = NULL)

Intersect(
  X,
  elements = LETTERS,
  comp = "sum",
  coeff = NULL,
  feat_space = FALSE
)
```

## Arguments

- X:

  Matrix (class "character") or data.frame (class "character", or
  columns = "factor"). The elements in X are assumed to be categorical
  in nature.

- elements:

  All potential elements (symbols) that can appear in the sets. If there
  are some elements that are not of interest, they can be excluded so
  they are not taken into account by these kernels. (Defaults: LETTERS).

- comp:

  When *D*\>1, this argument indicates how the variables of the dataset
  are combined. Options are: "mean", "sum" and "weighted". (Defaults:
  "mean")

  - "sum" gives the same importance to all variables, and returns an
    unnormalized kernel matrix.

  - "mean" gives the same importance to all variables, and returns a
    normalized kernel matrix (all its elements range between 0 and 1).

  - "weighted" weights each variable according to the \`coeff\`
    parameter, and returns a normalized kernel matrix.

- coeff:

  (optional) A vector of weights with length *D*.

- feat_space:

  (not available for the Jaccard kernel). If FALSE, only the kernel
  matrix is returned. Otherwise, the feature space is returned too.
  (Defaults: FALSE).

## Value

Kernel matrix (dimension: *NxN*), or a list with the kernel matrix and
the feature space.

## Details

Let \\A,B\\ be two sets. Then, the Intersect kernel is defined as:

\$\$K\_{Intersect}(A,B)=\|A \cap B\| \$\$

And the Jaccard kernel is defined as:

\$\$K\_{Jaccard}(A,B)=\|A \cap B\| / \|A \cup B\|\$\$

This specific implementation of the Intersect and Jaccard kernels
expects that the set members (elements) are character symbols
(length=1). In case the set data is multivariate (*D*\>1 columns, and
each one contains a set feature), elements for the *D* sets should come
from the same domain (universe). For instance, a dataset with two
variables, so the elements in the first one are colors
c("green","black","white","red") and the second are names
c("Anna","Elsa","Maria") is not allowed. In that case, set factors
should be recoded to colors c("g","b","w","r") and names c("A","E","M")
and, if necessary, 'Intersect()' (or \`Jaccard()\`) should be called
twice.

## References

Bouchard, M., Jousselme, A. L., and Doré, P. E. (2013). A proof for the
positive definiteness of the Jaccard index matrix. International Journal
of Approximate Reasoning, 54(5), 615-626.

Ruiz, F., Angulo, C., and Agell, N. (2008). Intersection and
Signed-Intersection Kernels for Intervals. Frontiers in Artificial
Intelligence and Applications. 184. 262-270. doi:
10.3233/978-1-58603-925-7-262.

## Examples

``` r
# Sets data
## Generating a dataset with sets containing uppercase letters
random_set <- function(x)paste(sort(sample(LETTERS,x,FALSE)),sep="",collapse = "")
max_setsize <- 4
setsdata <- matrix(replicate(20,random_set(sample(2:max_setsize,1))),nrow=4,ncol=5)

## Computing the Intersect kernel:
Intersect(setsdata,elements=LETTERS,comp="sum")
#>    1  2  3  4
#> 1 13  2  1  2
#> 2  2 16  1  3
#> 3  1  1 15  0
#> 4  2  3  0 16

## Computing the Jaccard kernel weighting the variables:
coeffs <- c(0.1,0.15,0.15,0.4,0.20)
Jaccard(setsdata,elements=LETTERS,comp="weighted",coeff=coeffs)
#>            1          2     3          4
#> 1 1.00000000 0.06000000 0.025 0.05357143
#> 2 0.06000000 1.00000000 0.080 0.08666667
#> 3 0.02500000 0.08000000 1.000 0.00000000
#> 4 0.05357143 0.08666667 0.000 1.00000000
```
