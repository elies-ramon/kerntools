# Chi-squared kernel

\`Chi2()\` computes the basic \\\chi^2\\ kernel for bag-of-words (BoW)
or bag-of-visual-words data. This kernel computes the similarity between
two nonnegative vectors that represent the occurrence counts of words in
two different documents.

## Usage

``` r
Chi2(X, g = NULL)
```

## Arguments

- X:

  Matrix or data.frame (dimension *NxD*) that contains nonnegative
  numbers. Each row represents the counts of words of *N* documents,
  while each column is a word.

- g:

  Gamma hyperparameter. If g=0 or NULL, \`Chi2()\` returns the LeCam
  distances between the documents instead of the \\\chi^2\\ kernel
  matrix. (Defaults=NULL).

## Value

Kernel matrix (dimension: *NxN*).

## References

Zhang, Jianguo, et al. Local features and kernels for classification of
texture and object categories: A comprehensive study. International
journal of computer vision 73 (2007): 213-238.
[Link](https://inria.hal.science/inria-00548574/document)

## Examples

``` r
## Example dataset: word counts in 4 documents
documents <- matrix( c(0, 1, 3, 2, 1, 0,  1, 1, 6,4,3,1,3,5,6,2), nrow=4,byrow=TRUE)
rownames(documents) <- paste0("doc",1:4)
colnames(documents) <- c("animal","life","tree","ecosystem")
documents
#>      animal life tree ecosystem
#> doc1      0    1    3         2
#> doc2      1    0    1         1
#> doc3      6    4    3         1
#> doc4      3    5    6         2
Chi2(documents,g=NULL)
#>          1        2        3        4
#> 1 0.000000 1.290994 2.016598 1.825742
#> 2 1.290994 0.000000 2.070197 2.225395
#> 3 2.016598 2.070197 0.000000 1.105542
#> 4 1.825742 2.225395 1.105542 0.000000
```
