# Contributions of the variables to the Principal Components ("loadings")

\`kPCA_imp()\` performs a PCA and a kernel PCA simultaneously and
returns the contributions of the variables to the Principal Components
(sometimes, these contributions are called "loadings") in Feature Space.
Optionally, it can also return the samples' projection (cropped to the
relevant PCs) and the values used to centering the variables in Feature
Space. It does not return any plot, nor it projects test data. To do so,
please use \`kPCA()\`.

## Usage

``` r
kPCA_imp(DATA, center = TRUE, projected = NULL, secure = FALSE)
```

## Arguments

- DATA:

  A matrix or data.frame (NOT a kernel matrix) containing the data in
  feature space. Please note that nrow(DATA) should be higher than
  ncol(DATA). If the Linear kernel is used, this feature space is simply
  the original space.

- center:

  A logical value. If TRUE, the variables are zero-centered. (Defaults:
  TRUE).

- projected:

  (optional) If desired, the PCA projection (generated, for example, by
  \`kPCA()\`) can be included. If DATA is big (especially in the number
  of rows) this may save some computation time.

- secure:

  (optional) If TRUE, it tests the quality of the loadings This may be
  slow. (Defaults: FALSE).

## Value

A list with three objects:

\* The PCA projection (class "matrix") using only the relevant Principal
Components.

\* The loadings.

\* The values used to center each variable in Feature Space.

## Details

This function may be not valid for all kernels. Do not use it with the
RBF, Laplacian, Bray-Curtis, Jaccard/Ruzicka, or Kendall's tau kernels
unless you know exactly what you are doing.

## Examples

``` r
dat <- matrix(rnorm(150),ncol=30,nrow=50)
contributions <- kPCA_imp(dat)
#> Do not use this function if the PCA was created with the RBF,
#>           Laplacian, Bray-Curtis, Jaccard/Ruzicka, or Kendall's tau kernels
contributions$loadings[c("PC1","PC2"),1:5]
#>           [,1]       [,2]      [,3]       [,4]       [,5]
#> PC1  0.2212261 -0.1874640 0.1261597  0.2212261 -0.1874640
#> PC2 -0.2257706 -0.1761048 0.1342190 -0.2257706 -0.1761048
```
