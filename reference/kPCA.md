# Kernel PCA

\`kPCA()\` computes the kernel PCA from a kernel matrix and, if desired,
produces a plot. The contribution of the original variables to the
Principal Components (PCs), sometimes referred as "loadings", is NOT
returned (to do so, go to \`kPCA_imp()\`).

## Usage

``` r
kPCA(
  K,
  center = TRUE,
  Ktest = NULL,
  plot = NULL,
  y = NULL,
  colors = "black",
  na_col = "grey70",
  title = "Kernel PCA",
  pos_leg = "right",
  name_leg = "",
  labels = NULL,
  ellipse = NULL
)
```

## Arguments

- K:

  Kernel matrix (class "matrix").

- center:

  A logical value. If TRUE, the variables are zero-centered before the
  PCA. (Defaults: TRUE).

- Ktest:

  (optional) An additional kernel matrix corresponding to test samples,
  with dimension *Ntest x Ntraining*. These new samples are projected
  (using the color defined by \`na_col\`) over the kernel PCA computed
  from K. Remember than the data that generated \`Ktest\` should be
  centered beforehand, using the same values used for centering \`K\`.

- plot:

  (optional) A \`ggplot2\` is displayed. The input should be a vector of
  integers with length 2, corresponding to the two Principal Components
  to be displayed in the plot.

- y:

  (optional) A factor, or a numeric vector, with length equal to
  \`nrow(K)\` (number of samples). This parameter allows to paint the
  points with different colors.

- colors:

  A single color, or a vector of colors. If \`y\` is numeric, a gradient
  of colors between the first and the second entry will be used to paint
  the points. (Defaults: "black").

- na_col:

  Color of the entries that have a NA in the parameter \`y\`, or the
  entries corresponding to \`Ktest\` (when \`Ktest\` is not NULL).
  Otherwise, this parameter is ignored.

- title:

  Plot title.

- pos_leg:

  Position of the legend.

- name_leg:

  Title of the legend. (Defaults: blank)

- labels:

  (optional) A vector of the same length than nrow(K). A name will be
  displayed next to each point.

- ellipse:

  (optional) A float between 0 and 1. An ellipse will be drawn for each
  group of points defined by \`y\`. Here \`y\` should be of class
  "factor." This parameter will indicate the spread of the ellipse.

## Value

A list with two objects:

\* The PCA projection (class "matrix"). Please note that if K was
computed from a *NxD* table with *N \> D*, only the first *N-D* PCs may
be useful.

\* (optional) A \`ggplot2\` plot of the selected PCs.

## Details

As the ordinary PCA, kernel PCA can be used to summarize, visualize
and/or create new features of a dataset. Data can be projected in a
linear or nonlinear way, depending on the kernel used. When the kernel
is \`Linear()\`, kernel PCA is equivalent to ordinary PCA.

## Examples

``` r
dat <- matrix(rnorm(150),ncol=50,nrow=30)
K <- Linear(dat)

## Projection's coordinates only:
pca <- kPCA(K)

## Coordinates + plot of the two first principal components (PC1 and PC2):
pca <- kPCA(K,plot=1:2, colors = "coral2")
pca$plot
```
