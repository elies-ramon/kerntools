# Kernel matrix heatmap

\`heatK()\` plots the heatmap of a kernel matrix.

## Usage

``` r
heatK(
  K,
  cos.norm = FALSE,
  title = NULL,
  color = c("red", "yellow"),
  name_leg = NULL,
  raster = FALSE
)
```

## Arguments

- K:

  Kernel matrix (class "matrix").

- cos.norm:

  If TRUE, the cosine normalization is applied to the kernel matrix so
  its elements have a maximum value of 1. (Defaults: FALSE).

- title:

  Heatmap title (optional).

- color:

  A vector of length 2 containing two colors. The first color will be
  used to represent the minimum value and the second the maximum value
  of the kernel matrix.

- name_leg:

  Title of the legend.

- raster:

  In large kernel matrices, raster = TRUE will draw quicker and
  better-looking heatmaps. (Defaults=FALSE).

## Value

A \`ggplot2\` heatmap.

## Examples

``` r
data <- matrix(rnorm(150),ncol=50,nrow=30)
K <- Linear(data)
heatK(K)
```
