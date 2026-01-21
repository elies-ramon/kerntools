# Compositional kernels

\`cLinear()\` is the compositional-linear kernel, which is useful for
compositional data (relative frequencies or proportions).
\`Aitchison()\` is akin to the RBF kernel for this type of data. Thus,
the expected input for both kernels is a matrix or data.frame containing
strictly non-negative or (even better) positive numbers. This input has
dimension *NxD*, with *N*\>1 samples and *D*\>1 compositional features.

## Usage

``` r
cLinear(X, cos.norm = FALSE, feat_space = FALSE, zeros = "none")

Aitchison(X, g = NULL, zeros = "none")
```

## Arguments

- X:

  Matrix or data.frame that contains the compositional data.

- cos.norm:

  Should the resulting kernel matrix be cosine normalized? (Defaults:
  FALSE).

- feat_space:

  If FALSE, only the kernel matrix is returned. Otherwise, the feature
  space is also returned. (Defaults: FALSE).

- zeros:

  "none" to warrant that there are no zeroes in X, "pseudo" to replace
  zeroes by a pseudocount. (Defaults="none").

- g:

  Gamma hyperparameter. If g=0 or NULL, the matrix of squared Aitchison
  distances is returned instead of the Aitchison kernel matrix.
  (Defaults=NULL).

## Value

Kernel matrix (dimension: *NxN*).

## Details

In compositional data, samples (rows) sum to an arbitrary or irrelevant
number. This is most clear when working with relative frequencies, as
all samples add to 1 (or 100, or other uninformative value). Zeroes are
a typical challenge when using compositional approaches. They introduce
ambiguity because they can have multiple causes; a zero may signal a
true absence, or a value so small that it is below the detection
threshold of an instrument. A simple approach to deal with zeroes is
replacing them by a pseudocount. More sophisticated approaches are
reviewed elsewhere; see for instance the R package \`zCompositions\`.

## References

Ramon, E., Belanche-Muñoz, L. et al (2021). kernInt: A kernel framework
for integrating supervised and unsupervised analyses in spatio-temporal
metagenomic datasets. Frontiers in microbiology 12 (2021): 609048. doi:
10.3389/fmicb.2021.609048

## Examples

``` r
data <- soil$abund

## This data is sparse and contains a lot of zeroes. We can replace them by pseudocounts:
Kclin <- cLinear(data,zeros="pseudo")
Kclin[1:5,1:5]
#>            X103.CA2  X103.CO3   X103.SR3   X103.IE2   X103.BP1
#> X103.CA2 10275.8219  1904.924  2674.1133   416.1493  2110.8775
#> X103.CO3  1904.9238 10322.174  2081.3274  1388.0667  2174.0920
#> X103.SR3  2674.1133  2081.327 11933.4557   688.7321  1751.3419
#> X103.IE2   416.1493  1388.067   688.7321 13005.5436   805.6387
#> X103.BP1  2110.8775  2174.092  1751.3419   805.6387 10272.4113

## With the feature space:
Kclin <- cLinear(data,zeros="pseudo",feat_space=TRUE)

## With cosine normalization:
Kcos <- cLinear(data,zeros="pseudo",cos.norm=TRUE)
Kcos[1:5,1:5]
#>            X103.CA2  X103.CO3   X103.SR3   X103.IE2   X103.BP1
#> X103.CA2 1.00000000 0.1849625 0.24148402 0.03599786 0.20545587
#> X103.CO3 0.18496253 1.0000000 0.18753041 0.11980101 0.21113301
#> X103.SR3 0.24148402 0.1875304 1.00000000 0.05528444 0.15818002
#> X103.IE2 0.03599786 0.1198010 0.05528444 1.00000000 0.06970113
#> X103.BP1 0.20545587 0.2111330 0.15818002 0.06970113 1.00000000

## Aitchison kernel:
Kait <- Aitchison(data,g=0.0001,zeros="pseudo")
Kait[1:5,1:5]
#>           X103.CA2  X103.CO3   X103.SR3   X103.IE2  X103.BP1
#> X103.CA2 1.0000000 0.1865950 0.18523961 0.10593743 0.1954115
#> X103.CO3 0.1865950 1.0000000 0.16376915 0.12807255 0.1969826
#> X103.SR3 0.1852396 0.1637692 1.00000000 0.09478411 0.1540746
#> X103.IE2 0.1059374 0.1280725 0.09478411 1.00000000 0.1145587
#> X103.BP1 0.1954115 0.1969826 0.15407461 0.11455872 1.0000000
```
