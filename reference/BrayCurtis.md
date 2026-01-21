# Kernels for count data

Ruzicka and Bray-Curtis are kernel functions for absolute or relative
frequencies and count data. Both kernels have as input a matrix or
data.frame with dimension *NxD* and *N*\>1, *D*\>1, containing strictly
non-negative real numbers. Samples should be in the rows. Thus, when
working with relative frequencies, \`rowSums(X)\` should be 1 (or 100,
or another arbitrary number) for *all* rows (samples) of the dataset.

## Usage

``` r
BrayCurtis(X)

Ruzicka(X)
```

## Arguments

- X:

  Matrix or data.frame that contains absolute or relative frequencies.

## Value

Kernel matrix (dimension: *NxN*).

## Details

For more info about these measures, please check Details in
?vegan::vegdist(). Note that, in the vegan help page, "Ruzicka"
corresponds to "quantitative Jaccard". \`BrayCurtis(X)\` gives the same
result than \`1-vegan::vegdist(X,method="bray")\`, and the same with
\`Ruzicka(data)\` and \`1-vegan::vegdist(data,method="jaccard")\`.

## Examples

``` r
data <- soil$abund
Kruz <- Ruzicka(data)
Kbray <- BrayCurtis(data)
Kruz[1:5,1:5]
#>            1          2          3          4          5
#> 1 1.00000000 0.12115505 0.11891559 0.03207739 0.12395094
#> 2 0.12115505 1.00000000 0.09132161 0.05955335 0.12264724
#> 3 0.11891559 0.09132161 1.00000000 0.04538870 0.07617411
#> 4 0.03207739 0.05955335 0.04538870 1.00000000 0.04447776
#> 5 0.12395094 0.12264724 0.07617411 0.04447776 1.00000000
Kbray[1:5,1:5]
#>            1         2          3          4          5
#> 1 1.00000000 0.2161254 0.21255507 0.06216083 0.22056289
#> 2 0.21612542 1.0000000 0.16735967 0.11241218 0.21849648
#> 3 0.21255507 0.1673597 1.00000000 0.08683603 0.14156466
#> 4 0.06216083 0.1124122 0.08683603 1.00000000 0.08516746
#> 5 0.22056289 0.2184965 0.14156466 0.08516746 1.00000000
```
