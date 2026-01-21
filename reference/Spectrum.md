# Spectrum kernel

\`Spectrum()\` computes the basic Spectrum kernel between strings. This
kernel computes the similarity of two strings by counting how many
matching substrings of length *l* are present in each one.

## Usage

``` r
Spectrum(
  x,
  alphabet,
  l = 1,
  group.ids = NULL,
  weights = NULL,
  feat_space = FALSE,
  cos.norm = FALSE
)
```

## Arguments

- x:

  Vector of strings (length *N*).

- alphabet:

  Alphabet of reference.

- l:

  Length of the substrings.

- group.ids:

  (optional) A vector with ids. It allows to compute the kernel over
  groups of strings within x, instead of the individual strings.

- weights:

  (optional) A numeric vector as long as x. It allows to weight
  differently each one of the strings.

- feat_space:

  If FALSE, only the kernel matrix is returned. Otherwise, the feature
  space (i.e. a table with the number of times that a substring of
  length *l* appears in each string) is also returned (Defaults: FALSE).

- cos.norm:

  Should the resulting kernel matrix be cosine normalized? (Defaults:
  FALSE).

## Value

Kernel matrix (dimension: *NxN*), or a list with the kernel matrix and
the feature space.

## Details

In large datasets this function may be slow. In that case, you may use
the \`stringdot()\` function of the \`kernlab\` package, or the
\`spectrumKernel()\` function of the \`kebabs\` package.

## References

Leslie, C., Eskin, E., and Noble, W.S. The spectrum kernel: a string
kernel for SVM protein classification. Pac Symp Biocomput. 2002:564-75.
PMID: 11928508.
[Link](http://psb.stanford.edu/psb-online/proceedings/psb02/abstracts/p564.md)

## Examples

``` r
## Examples of alphabets. _ stands for a blank space, a gap, or the
## start or the end of sequence)
NT <- c("A","C","G","T","_") # DNA nucleotides
AA <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T",
"V","W","Y","_") ##canonical aminoacids
letters_ <- c(letters,"_")
## Example of data
strings <- c("hello_world","hello_word","hola_mon","kaixo_mundua",
"saluton_mondo","ola_mundo", "bonjour_le_monde")
names(strings) <- c("english1","english_typo","catalan","basque",
"esperanto","galician","french")
## Computing the kernel:
Spectrum(strings,alphabet=letters_,l=2)
#>              english1 english_typo catalan basque esperanto galician french
#> english1           10            8       0      1         0        0      0
#> english_typo        8            9       0      1         0        0      0
#> catalan             0            0       7      1         4        4      4
#> basque              1            1       1     11         2        4      2
#> esperanto           0            0       4      2        14        3      7
#> galician            0            0       4      4         3        8      2
#> french              0            0       4      2         7        2     17
```
