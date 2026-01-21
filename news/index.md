# Changelog

## kerntools 1.2.0

CRAN release: 2025-02-19

### Major changes

- Two new kernels, ‘Aitchison()’ and ‘cLinear()’, for compositional
  data.
- A new kernel, ‘Chi2()’, for computing the chi-squared kernel for
  bags-of-words.
- New function, ‘aggregate_imp()’, that aggregates a table of features
  following a user-defined grouping criterion.

### Minor improvements and bug fixes

- ‘plotImp()’ now throws an error when the user requests more features
  than the actual number of features present in x. Also, problems when
  displaying the correct xlim are fixed.
- ‘heatK()’ allows the user to input a legend title.
- ‘KTA()’ now only accepts binary variables for its parameter ‘y’.

## kerntools 1.1.0

CRAN release: 2024-10-25

- New vignettes explaining in depth kernel PCA and the kernel functions
  implemented in this package.

- New example dataset (‘soil’) that contains bacterial counts.

## kerntools 1.0.2

CRAN release: 2024-09-04

- Fixed errors that arose during CRAN Package Check when using
  alternative BLAS/LAPACK implementations.

## kerntools 1.0.1

CRAN release: 2024-08-27

- Fixed some typos in the DESCRIPTION file and the ‘Acc()’ function.

## kerntools 1.0.0

- Initial CRAN submission.
