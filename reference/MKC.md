# Multiple Kernel (Matrices) Combination

Combination of kernel matrices coming from different datasets / feature
types into a single kernel matrix.

## Usage

``` r
MKC(K, coeff = NULL)
```

## Arguments

- K:

  A three-dimensional *NxDxM* array containing *M* kernel matrices.

- coeff:

  A vector of length *M* with the weight of each kernel matrix. If NULL,
  all kernel matrices have the same weight. (Defaults: NULL)

## Value

A kernel matrix.

## Examples

``` r
# For illustrating a possible use of this function, we work with a dataset
# that contains numeric and categorical features.

summary(mtcars)
#>       mpg             cyl             disp             hp       
#>  Min.   :10.40   Min.   :4.000   Min.   : 71.1   Min.   : 52.0  
#>  1st Qu.:15.43   1st Qu.:4.000   1st Qu.:120.8   1st Qu.: 96.5  
#>  Median :19.20   Median :6.000   Median :196.3   Median :123.0  
#>  Mean   :20.09   Mean   :6.188   Mean   :230.7   Mean   :146.7  
#>  3rd Qu.:22.80   3rd Qu.:8.000   3rd Qu.:326.0   3rd Qu.:180.0  
#>  Max.   :33.90   Max.   :8.000   Max.   :472.0   Max.   :335.0  
#>       drat             wt             qsec             vs        
#>  Min.   :2.760   Min.   :1.513   Min.   :14.50   Min.   :0.0000  
#>  1st Qu.:3.080   1st Qu.:2.581   1st Qu.:16.89   1st Qu.:0.0000  
#>  Median :3.695   Median :3.325   Median :17.71   Median :0.0000  
#>  Mean   :3.597   Mean   :3.217   Mean   :17.85   Mean   :0.4375  
#>  3rd Qu.:3.920   3rd Qu.:3.610   3rd Qu.:18.90   3rd Qu.:1.0000  
#>  Max.   :4.930   Max.   :5.424   Max.   :22.90   Max.   :1.0000  
#>        am              gear            carb      
#>  Min.   :0.0000   Min.   :3.000   Min.   :1.000  
#>  1st Qu.:0.0000   1st Qu.:3.000   1st Qu.:2.000  
#>  Median :0.0000   Median :4.000   Median :2.000  
#>  Mean   :0.4062   Mean   :3.688   Mean   :2.812  
#>  3rd Qu.:1.0000   3rd Qu.:4.000   3rd Qu.:4.000  
#>  Max.   :1.0000   Max.   :5.000   Max.   :8.000  
cat_feat_idx <- which(colnames(mtcars) %in% c("vs", "am"))

# vs and am are categorical variables. We make a list, with the numeric features
# in the first element and the categorical features in the second:
DATA <- list(num=mtcars[,-cat_feat_idx], cat=mtcars[,cat_feat_idx])
# Our N, D and M dimensions are:
N <- nrow(mtcars); D <- ncol(mtcars); M <- length(DATA)

# Now we prepare a kernel matrix:
K <- array(dim=c(N,N,M))
K[,,1] <- Linear(DATA[[1]],cos.norm = TRUE) ## Kernel for numeric data
K[,,2] <- Dirac(DATA[[2]]) ## Kernel for categorical data

# Here, K1 has the same weight than K2 when computing the final kernel, although
# K1 has 9 variables and K2 has only 2.
Kconsensus <- MKC(K)
Kconsensus[1:5,1:5]
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 1.0000000 0.9999976 0.7459507 0.4898233 0.7429409
#> [2,] 0.9999976 1.0000000 0.7460135 0.4898009 0.7428780
#> [3,] 0.7459507 0.7460135 1.0000000 0.7244304 0.4786498
#> [4,] 0.4898233 0.4898009 0.7244304 1.0000000 0.7489934
#> [5,] 0.7429409 0.7428780 0.4786498 0.7489934 1.0000000

# If we want to weight equally each one of the 11 variables in the final
# kernel, K1 will weight 9/11 and K2 2/11.
coeff <- sapply(DATA,ncol)
coeff
#> num cat 
#>   9   2 
Kweighted <- MKC(K,coeff=coeff)
Kweighted[1:5,1:5]
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 1.0000000 0.9999960 0.9024648 0.8015291 0.8975396
#> [2,] 0.9999960 1.0000000 0.9025676 0.8014924 0.8974368
#> [3,] 0.9024648 0.9025676 1.0000000 0.8672498 0.7832451
#> [4,] 0.8015291 0.8014924 0.8672498 1.0000000 0.9074437
#> [5,] 0.8975396 0.8974368 0.7832451 0.9074437 1.0000000
```
