# Convert categorical data to dummies.

Given a matrix or data.frame containing character/factors, this function
performs one-hot-encoding.

## Usage

``` r
dummy_data(X, lev = NULL)
```

## Arguments

- X:

  A matrix, or a data.frame containing factors. (If the columns are of
  any other class, they will be coerced into factors anyway).

- lev:

  (optional) A vector with the categories ("levels") of each factor.

## Value

X (class: "matrix") after performing one-hot-encoding.

## Examples

``` r
summary(CO2)
#>      Plant             Type         Treatment       conc          uptake     
#>  Qn1    : 7   Quebec     :42   nonchilled:42   Min.   :  95   Min.   : 7.70  
#>  Qn2    : 7   Mississippi:42   chilled   :42   1st Qu.: 175   1st Qu.:17.90  
#>  Qn3    : 7                                    Median : 350   Median :28.30  
#>  Qc1    : 7                                    Mean   : 435   Mean   :27.21  
#>  Qc3    : 7                                    3rd Qu.: 675   3rd Qu.:37.12  
#>  Qc2    : 7                                    Max.   :1000   Max.   :45.50  
#>  (Other):42                                                                  
CO2_dummy <- dummy_data(CO2[,1:3],lev=dummy_var(CO2[,1:3]))
CO2_dummy[1:10,1:5]
#>    Plant_Qn1 Plant_Qn2 Plant_Qn3 Plant_Qc1 Plant_Qc3
#> 1          1         0         0         0         0
#> 2          1         0         0         0         0
#> 3          1         0         0         0         0
#> 4          1         0         0         0         0
#> 5          1         0         0         0         0
#> 6          1         0         0         0         0
#> 7          1         0         0         0         0
#> 8          0         1         0         0         0
#> 9          0         1         0         0         0
#> 10         0         1         0         0         0
```
