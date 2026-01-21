# Levels per factor variable

This function gives the categories ("levels") per categorical variable
("factor").

## Usage

``` r
dummy_var(X)
```

## Arguments

- X:

  A matrix, or a data.frame containing factors. (If the columns are of
  any other class, they will be coerced into factors anyway).

## Value

A list with the levels.

## Examples

``` r
summary(showdata)
#>    Favorite.color             Favorite.actress        Favorite.actor
#>  blue     :20     Sophie Turner       :22      Peter Dinklage:20    
#>  red      :20     Emilia Clarke       :19      Kit Harington :17    
#>  black    :14     Anya Chalotra       :12      Henry Cavill  :15    
#>  purple   :12     Freya Allan         :10      Lee Jung-jae  :15    
#>  green    : 8     Helena Bonham Carter: 8      Hyun Bin      : 8    
#>  lightblue: 8     Lorraine Ashbourne  : 7      Josh O'Connor : 7    
#>  (Other)  :18     (Other)             :22      (Other)       :18    
#>           Favorite.show Liked.new.show
#>  Game of Thrones :22    No :48        
#>  The witcher     :17    Yes:52        
#>  Bridgerton      :14                  
#>  Squid game      :11                  
#>  The crown       :11                  
#>  La casa de papel: 8                  
#>  (Other)         :17                  
dummy_var(showdata)
#> $Favorite.color
#>  [1] "black"     "blue"      "green"     "grey"      "lightblue" "orange"   
#>  [7] "purple"    "red"       "white"     "yellow"   
#> 
#> $Favorite.actress
#>  [1] "Alba Flores"          "Anya Chalotra"        "Diana Rigg"          
#>  [4] "Elisabet Casanovas"   "Emilia Clarke"        "Freya Allan"         
#>  [7] "Helena Bonham Carter" "Lorraine Ashbourne"   "Mj Rodriguez"        
#> [10] "Soo Ye-jin"           "Sophie Turner"       
#> 
#> $Favorite.actor
#>  [1] "Alvaro Morte"       "David Solans"       "Francesc Orella"   
#>  [4] "Henry Cavill"       "Hyun Bin"           "Jonathan Bailey"   
#>  [7] "Josh O'Connor"      "Kit Harington"      "Lee Jung-jae"      
#> [10] "Michael K Williams" "Peter Dinklage"    
#> 
#> $Favorite.show
#> [1] "Bridgerton"       "Game of Thrones"  "La casa de papel" "Merli"           
#> [5] "Pose"             "Squid game"       "The crown"        "The wire"        
#> [9] "The witcher"     
#> 
#> $Liked.new.show
#> [1] "No"  "Yes"
#> 
```
