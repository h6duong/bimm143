Class7: R functions and packages
================
Han Duong
April 23, 2019

Functions revisted
==================

We will source a file from online with our function from last slide.

``` r
source("http://tinyurl.com/rescale-R")
```

Try out the last day's rescale() function

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

warning() and stop() functions are used inside functions to handle and report on unexpected situations

"Fail early and loudly"

Try rescale2() function that catches string inputs

``` r
rescale2(1:10, "string")
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

Find missing NA values in two vectors
=====================================

Now we're going to write a function both\_na(), that counts how many positions in two input vectors, x and y, both have a missing value 1. Always with a simple definition of the problem

``` r
# Lets define an example x and y
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)

#Use the is.na() function (can be found by simply searching "r find missing values in a vector")
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
#this function reports TRUE when it encounters NA in a vector, but how do you get it to report only NA present in both vectors
```

Try putting these together with an AND

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

Take the sum() to find out how many TRUE NAs we have and thus how many NAs we had in both x and y

``` r
#our working snippet
sum(is.na(x) & is.na(y))
```

    ## [1] 1

Now I can make this into our first function

``` r
both_na <- function(x, y){
  sum(is.na(x) & is.na(y))
}
```

``` r
both_na(x, c(NA, 3, NA, 2, NA))
```

    ## [1] 2

Test, test, test

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
```

This function recycles the first value of x to make it the same length as y2. Which gives x &lt;- c(NA, NA, NA, NA), resulting in "3". Different vector lengths might result in an inaccurate value.

``` r
both_na(x, y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

``` r
3==3
```

    ## [1] TRUE

Now let's try the bothna2() function on our different length input vectors

``` r
which(c(F, F, T, F, T))
```

    ## [1] 3 5

``` r
#which (is.na(c(1, 2, NA, 4)))
```

intersect function
------------------

``` r
df2
```

    ##     IDs exp
    ## 1 gene2  -2
    ## 2 gene4  NA
    ## 3 gene3   1
    ## 4 gene5   2

Make things simple, turn them into single vectors

``` r
x <- df1$IDs
y <- df2$IDs

x
```

    ## [1] "gene1" "gene2" "gene3"

``` r
y
```

    ## [1] "gene2" "gene4" "gene3" "gene5"

R intersection of two vectors

``` r
intersect(x, y)
```

    ## [1] "gene2" "gene3"

``` r
which (x %in% y)
```

    ## [1] 2 3

``` r
x[x %in% y]
```

    ## [1] "gene2" "gene3"

Use the RStudio shortcut 'CODE &gt; EXTRACT FUNCTION' to turn our snippet into a working function (highlighted)

``` r
gene_intersect <- function(x, y) {
  cbind(x[x %in% y],
        y[ y%in% x])
}
```

``` r
gene_intersect(df1$IDs, df2$IDs)
```

    ##      [,1]    [,2]   
    ## [1,] "gene2" "gene2"
    ## [2,] "gene3" "gene3"

``` r
gene_intersect2(df1, df2)
```

    ##     IDs exp df2[df2$IDs %in% df1$IDs, "exp"]
    ## 2 gene2   1                               -2
    ## 3 gene3   1                                1

``` r
gene_intersect3(df1, df2)
```

    ##     IDs exp exp2
    ## 2 gene2   1   -2
    ## 3 gene3   1    1

Create a function
=================

grade() that finds average grade dropping the worst homework score

``` r
x <- c(100,100,100,100,100,100,100,90)
y <- c(100,90,90,90,90,90,97,80)

grade <- function(x) {
  score <- (sum(x) - min(x))/ (length(x)-1)
if (score < 90){
  paste("You got",score,".", " Work harder!")
}
}
```

``` r
grade(c(70,50,40,80,90,100))
```

    ## [1] "You got 78 .  Work harder!"
