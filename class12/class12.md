Class 12: NMA and View
================
Han Duong
June 5, 2019

``` r
#install.packages("devtools")
#devtools::install_bitbucket("Grantlab/bio3d-view")
#install.packages("rgl")
```

Normal Mode Analysis
====================

A bioinfomatics method to predict the intrinsic dynamics of biomolecules

``` r
library(bio3d)
library(bio3d.view)
library(rgl)
pdb <- read.pdb("1hel")
```

    ##   Note: Accessing on-line PDB file

``` r
modes <- nma(pdb)
```

    ##  Building Hessian...     Done in 0.03 seconds.
    ##  Diagonalizing Hessian...    Done in 0.19 seconds.

``` r
m7 <- mktrj(modes, mode=7, file="mode_7.pdb")

view(m7, col=vec2color(rmsf(m7)))
```

    ## Potential all C-alpha atom structure(s) detected: Using calpha.connectivity()

``` r
rglwidget(width=500, height=500)
```

    ## PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.

<!--html_preserve-->

<!--/html_preserve-->
The created mode\_7.pdb file can be viewed on VMD. use "TUBE" and hit play

``` r
pdb <- read.pdb("5p21")
```

    ##   Note: Accessing on-line PDB file

``` r
view(pdb)
```

    ## Computing connectivity from coordinates...

``` r
view(pdb, "overview", col= "sse")
```

    ## Computing connectivity from coordinates...

``` r
view(m7)
```

    ## Potential all C-alpha atom structure(s) detected: Using calpha.connectivity()