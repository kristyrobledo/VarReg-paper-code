
# VarReg paper code

This repository is for the data and code presented in the paper
published here (to be added).

## Install packages

Firstly, you will need to install the `VarReg` package and a few others
to perform the code below.

``` r
library(VarReg)
library(gamlss) #for the CD4 dataset
library(ggplot2)
library(plotly)
```

## CD4 dataset

This dataset is located within the `gamlss` package. CD4 is a type of
white blood cell, and in this dataset, it has been measured in
uninfected children born from HIV-1 infected women \[@Wade1994\]. The
dataset contains 609 measurements of CD4 cell counts and the child’s age
at which the measurements were taken.

``` r
data("CD4")

str(CD4)
```

    ## 'data.frame':    609 obs. of  2 variables:
    ##  $ cd4: num  387 2183 904 1681 656 ...
    ##  $ age: num  4.5 0.83 2.06 1.44 2.67 1.17 1.94 1.72 2.54 1.66 ...

Lets reproduce the graphic in the paper, showing that at younger ages
there is more variation in the CD4 counts than at older ages,
demonstrating heteroscedasticity.

``` r
p<-ggplot(data=CD4, aes(y=cd4, x=age)) +
  geom_point()+
  xlab("Age of child in years")+
  ylab("CD4 cell count")+
  theme_minimal()
ggplotly(p)
```

![](readme_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

## Viral load dataset
