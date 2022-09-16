
<!-- README.md is generated from README.Rmd. Please edit that file -->

# candiSNP

<!-- badges: start -->
<!-- badges: end -->

candiSNP is now an `R` package that you can run on your own computer.
You will need to install `R` and some packages to make it work. You will
also need a Java installation but that *should* already be present on
your machine.

## Installation

### Install R

The `R` statistical programming language can be installed from one of
the CRAN mirrors listed here <https://cran.r-project.org/mirrors.html>.
You will need version 4.1.0 or later.

### Install the `candiSNP` package

Once you have installed R, use the R console to install the development
version of candiSNP using the `devtools` package, install that first.

``` r
install.packages('devtools')
```

Then you can install `candiSNP`

``` r
devtools::install_github("TeamMacLean/candiSNP")
```

## Installing genome annotations for SNPEff

`candiSNP` brings with it version 3.6 of `SNPEff`, but not any genome
annotations. To install the default genomes used in the original web
version of `candiSNP` you can use the `install_default_genomes()`
function.

At the `R` console, type

``` r
library("candiSNP")
install_default_genomes()
```

This operation can take a few minutes, depending on your internet speed.
You should only need to do it once, though.

## Starting the app

Once the genomes are installed, you can proceed to start `candiSNP` as
follows.

``` r
app()
```

Note that the app starts in a browser window. The `R` console remains
busy while the app is running. Use `Ctrl-C` to quit the process and get
your console back.