
<!-- README.md is generated from README.Rmd. Please edit that file -->

# similaritystructure

<!-- badges: start -->

<!-- badges: end -->

`similaritystructure` estimates similarity structure between two samples
using repeated subsampling, statistical testing, smoothing and
parametric fitting.

## Installation

You can install the development version of similaritystructure like so:

``` r
devtools::install_github("serbiodh/similaritystructure")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(similaritystructure)

seed <- 12345
set.seed(seed)
n1 <- rnorm(20000, 0, 1)
n2 <- n1 + 0.8

res <- similarity_structure(
  n1 = n1,
  n2 = n2,
  N_init = 2,
  N_fin = round((80/0.8^2)),
  num_N = 60,
  num_repet = 300,
  test = "t-test",
  alpha = 0.05,
  seed <- seed,
  plotting = TRUE
)
#> [1] 2
#> [1] 4
#> [1] 6
#> [1] 8
#> [1] 10
#> [1] 12
#> [1] 14
#> [1] 17
#> [1] 19
#> [1] 21
#> [1] 23
#> [1] 25
#> [1] 27
#> [1] 29
#> [1] 31
#> [1] 33
#> [1] 35
#> [1] 37
#> [1] 40
#> [1] 42
#> [1] 44
#> [1] 46
#> [1] 48
#> [1] 50
#> [1] 52
#> [1] 54
#> [1] 56
#> [1] 58
#> [1] 60
#> [1] 62
#> [1] 64
#> [1] 67
#> [1] 69
#> [1] 71
#> [1] 73
#> [1] 75
#> [1] 77
#> [1] 79
#> [1] 81
#> [1] 83
#> [1] 85
#> [1] 87
#> [1] 90
#> [1] 92
#> [1] 94
#> [1] 96
#> [1] 98
#> [1] 100
#> [1] 102
#> [1] 104
#> [1] 106
#> [1] 108
#> [1] 110
#> [1] 112
#> [1] 115
#> [1] 117
#> [1] 119
#> [1] 121
#> [1] 123
#> [1] 125
```

<img src="man/figures/README-example-1.png" width="100%" />
