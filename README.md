<!-- README.md is generated from README.Rmd. Please edit README.Rmd. -->

# similaritystructure

The **similaritystructure** package estimates the similarity structure (SS) from repeated subsampling and hypothesis tests. It returns a normalized distribution over subsample sizes, the similarity size (`nS`), and the similar effect probability (`pS`) for comparisons between effects.

Two approximations are available: parametric fitting and binomial GAM smoothing. The first is normalized over an infinite range; the second over the evaluated range of sample sizes.

## Installation

```r
remotes::install_github("serbiodh/similaritystructure")
```

The current development version is in the `update-paper` branch:

```r
remotes::install_github("serbiodh/similaritystructure", ref = "update-paper")
```

## A first similarity structure

```r
library(similaritystructure)

set.seed(12345)
A <- rnorm(600)
B <- rnorm(600, mean = 0.5)
C <- rnorm(600, mean = 1)

AB <- similarity_structure(
  A, B, N_fin = 100, num_N = 20, num_repet = 50,
  replace = TRUE, test = "t-test", alpha = 0.05,
  seed = 12345, plotting = FALSE
)

AC <- similarity_structure(
  A, C, N_fin = 100, num_N = 20, num_repet = 50,
  replace = TRUE, test = "t-test", alpha = 0.05,
  seed = 12346, plotting = FALSE
)
```

The empirical non-rejection curve is in `scaled_bn`. The normalized SS and similarity size are provided by `smooth` and, when a valid parametric fit is found, by `fitting`.

```r
AB$scaled_bn

AB$simil_struct_Ns$smooth$Ns_smooth
plot(AB$simil_struct_Ns$smooth$vect_N_smooth,
     AB$simil_struct_Ns$smooth$Sn,
     type = "l", xlab = "Subsample size (n)", ylab = "S(n)")

if (!is.null(AB$simil_struct_Ns$fitting))
  AB$simil_struct_Ns$fitting$Ns_fitting
```

## Compare effects with pS

Use `prob_effect()` to compare two normalized similarity structures.

```r
prob_effect(smooth_A = AB$simil_struct_Ns$smooth,
            smooth_B = AC$simil_struct_Ns$smooth)

if (!is.null(AB$simil_struct_Ns$fitting) &&
    !is.null(AC$simil_struct_Ns$fitting)) {
  prob_effect(fitting_A = AB$simil_struct_Ns$fitting,
              fitting_B = AC$simil_struct_Ns$fitting)
}
```

`pS` is a descriptive comparison probability, not a hypothesis-test p-value. Interpret comparisons made using the same statistical test, significance level and compatible sample-size supports. In particular, examine the empirical curve and fit before interpreting extrapolated tails.

For more information, see `?similarity_structure` and `?prob_effect`.

## Reproducible examples

The package includes three R Markdown analyses:

| Example | Data | Test |
| --- | --- | --- |
| `example_01_mean.Rmd` | Normal data with 24 different mean shifts | KS test (t-test option) |
| `example_02_bmi.Rmd` | BMI by smoking status, NHANES 2017–2018 | Welch's t-test |
| `example_03_betadivers.Rmd` | Gut microbiome beta diversity by sex and smoking | PERMANOVA |

Each example runs quickly by default using results supplied with the package. Set `calculate_ss: true` to repeat the full subsampling analysis. The latter can take a long time, especially for the microbiome.

Copy the examples to a **writable folder** before running them:

```r
examples <- c("example_01_mean.Rmd",
              "example_02_bmi.Rmd",
              "example_03_betadivers.Rmd")

dir.create("examples", showWarnings = FALSE)

file.copy(
  file.path(system.file("examples", package = "similaritystructure"), examples),
  file.path("examples", examples),
  overwrite = TRUE
)

rmarkdown::render("examples/example_01_mean.Rmd")
rmarkdown::render("examples/example_02_bmi.Rmd")
rmarkdown::render("examples/example_03_betadivers.Rmd")
```

To perform the full analysis, for example for BMI:

```r
rmarkdown::render("examples/example_02_bmi.Rmd",
                  params = list(calculate_ss = TRUE))
```

The normal-data example also accepts `test = "t-test"` in its R Markdown parameters. The BMI analysis downloads public NHANES data when needed; the microbiome distance matrix and metadata are supplied with the package.

For the examples, install `ggplot2`, `haven`, `knitr` and `rmarkdown`; `vegan` is also needed to repeat the microbiome PERMANOVA subsampling.

## Notes on interpretation

Subsampling with replacement can include the same observed participant more than once. It does not create new independent observations, even when the evaluated subsample size exceeds the size of the original data.

For GAM smoothing, the normalized SS is defined on the finite evaluated interval. For the parametric fit, normalization uses the fitted infinite tail. Consequently, their `nS` and `pS` values need not agree exactly, particularly if the empirical curve has not decayed by the selected upper sample size.
