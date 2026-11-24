<!-- README.md is generated from README.Rmd. Please edit that file -->

# similaritystructure

Estimate the similarity structure (SS) and similarity size of two samples
using repeated subsampling, hypothesis testing, smoothing and parametric
fitting. Compare the sizes of effects through the similar effect probability
(pS).

## Installation

```r
remotes::install_github("serbiodh/similaritystructure")
```

## Calculate and extract the similarity structure

```r
library(similaritystructure)

set.seed(12345)
A <- rnorm(1000)
B <- rnorm(1000, 0.5)
C <- rnorm(1000, 0.9)

AB <- similarity_structure(A, B, N_fin = 300,
                           test = "t-test", seed = 12345,
                           plotting = FALSE)
AC <- similarity_structure(A, C, N_fin = 300,
                           test = "t-test", seed = 12346,
                           plotting = FALSE)

# Empirical scaled non-rejection curve
AB$scaled_bn

# Similarity structure and similarity size: smoothing
AB$simil_struct_Ns$smooth$Sn
AB$simil_struct_Ns$smooth$Ns_smooth

plot(AB$simil_struct_Ns$smooth$vect_N_smooth,
     AB$simil_struct_Ns$smooth$Sn,
     type = "l", xlab = "Subsample size", ylab = "S(n)")

# Similarity structure and similarity size: parametric fitting
if (!is.null(AB$simil_struct_Ns$fitting)) {
  AB$simil_struct_Ns$fitting$Sn
  AB$simil_struct_Ns$fitting$Ns_fitting
  plot(AB$simil_struct_Ns$fitting$vect_N_fit,
       AB$simil_struct_Ns$fitting$Sn,
       type = "l", xlab = "Subsample size", ylab = "S(n)")
}
```

The empirical non-rejection curve is scaled by its first value for
parametric fitting. Smoothing uses the original counts in a binomial GAM
with a logit link and an unconstrained thin plate spline (REML), then
rescales the predictions for graphical comparison. Neither scaled curve
is the similarity structure itself.
The smoothed SS is normalized over the finite interval; the fitted SS
is normalized over the infinite support starting at N_init. The returned
fitted vector only covers a finite display range and need not sum to one.
Choose N_fin to cover the decay of the empirical curve; inspect the fit
before interpreting its extrapolated tail.

## Compare effect sizes using pS

```r
# Smoothing
prob_effect(smooth_A = AB$simil_struct_Ns$smooth,
            smooth_B = AC$simil_struct_Ns$smooth)

# Fitting, or both methods when the fits are available
if (!is.null(AB$simil_struct_Ns$fitting) &&
    !is.null(AC$simil_struct_Ns$fitting)) {
  prob_effect(fitting_A = AB$simil_struct_Ns$fitting,
              fitting_B = AC$simil_struct_Ns$fitting)

  prob_effect(fitting_A = AB$simil_struct_Ns$fitting,
              fitting_B = AC$simil_struct_Ns$fitting,
              smooth_A = AB$simil_struct_Ns$smooth,
              smooth_B = AC$simil_struct_Ns$smooth)
}
```

pS is the minimum of the two cross-tail probabilities evaluated at the
similarity size of the other comparison. It is not a hypothesis-test
p-value; a threshold such as 0.05 is a user-defined decision criterion.
Compare structures computed with the same test, significance level and
compatible sample-size supports.

For more information, run `?similarity_structure` and `?prob_effect`.

## Example files

The main-text examples are:

- `example_01_mean.Rmd`: normal samples with different means (Figure 3).
- `example_02_bmi.Rmd`: NHANES BMI and smoking status (Figure 4; original parameters and verified cohort).
- `example_03_betadivers.Rmd`: gut microbiome composition (Figure 5).

The BMI example includes 180 original non-rejection counts in `inst/extdata/bmi_counts.csv`.
It refits the binomial GAM and parametric model without repeating subsampling;
set `calculate_ss: true` in the Rmd to run missing comparisons and save their RDS files. Delete any old local RDS files to force a rerun.
Copy the Rmd to a writable working folder before knitting so NHANES XPT files
can be downloaded into `data/` and results saved into `results/`.

Historical cached results for the normal and microbiome examples retain
original Monte Carlo counts. They are updated to the binomial GAM on loading,
without repeating tests or refitting the parametric models. To reproduce the
current subsample-size grid, repeat the underlying subsampling.
The standard-deviation example is retained separately under `supplementary/`
for the supporting normal-distribution analyses (Supplementary Figure S5.B).

```r
system.file("extdata", package = "similaritystructure")
system.file("examples", package = "similaritystructure")

file <- system.file("examples", "example_01_mean.Rmd",
                    package = "similaritystructure")
file.edit(file)
```

The examples may require additional packages (tictoc, data.table, dplyr,
tidyr, pbapply, ggplot2, stringr, patchwork, viridis, vegan and haven).
Install the relevant packages before running the examples.

To update your own saved results without repeating subsampling:

```r
source(system.file("scripts", "upgrade_cached_results.R",
                   package = "similaritystructure"))
# Supply the VERIFIED number of tests per subsample size from the original run.
upgrade_cached_file("original.rds", "updated_binomial.rds", num_repet = 300)
```

The updater preserves the original empirical grid and parametric fitting.
It writes a separate file and never overwrites the original results.
