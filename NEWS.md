# similaritystructure 0.2.0 (development)

* Rename internal functions pnalpha_scaled, fitting_pnalpha and fitting_pnalpha_norm to sn_norm, fitting_sn and fitting_sn_norm.
* Rename output fields to scaled_bn, Sn and param_Sn.
* Compute pS as the minimum of both cross-tail probabilities, using actual integer sample sizes and including the endpoint.
* Supply documented defaults and an approximately equidistant grid including both endpoints.
* Retain smoothing when the parametric fit fails, validate inputs and restore the previous future plan.
* Export adonis_ss for custom PERMANOVA tests.
* Replace Gaussian smoothing with an unconstrained binomial GAM (logit link, REML), using original Monte Carlo counts. Retain the number of tests as an attribute of the empirical curve.
* Document extraction of the similarity structure and comparison of effects.

This version changes numerical results for pS, smoothing and the subsampling grid.
* Restrict main examples to normal mean differences, BMI and gut microbiome. Remove the skewness and eye-tracking examples and associated data. Retain standard-deviation material separately for the supplement.
* Add a cached-result updater that recalculates binomial smoothing from original counts without rerunning subsampling or fitting. Add a provisional NHANES BMI reconstruction pending comparison with the original analysis.

Historical cached analyses should be preserved and new article results
recomputed and checked before publication.

* Remove unused Imports and the vignette builder declaration (example Rmd files are not package vignettes); allow roxygen2 to regenerate prob_effect help.

* Replace the provisional NHANES BMI example with a verified, fully commented R Markdown analysis. Store original empirical non-rejection counts for fast refitting without MATLAB or repeated subsampling.

* Correct sample-ID punctuation in the gut microbiome example (1,011 samples). Replace four outdated cached results with original empirical counts and parametric coefficients from the corrected, with-replacement analysis; update binomial GAM smoothing and pS figures.
