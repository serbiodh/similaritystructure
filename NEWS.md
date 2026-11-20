# similaritystructure 0.2.0 (development)

* Rename internal functions pnalpha_scaled, fitting_pnalpha and fitting_pnalpha_norm to sn_norm, fitting_sn and fitting_sn_norm.
* Rename output fields to scaled_bn, Sn and param_Sn.
* Compute pS as the minimum of both cross-tail probabilities, using actual integer sample sizes and including the endpoint.
* Supply documented defaults and an approximately equidistant grid including both endpoints.
* Retain smoothing when the parametric fit fails, validate inputs and restore the previous future plan.
* Export adonis_ss for custom PERMANOVA tests.
* Document extraction of the similarity structure and comparison of effects.

This version changes numerical results for pS and the subsampling grid.
Historical cached analyses should be preserved and new article results
recomputed and checked before publication.
