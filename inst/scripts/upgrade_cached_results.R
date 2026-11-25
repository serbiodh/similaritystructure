# Migrate stored names and optionally refit ONLY the smoothing component.
# num_repet must be verified from the original script or recorded attribute.
# No resampling or parametric refitting is performed.
upgrade_ss <- function(x, num_repet = NULL) {
  if (is.null(x$scaled_bn)) {
    x$scaled_bn <- if (!is.null(x$scaled_gn)) x$scaled_gn else x$empiric_PNalpha
  }
  x$scaled_gn <- NULL
  x$empiric_PNalpha <- NULL
  if (!is.null(x$simil_struct_Ns$smooth)) {
    if (is.null(x$simil_struct_Ns$smooth$Sn))
      x$simil_struct_Ns$smooth$Sn <- x$simil_struct_Ns$smooth$fun_smooth_norm
    x$simil_struct_Ns$smooth$fun_smooth_norm <- NULL
  }
  if (!is.null(x$simil_struct_Ns$fitting)) {
    if (is.null(x$simil_struct_Ns$fitting$Sn))
      x$simil_struct_Ns$fitting$Sn <- x$simil_struct_Ns$fitting$P_N_alpha
    if (is.null(x$simil_struct_Ns$fitting$param_Sn))
      x$simil_struct_Ns$fitting$param_Sn <- x$simil_struct_Ns$fitting$param_P_N_alpha
    x$simil_struct_Ns$fitting$P_N_alpha <- NULL
    x$simil_struct_Ns$fitting$param_P_N_alpha <- NULL
  }
  if (is.null(x$scaled_bn)) stop("No empirical curve found in saved result.")
  if (is.null(num_repet)) num_repet <- attr(x$scaled_bn, "num_repet")
  if (is.null(num_repet)) {
    stop("Supply the original number of tests per subsample size (num_repet).")
  }
  recorded <- attr(x$scaled_bn, "num_repet")
  if (!is.null(recorded) && recorded != num_repet) {
    stop("num_repet disagrees with the recorded number of tests.")
  }

  # Preserve the original smoothing interval, including any original
  # endpoint not represented in the older empirical grid.
  n <- x$simil_struct_Ns$smooth$vect_N_smooth
  if (!length(n)) {
    n <- range(x$scaled_bn$N)
    warning("No stored smoothing grid; using empirical endpoints.")
  }
  x$simil_struct_Ns$smooth <- similaritystructure:::smooth_binomial(
    min(n), max(n), x$scaled_bn, num_repet
  )
  attr(x$scaled_bn, "num_repet") <- num_repet
  fit <- x$simil_struct_Ns$fitting
  x$simil_struct_Ns$Ns_diff <- if (is.null(fit)) NA_real_ else
    100 * (x$simil_struct_Ns$smooth$Ns_smooth - fit$Ns_fitting) / fit$Ns_fitting
  x
}

# Update one saved object or a list of saved similarity structures.
# Write a new file, never overwrite the input.
upgrade_cached_file <- function(input, output, num_repet) {
  if (!file.exists(input)) stop("Input file does not exist: ", input)
  if (file.exists(output)) stop("Output already exists: ", output)
  obj <- readRDS(input)
  updated <- if (!is.null(obj$simil_struct_Ns)) {
    upgrade_ss(obj, num_repet = num_repet)
  } else {
    lapply(obj, upgrade_ss, num_repet = num_repet)
  }
  saveRDS(updated, output)
  invisible(updated)
}
