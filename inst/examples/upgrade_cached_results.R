# Convert cached results to the current field names without changing values.
# This also allows examples to read results saved by earlier package versions.
upgrade_ss <- function(x) {
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
  x
}
