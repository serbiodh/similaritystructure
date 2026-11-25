# Load the analysis results distributed with the examples.
# This is an internal helper for the R Markdown files, not part of the public API.
# The expensive subsampling tests are not repeated.

read_example_cache <- function(example, test = "ks-test") {

  if (example == "mean") {
    if (!test %in% c("ks-test", "t-test"))
      stop("Choose 'ks-test' or 't-test'.")

    name <- if (test == "ks-test") "example_01_ks.RData" else "example_01_ttest.RData"
    file <- system.file("extdata", name, package = "similaritystructure")
    if (!nzchar(file)) stop("The saved normal-sample results are missing.")

    # Convert the stored empirical counts into the package's current output.
    source(system.file("scripts", "upgrade_cached_results.R",
                       package = "similaritystructure"), local = TRUE)
    return(lapply(readRDS(file), upgrade_ss, num_repet = 300L))
  }

  if (example == "bmi") {
    file <- system.file("extdata", "bmi_counts.csv",
                        package = "similaritystructure")
    if (!nzchar(file)) stop("The saved BMI results are missing.")
    counts <- read.csv(file)

    n_fin <- c("N-C" = 4500L, "F-N" = 1500L, "F-C" = 800L)
    out <- list()

    for (name in names(n_fin)) {
      empirical <- counts[counts$comparison == name, c("N", "num_nosig")]
      if (nrow(empirical) != 60L) stop("Invalid cached results: ", name)
      empirical$scaled <- empirical$num_nosig / empirical$num_nosig[1]
      attr(empirical, "num_repet") <- 300L

      estimate <- similaritystructure:::smooth_fitting_Ns(
        2L, n_fin[[name]], empirical, plotting = FALSE, num_repet = 300L
      )
      out[[name]] <- list(scaled_bn = empirical, simil_struct_Ns = estimate)
    }
    return(out)
  }

  if (example == "microbiome") {
    count_file <- system.file("extdata", "betadivers_counts.csv",
                              package = "similaritystructure")
    fit_file <- system.file("extdata", "betadivers_fitting.csv",
                            package = "similaritystructure")
    if (!nzchar(count_file) || !nzchar(fit_file))
      stop("The saved microbiome results are missing.")

    counts <- read.csv(count_file)
    fits <- read.csv(fit_file)
    n_fin <- c(Female = 70L, Male = 200L, no = 100L, yes = 100L)
    out <- list()

    for (name in names(n_fin)) {
      empirical <- counts[counts$comparison == name, c("N", "num_nosig")]
      if (nrow(empirical) < 3L) stop("Invalid cached results: ", name)
      empirical$scaled <- empirical$num_nosig / empirical$num_nosig[1]
      attr(empirical, "num_repet") <- 300L

      smooth <- similaritystructure:::smooth_binomial(
        N_init = 2L, N_fin = n_fin[[name]],
        propsig = empirical, num_repet = 300L
      )

      # Restore the original parametric curve and normalized coefficients.
      fitrow <- fits[fits$comparison == name, , drop = FALSE]
      if (nrow(fitrow) != 1L) stop("Missing parametric fit: ", name)
      raw <- as.numeric(unlist(fitrow[paste0("raw_p", 1:5)]))
      coef <- as.numeric(unlist(fitrow[paste0("ss_p", 1:5)]))
      n <- seq.int(2L, fitrow$N_fit_max)

      fitting <- list(
        vect_N_fit = n,
        fun_fit = (raw[1] + raw[2]*n + raw[3]*n^2 + raw[5]*n^3)*exp(-raw[4]*n),
        Sn = (coef[1] + coef[2]*n + coef[3]*n^2 + coef[5]*n^3)*exp(-coef[4]*n),
        param_Sn = coef,
        Ns_fitting = fitrow$Ns_fitting,
        Gof_fit = fitrow$Gof_fit
      )

      out[[name]] <- list(scaled_bn = empirical,
                          simil_struct_Ns = list(smooth = smooth, fitting = fitting))
    }
    return(out)
  }

  stop("Unknown example: ", example)
}
