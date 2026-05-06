#' PERMANOVA test for subsampled distance matrices
#'
#' @param n1sub First subsample distance matrix.
#' @param n2sub Second subsample distance matrix.
#' @param subset_choice Character. One of `"Female"`, `"Male"`, `"yes"` or `"no"`.
#' @param permutations Number of permutations passed to `vegan::adonis2()`.
#'
#' @return First row of the `adonis2()` table, including a `p.value` column.
#'
#' @keywords internal

adonis_ss <- function(n1sub, n2sub,
                      subset_choice = c("Female", "Male", "yes", "no"),
                      permutations = 999) {

  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for adonis_ss(). Please install it.")
  }

  subset_choice <- match.arg(subset_choice)

  D <- rbind(n1sub, n2sub)
  D <- D[, rownames(D)]

  if (subset_choice %in% c("Female", "Male")) {
    meta <- data.frame(
      smoker = c(rep("yes", nrow(n1sub)), rep("no", nrow(n2sub)))
    )
  } else {
    meta <- data.frame(
      gender = c(rep("Female", nrow(n1sub)), rep("Male", nrow(n2sub)))
    )
  }

  form <- stats::as.formula(paste0("D ~ ", colnames(meta)))

  fit <- vegan::adonis2(
    formula = form,
    data = meta,
    permutations = permutations
  )

  tab <- as.data.frame(fit)
  tab$term <- rownames(tab)
  tab <- tab[1, ]
  tab$p.value <- tab$`Pr(>F)`

  tab
}
