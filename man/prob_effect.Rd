% Generated from documentation in R/prob_effect.R; regenerate with roxygen2.
\name{prob_effect}
\alias{prob_effect}
\title{Calculate the similar effect probability (pS)}
\usage{
prob_effect(fitting_A = NULL, fitting_B = NULL, smooth_A = NULL, smooth_B = NULL)
}
\arguments{
\item{fitting_A, fitting_B}{Fitting components from two results, obtained with \code{$simil_struct_Ns$fitting}.}
\item{smooth_A, smooth_B}{Smoothing components from two results, obtained with \code{$simil_struct_Ns$smooth}.}
}
\value{A numeric probability for one method, or a list with \code{ps_fit} and \code{ps_smooth} for both.}
\description{Calculate the minimum of the two tail probabilities evaluated at the similarity size of the other comparison.}
\details{
Integer subsample sizes are compared to the unrounded similarity sizes.
Use the same test, significance level and compatible sample-size supports.
Smoothing has finite support; fitting uses the infinite parametric tail.
pS is not an inferential p-value. A threshold such as 0.05 is a user-defined
decision criterion.
}
\examples{
\dontrun{
set.seed(123)
A <- rnorm(1000)
B <- rnorm(1000, 0.5)
C <- rnorm(1000, 0.9)
AB <- similarity_structure(A, B, N_fin = 300, test = "t-test", seed = 123)
AC <- similarity_structure(A, C, N_fin = 300, test = "t-test", seed = 124)
prob_effect(smooth_A = AB$simil_struct_Ns$smooth,
            smooth_B = AC$simil_struct_Ns$smooth)
if (!is.null(AB$simil_struct_Ns$fitting) &&
    !is.null(AC$simil_struct_Ns$fitting)) {
  prob_effect(fitting_A = AB$simil_struct_Ns$fitting,
              fitting_B = AC$simil_struct_Ns$fitting)
}
}
}
\seealso{\code{\link{similarity_structure}}}
