#' essinfor_per_replicate
#'
#' This function computes the effective sample size of one fitted replicate after the
#' importance-sampling correction for the simplex probability
#' \eqn{\mathbb{C}(\theta)}, expressed per second of total run time.
#'
#' @import rstan
#'
#' @param result The \code{stanfit} output of LMGP function.
#' @param w_norm Numeric vector of normalised importance weights
#'   \eqn{\tilde{w}_i} for that replicate, one per posterior draw, as produced
#'   by \code{\link{wnorm_calc}}. Its length is the number of draws \eqn{L}.
#' @param par Character vector of parameter names to report, passed to
#'   \code{\link[rstan]{summary}} (for example \code{c("b", "g", "r", "R0")},
#'   with \code{"p"} added under the reporting-probability model).
#'
#'
#' @details
#' MCMC effective sample size alone overstates the information available here,
#' because the draws are subsequently reweighted by
#' \eqn{w_i = 1 / \widehat{\mathbb{C}}(\theta_i)} and unequal weights discard
#' part of the sample. The adjustment multiplies Stan's per-parameter
#' \code{n_eff} by the weight-degeneracy factor
#' \deqn{\frac{1}{L}\left(\sum_{i=1}^{L} \tilde{w}_i^2\right)^{-1},}
#' the Kish effective sample size of the normalised weights
#' \eqn{\tilde{w}_i} as a fraction of the number of draws \eqn{L}, giving
#' \deqn{\mathrm{ESS} = \frac{\mathrm{ESS}_{\mathrm{MCMC}}}{L}
#'       \cdot \frac{1}{\sum_{i=1}^{L} \tilde{w}_i^2}.}
#' Weights that are close to uniform leave \code{n_eff} nearly unchanged;
#' weights that degenerate onto a few draws shrink it sharply.
#'
#' The result is then divided by \code{result$total_time}, the elapsed time of
#' the whole corrected pipeline -- Stan sampling including warmup, plus the
#' per-draw \eqn{\widehat{\mathbb{C}}(\theta_i)} loop. It is therefore the
#' throughput of the corrected procedure end to end, and is not comparable with
#' an ESS per second computed from post-warmup sampling time alone; the two
#' should be reported separately rather than merged.
#'
#'  @seealso \code{\link{wnorm_calc}}, \code{\link{compute_Z_mc}},
#'   \code{\link{summary_stan_list_adj}}
#'
#'   @export


essinfo_per_replicate <- function(result, w_norm, par) {
  total_time <- result$total_time

  L        <- length(w_norm)
  ess_is   <- 1 / sum(w_norm^2)
  ess_mcmc <- rstan::summary(result$fit_full, pars = par)$summary[, "n_eff"]
  ess_comb <- (ess_mcmc / L) * ess_is
  ess_comb / total_time
}
