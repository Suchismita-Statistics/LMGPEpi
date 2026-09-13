#' zhat_stability_replicate
#' Assess Monte Carlo stability of estimated normalizing constants
#'
#' Summarizes the estimated values \code{Z_hat} from one simulation replicate.
#' The function calculates the mean, standard deviation, minimum, Monte Carlo
#' standard error, coefficient of variation at the minimum estimate, effective
#' sample size fraction, and proportion of valid estimates.
#'
#' @param x A list containing a numeric component named \code{Z_hat}.
#' @param n_mc Positive integer specifying the number of Monte Carlo samples
#'   used to calculate the Monte Carlo standard error.
#'
#' @return A named numeric vector containing:
#' \describe{
#'   \item{\code{mean_Zhat}}{Mean of the valid \code{Z_hat} values.}
#'   \item{\code{sd_Zhat}}{Standard deviation of the valid \code{Z_hat} values.}
#'   \item{\code{mean_se_mc}}{Mean Monte Carlo standard error.}
#'   \item{\code{min_Zhat}}{Minimum valid \code{Z_hat} value.}
#'   \item{\code{ess_frac}}{Effective sample size as a fraction of the total
#'   number of estimates.}
#' }
#'
#' @export



zhat_stability_replicate <- function(x, n_mc = 2000) {
  Zhat <- x$Z_hat

  ok      <- !is.na(Zhat) & Zhat > 0

  Zok    <- Zhat[ok]
  mean_Z <- mean(Zok)
  sd_Z   <- sd(Zok)
  min_Z  <- min(Zok)

  se_mc   <- sqrt(Zok * (1 - Zok) / n_mc)
  mean_se <- mean(se_mc)


  # recompute normalized weights from Z_hat directly -- same procedure as
  # inside much_bhi(), no dependence on a stored w_norm field
  log_w     <- rep(NA_real_, length(Zhat))
  log_w[ok] <- -log(Zhat[ok])
  log_w[ok] <- log_w[ok] - max(log_w[ok])
  w         <- rep(0, length(Zhat))
  w[ok]     <- exp(log_w[ok])
  w_norm    <- w / sum(w)

  ess_frac <- (1 / sum(w_norm^2, na.rm = TRUE)) / length(Zhat)

  c(
    mean_Zhat = mean_Z,
    sd_Zhat = sd_Z,
    mean_se_mc = mean_se,
    min_Zhat = min_Z,
    ess_frac = ess_frac
  )
}

