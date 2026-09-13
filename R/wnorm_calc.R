#' wnorm_calc
#' Taking the estimates of estimates of the simplex probabilities, it calculates the normalized weights which are inversely proportional to the probabilities.
#' @param Z_hat The vector of Monte-carlo estimates of the simplex probability.
#' @param M Positive integer specifying the number of weights to return.
#'
#' @return A list with two components:
#' \describe{
#'   \item{ok}{Logical vector indicating which entries of \code{Z_hat} are
#'   non-missing and positive.}
#'   \item{w_norm}{Numeric vector of normalized weights. Invalid entries are
#'   assigned weight zero.}
#' }
#' @export

wnorm_calc = function(Z_hat, M)
{
  ok <- !is.na(Z_hat) & Z_hat > 0
  log_w <- rep(NA_real_, M)
  log_w[ok] <- -log(Z_hat[ok])
  log_w[ok] <- log_w[ok] - max(log_w[ok])
  w <- rep(0, M)
  w[ok] <- exp(log_w[ok])
  w_norm <- w / sum(w)

  return(list(ok = ok, w_norm = w_norm))
}
