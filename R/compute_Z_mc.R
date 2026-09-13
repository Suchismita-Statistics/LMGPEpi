#' compute_Z_mc
#' Using the mean and the covariance matrix of a multivariate normal, this function calculates the Monte Carlo estimates simplex probability.
#'
#' @param mn_vec The mean vector of the multivariate normal.
#' @param cov_mat The covariance matrix of the multivariate normal.
#' @param n_mc The number of Monte Carlo samples used to find the estimate. Default is 4000.
#' @param jitter A small non-negative value added to the diagonal of the multivariate normal for numerical stabilization. Default is \code{1e-10}.
#' @return A list with two components:
#' \describe{
#'   \item{Z}{The Monte Carlo estimate of the simplex probability.}
#'   \item{se}{The Monte Carlo standard error of the estimate.}
#' }
#' @export

compute_Z_mc <- function(mn_vec, cov_mat, n_mc = 4000, jitter = 1e-10) {
  cov_mat <- 0.5 * (cov_mat + t(cov_mat))
  diag(cov_mat) <- diag(cov_mat) + jitter
  L <- tryCatch(chol(cov_mat), error = function(e) NULL)
  if (is.null(L)) {
    diag(cov_mat) <- diag(cov_mat) + 1e-6
    L <- tryCatch(chol(cov_mat), error = function(e) NULL)
  }
  if (is.null(L)) return(list(Z = NA_real_, se = NA_real_))

  z   <- matrix(rnorm(n_mc * length(mn_vec)), n_mc, length(mn_vec))
  smp <- sweep(z %*% L, 2, mn_vec, "+")
  in_R <- (rowSums(smp > 0) == ncol(smp)) & (rowSums(smp) < 1)
  Zhat <- mean(in_R)
  list(Z = Zhat, se = sqrt(Zhat * (1 - Zhat) / n_mc))
}
