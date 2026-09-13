#' summary_stan_list_adj
#'
#' It runs summary_stan_adj functions for all parameters - \eqn{\beta}{beta}, \eqn{\gamma}{gamma}, \eqn{\rho}{rho} (for probability model also \eqn{p}) and returns in the form of a matrix. It can also return the CIs that do not contain the true parameters if notcap = TRUE.
#' @param result_list List of \code{stanfit} objects that we want to summarize.
#' @param true_values  Vector of true values of the parameters in the order \eqn{\beta}{beta}, \eqn{\gamma}{gamma}, \eqn{\rho}{rho} respectively.
#' @param prob Logical. If TRUE, it will consider the probability adjusted model and give estimates of \eqn{p}. Default is False.
#' @param not_cap Logical. If TRUE, it returns the pairs of 95% CI for all four parameters, which do not contain the true parameter. Otherwise, it returns mean of the means, standard deviations and the proportion of 95% CIs captures the true parameter for all four parameters.
#' @return Mean of posterior means, sd of posterior means, mean of the posterior sds, \eqn{95%} frequentist coverage, mean of ESS/sec of the list for all parameters.
#' @seealso \code{\link{summary_stan_adj}}
#'
#' @export


summary_stan_list_adj <- function(result_list, true_value, prob = FALSE, notcap = FALSE)
{
  if(prob)
  {
    par = c("b", "g", "r", "R0", "p")
  }else{
    par = c("b", "g", "r", "R0")
  }
  draws_list <- lapply(result_list, function(x) rstan::extract(x$fit_full, pars = par))

  Z_hats <- lapply(result_list, function(x) x$Z_hat)
  M <- length(Z_hats[[1]])
  wts        <- lapply(Z_hats, function(x) wnorm_calc(x, M)$w_norm)

  ess_per_sec_mat <- t(mapply(essinfo_per_replicate, result_list, wts,
                              MoreArgs = list(par = par)))  # rows = replicates, cols = b,g,r,R0

  para_beta  <- summary_stan_adj(draws_list, wts, ess_per_sec_mat[, "b"],  index = 1, true_value = true_value[1], not_cap = notcap)
  para_gamma <- summary_stan_adj(draws_list, wts, ess_per_sec_mat[, "g"],  index = 2, true_value = true_value[2], not_cap = notcap)
  para_rho   <- summary_stan_adj(draws_list, wts, ess_per_sec_mat[, "r"],  index = 3, true_value = true_value[3], not_cap = notcap)
  para_R0    <- summary_stan_adj(draws_list, wts, ess_per_sec_mat[, "R0"], index = 4, true_value = true_value[1] / true_value[2], not_cap = notcap)

  if(!prob)
  {
    colnam <- c("mn", "sdofmean", "sd", "95%cvg", "ESS/s")
    if (notcap == FALSE) {
      mat <- matrix(c(para_beta, para_gamma, para_rho, para_R0), ncol = 5, byrow = TRUE)
      colnames(mat) <- colnam
      rownames(mat) <- c("beta", "gamma", "rho", "R0")
      return(mat)
    } else {
      not_captured <- list(beta = para_beta, gamma = para_gamma, rho = para_rho, R0 = para_R0)
      return(not_captured)
    }
  }else{
    para_p   <- summary_stan_adj(draws_list, wts, ess_per_sec_mat[, "p"],  index = 4, true_value = true_value[4], not_cap = notcap)
    colnam <- c("mn", "sdofmean", "sd", "95%cvg", "ESS/s")
    if (notcap == FALSE) {
      mat <- matrix(c(para_beta, para_gamma, para_rho, para_p, para_R0), ncol = 5, byrow = TRUE)
      colnames(mat) <- colnam
      rownames(mat) <- c("beta", "gamma", "rho", "p", "R0")
      return(mat)
    } else {
      not_captured <- list(beta = para_beta, gamma = para_gamma, rho = para_rho, p = para_p, R0 = para_R0)
      return(not_captured)
    }
  }


}
