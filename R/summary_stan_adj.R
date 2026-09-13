#' summary_stan_adj
#' Summarize weighted posterior draws and calculates the posterior coverage and the ESS/s for a single parameter.
#'
#' @param draws_list List of \code{stanfit} objects.
#' @param wts List of normalized weights.
#' @param ess_per_sec_par Numeric vector containing effective sample size per
#'   second values for the simulation replicates for a the parameter of interest.
#'  @param index It identifies the parameter we are interested to summarize:  \code{1} = \eqn{\beta}, \code{2} = \eqn{\gamma},
#'   \code{3} = \eqn{\rho}, and \code{4} = \eqn{R_0}.
#'   @param true_value The true value of the parameter.
#' @param not_cap Logical. If \code{TRUE}, it returns the pairs of 95% CI, which do not contain the true parameter. Otherwise, it returns mean of the means, standard deviations and the proportion of 95% CIs captures the true parameter.
#' @return Mean of posterior means, sd of posterior means, mean of the posterior sds, \eqn{95%} frequentist coverage, mean of ESS/sec summarized over the different \code{stanfit} results.
#' @seealso \code{\link{summary_stan_adj_list}}
#' @export


summary_stan_adj <- function(draws_list, wts, ess_per_sec_par, index, true_value, not_cap = FALSE)
{
  temp1 <- mapply(function(d, w) Hmisc::wtd.quantile(d[[index]], w, probs = 0.025, normwt = TRUE),
                  draws_list, wts)
  temp2 <- mapply(function(d, w) Hmisc::wtd.quantile(d[[index]], w, probs = 0.975, normwt = TRUE),
                  draws_list, wts)
  temp3 <- as.numeric(temp1 < true_value & true_value < temp2)

  if (not_cap == TRUE) {
    ret <- (1:length(draws_list))[which(temp3 == 0)]
    mat <- matrix(c(temp1[ret], temp2[ret]), ncol = 2)
    return(list(ret, mat))
  } else {
    foo  <- mapply(function(d, w) sum(w * d[[index]]), draws_list, wts)
    mn   <- mean(foo)
    sd_of_mn <- sd(foo)
    sd   <- mean(mapply(function(d, w) wtd_sd(d[[index]], w), draws_list, wts))
    covg <- mean(temp3)
    esspersec_val <- mean(ess_per_sec_par)
    return(c(mn, sd_of_mn, sd, covg, esspersec_val))
  }
}
