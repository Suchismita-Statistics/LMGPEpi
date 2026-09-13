#' noisy_data_sim
#'
#' Generates a data set with under-reported infection counts simulated via Sellke construction, under-reported by probability \eqn{p}.
#'
#' @import MASS
#' @import DSA.CountData
#'
#' @param N Integer. Intial number of susceptible.
#' @param beta Numeric. Infection rate parameter.
#' @param gamma Numeric. Recovery rate parameter.
#' @param rho Numeric. Limiting proportion of initially infected individuals relative to the initially susceptible population. Must satisfy 0 < rho < 1.
#' @param p Numeric. The reporting probability parameter. Must satisfy 0 < p < 1.
#' @param Tmax Numeric. Final observation time of epidemic.
#' @param nu Default is NULL, which simulate usual SIR CTMC dynamics. If specified some numerical value, then, it simulates data from frailty model. It is the parameter denoting the standard deviation of the frailty variable which is assumed to follow Gamma distribution.
#'
#' @return A data set with three columns - days, exact infection counts and under-reported counts with reporting probability \eqn{p}. The true trajectory is simulated via Sellke construction.
#' @seealso \code{\link[DSA.CountData]{sellke}} for reference of \eqn{\nu}.
#' @export

noisy_data_sim = function(N, beta, gamma, rho, p, Tmax, nu = NULL)
{
  count_data = sellke_counts(
    N = N,
    beta = beta,
    gamma = gamma,
    r = rho,
    Tmax = Tmax,
    nu = nu
  )
  true_noise = rbinom(nrow(count_data), count_data[, 2], p)

  final_data = cbind(count_data, underreported_counts = true_noise)

  return(final_data)
}
