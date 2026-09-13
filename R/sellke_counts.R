#' sellke_counts
#'
#' Generates a data set with exact infection counts simulated via Sellke construction.
#'
#' @import MASS
#' @import DSA.CountData
#'
#' @param N Integer. Intial number of susceptible.
#' @param beta Numeric. Infection rate parameter.
#' @param gamma Numeric. Recovery rate parameter.
#' @param rho Numeric. Limiting proportion of initially infected individuals relative to the initially susceptible population. Must satisfy 0 < rho < 1.
#' @param Tmax Numeric. Final observation time of epidemic
#' @param nu Default is NULL, which simulate usual SIR CTMC dynamics. If specified some numerical value, then, it simulates data from frailty model. It is the parameter denoting the standard deviation of the frailty variable which is assumed to follow Gamma distribution.
#'
#' @return a data set with two columns - days of infection and the infection counts.
#' @seealso \code{\link[DSA.CountData]{sellke}} for reference of \eqn{\nu}.
#' @export

sellke_counts = function(N, beta, gamma, rho, Tmax, nu = NULL)
{
  data_Sellke = DSA.CountData::sellke(
    n = N,
    beta = beta,
    gamma = gamma,
    rho = rho,
    Tmax = Tmax,
    nu = nu
  )
  initial_sus = data_Sellke[data_Sellke[, 1] != 0, ]

  M = nrow(data_Sellke) - N

  duplicates <- duplicated(initial_sus[which(initial_sus[, 1] < Tmax), 1])
  initial_sus = initial_sus[order(initial_sus[, 1]), ]
  infect_during_ep = subset(initial_sus, initial_sus[, 1] < Tmax)
  t = infect_during_ep[, 1]

  infection_days = as.numeric(names(table(floor(t)))) + 1
  infection_count = as.vector(table(floor(t)))
  inf_time_final = numeric(length = Tmax)


  inf_time_final[infection_days] = infection_count
  temp = matrix(c(infection_days, inf_time_final), ncol = 2)
  colnames(temp) = c("days_reporting", "true_counts")
  return(temp)
}
