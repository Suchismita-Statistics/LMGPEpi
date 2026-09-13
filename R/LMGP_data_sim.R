#' LMGP_data_sim
#'
#' Generates a data set of exact infection and recovery times according to the LMGP approximation, where the survival function of time to infection follows a Gaussian process.
#'
#' @import MASS
#' @param N  Integer. Intial number of susceptible.
#' @param beta Numeric. Infection rate parameter.
#' @param gamma Numeric. Recovery rate parameter.
#' @param rho Numeric. Limiting proportion of initially infected individuals relative to the initially susceptible population. Must satisfy 0 < rho < 1.
#' @param Tmax Numeric. Final observation time of epidemic
#' @param dt Numeric. time increment. Default is 0.1.
#' @param iter Number of dataset wish to simulate.
#'
#' @export

LMGP_data_sim = function(N,
                         beta,
                         gamma,
                         rho,
                         Tmax,
                         dt = 0.1,
                         iter = 1)
{
  full_grid = seq(0, Tmax, dt)
  len = rep(c(4, 2), length(full_grid) / 2)
  simp_coeff = c(1, len[-length(len)], 1)


  res = calc_K_sim(beta, gamma, rho, Tmax, dt = 0.1, simp_coeff, N)
  mn = res[[1]]
  cov_mat = res[[2]]

  simulation_mat = mvrnorm(iter, mn, cov_mat)

  data_gen = matrix(0, ncol = iter, nrow = N)

  for (i in 1:iter)
  {
    s_T = tail(simulation_mat[i, ], 1)

    u = runif(N)
    not_inf_during_epi = sum(as.numeric(u < s_T))

    u_infect = runif((N - not_inf_during_epi))
    temp = u_infect * (1 - s_T) + s_T
    samples = approx(x = simulation_mat[i, ], y = full_grid, xout = temp)

    infect_times = samples$y
    infect_times[is.na(infect_times)] = Tmax
    data_gen[, i] = sort(c(infect_times, rep(Tmax, not_inf_during_epi)))
  }

  return(data_gen)
}
