#' calc_K
#'
#' Calculate the mean trajectory and covariance matrix
#'
#' Uses the solution of \code{\link{fwd_ODE}} to calculate a mean trajectory
#' and its associated covariance matrix through numerical integration.
#'
#' @param beta Numeric. Infection rate parameter.
#' @param gamma Numeric. Recovery rate parameter.
#' @param rho Numeric. Limiting proportion of initially infected individuals relative to the initially susceptible population. Must satisfy 0 < rho < 1.
#' @param Tmax Numeric. Final observation time of epidemic
#' @param dt Numeric. time increment. Default is \code{0.1}.
#' @param simpson_coeff Numeric vector of Simpson's-rule integration
#'   coefficients for the specified odd-numbered time grid.
#' @param N Positive scaling factor, typically the population size. The
#'   resulting covariance matrix is divided by \code{N}.
#'
#' @return A list with two components:
#' \describe{
#'   \item{\code{mn}}{Numeric vector containing the mean infected trajectory.}
#'   \item{\code{cov}}{Covariance matrix corresponding to the time points of
#'   the ODE solution.}
#' }
#'
#' @details
#' The covariance matrix is based on the Gaussian-process approximation
#' \deqn{
#' \sqrt{N}\left\{\frac{S(t)}{N}-s(t)\right\}
#' \xrightarrow{d} \mathcal{GP}(0,\mathcal{K}),
#' }
#' where \eqn{S(t)/N} is the stochastic susceptible proportion and
#' \eqn{s(t)} is its deterministic limit.
#'
#' The covariance function \eqn{\mathcal{K}} is
#' \deqn{
#' \mathcal{K}(t,r)=
#' \int_0^{\min(t,r)}
#' \left[
#' \beta s(u)\iota(u)
#' \{\phi_{11}(t,u)-\phi_{12}(t,u)\}
#' \{\phi_{11}(r,u)-\phi_{12}(r,u)\}
#' +\gamma\iota(u)\phi_{12}(t,u)\phi_{12}(r,u)
#' \right]\,du.
#' }
#'
#' The functions \eqn{s(u)} and \eqn{\iota(u)} satisfy
#' \deqn{\frac{d s(u)}{d u}=-\beta s(u)\iota(u),}
#' \deqn{\frac{d\iota(u)}{d u}=\beta s(u)\iota(u)-\gamma\iota(u),}
#' with initial conditions
#' \deqn{s(0)=1,\qquad \iota(0)=\rho.}
#'
#' For each fixed \eqn{t}, the functions
#' \eqn{\phi_{11}(t,u)} and \eqn{\phi_{12}(t,u)} satisfy
#' \deqn{
#' \frac{\partial\phi_{11}(t,u)}{\partial u}
#' =\beta\iota(u)\phi_{11}(t,u)
#' -\beta\iota(u)\phi_{12}(t,u),
#' }
#' \deqn{
#' \frac{\partial\phi_{12}(t,u)}{\partial u}
#' =\beta s(u)\phi_{11}(t,u)
#' -\beta s(u)\phi_{12}(t,u)
#' +\gamma\phi_{12}(t,u),
#' }
#' with
#' \deqn{
#' \phi_{11}(t,t)=1,\qquad \phi_{12}(t,t)=0.
#' }
#'
#' Rather than solving this two-dimensional system separately for every
#' terminal time \eqn{t}, the function \code{\link{fwd_ODE}} solves the
#' associated fundamental-matrix system once. If
#' \eqn{F(u)} is the fundamental matrix satisfying
#' \deqn{F'(u)=B(u)F(u),\qquad F(0)=I_2,}
#' then
#' \deqn{
#' \begin{pmatrix}
#' \phi_{11}(t,u)\\
#' \phi_{12}(t,u)
#' \end{pmatrix}
#' =
#' F(u)F(t)^{-1}
#' \begin{pmatrix}1\\0\end{pmatrix}.
#' }
#'
#' The integral defining \eqn{\mathcal{K}(t,r)} is approximated numerically
#' using Simpson's one-third rule with step size \code{dt} and coefficients
#' supplied through \code{simpson_coeff}. The number of
#' subintervals between zero and the target time must be even. Consequently,
#' the target time must correspond to an odd-numbered grid point. For the resulting time grid, the
#' returned component \code{mn} contains the deterministic values
#' \eqn{s(t)}, while \code{cov} contains the covariance matrix with entries
#' approximately equal to \eqn{\mathcal{K}(t_i,t_j)/N}.
#'
#' @seealso \code{\link{fwd_ODE}}
#'
#' @export

calc_K = function(beta, gamma, rho, Tmax, dt, simpson_coeff, N)
{
  ode_solve = fwd_ODE(beta, gamma, rho, Tmax, dt)
  full_index = seq(0, Tmax, dt)
  ful_len = nrow(ode_solve)

  int_index = which(ode_solve[, 1] %% 1 == 0)[-1]

  inv_store = matrix(0, ncol = 2, nrow = ful_len)
  temp = ode_solve[, 4:7]
  det = temp[, 1]*temp[, 4] - temp[, 2]*temp[, 3]
  inv_store[, 1] = temp[, 4]/det
  inv_store[, 2] = - temp[, 3]/det

  phi11 = matrix(0, ncol = ful_len, nrow = nrow(ode_solve))
  phi12 = matrix(0, ncol = ful_len, nrow = nrow(ode_solve))
  for(i in 1:ful_len)
  {
    phi11[, i] = ode_solve[, 4]*inv_store[i, 1] + ode_solve[, 5]*inv_store[i, 2]
    phi12[, i] = ode_solve[, 6]*inv_store[i, 1] + ode_solve[, 7]*inv_store[i, 2]
  }

  cov_mat = matrix(0, ncol = ful_len, nrow = ful_len)

  for(i in 2:ful_len)
  {

    upto = i
    s_i = ode_solve[1:upto, 2]
    iota_i = ode_solve[1:upto, 3]
    phi11_i = phi11[1:upto, i]
    phi12_i = phi12[1:upto, i]

    diff_i = phi11_i - phi12_i

    simp_coeff = rep(1, upto)
    simp_coeff[1:(upto - 1)] = simpson_coeff[1:(upto - 1)]

    cov_mat[i, i] = sum(beta*s_i*iota_i*diff_i*diff_i*dt*simp_coeff)/3 + sum(dt*gamma*phi12_i*phi12_i*iota_i*simp_coeff)/3
    if(i < ful_len)
    {
      for(j in (i+1):ful_len)
      {
        diff_j = phi11[1:upto, j] - phi12[1:upto, j]
        phi12_j = phi12[1:upto, j]
        cov_mat[i, j] = sum(dt*beta*s_i*iota_i*diff_i*diff_j*simp_coeff)/3 + sum(dt*gamma*phi12_i*phi12_j*iota_i*simp_coeff)/3

        cov_mat[j, i] = cov_mat[i, j]
      }
    }
  }
  return(list(mn = ode_solve[, 2], cov = cov_mat/N))
}
