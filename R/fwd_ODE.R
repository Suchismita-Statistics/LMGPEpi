#' fwd_ODE
#' Solve the forward deterministic and variational ODE system
#'
#' Solves a deterministic SIR-type system together with an auxiliary
#' matrix-valued system used to calculate the covariance structure of the
#' infected trajectory.
#'
#' @param beta Numeric. Infection rate parameter.
#' @param gamma Numeric. Recovery rate parameter.
#' @param rho Numeric. Limiting proportion of initially infected individuals relative to the initially susceptible population. Must satisfy 0 < rho < 1.
#' @param Tmax Numeric. Final observation time of epidemic
#' @param dt Numeric. time increment. Default is \code{0.1}.
#'
#' @return An object returned by \code{\link[deSolve]{ode}} containing the
#'   solution at each time point. The columns are \code{time}, \code{S},
#'   \code{I}, \code{f11}, \code{f12}, \code{f21}, and \code{f22}.
#'
#' @details
#' The function solves the following six-dimensional system over the interval
#' \eqn{0 \leq u \leq T_{\max}}:
#'
#' \deqn{\frac{d s(u)}{d u} = -\beta s(u)\iota(u),}
#' \deqn{\frac{d \iota(u)}{d u} =
#'   \beta s(u)\iota(u) - \gamma\iota(u),}
#' \deqn{\frac{d F_{11}(u)}{d u} =
#'   \beta\iota(u)\{F_{11}(u)-F_{21}(u)\},}
#' \deqn{\frac{d F_{12}(u)}{d u} =
#'   \beta\iota(u)\{F_{12}(u)-F_{22}(u)\},}
#' \deqn{\frac{d F_{21}(u)}{d u} =
#'   \beta s(u)\{F_{11}(u)-F_{21}(u)\}+\gamma F_{21}(u),}
#' \deqn{\frac{d F_{22}(u)}{d u} =
#'   \beta s(u)\{F_{12}(u)-F_{22}(u)\}+\gamma F_{22}(u).}
#'
#' The initial conditions are
#' \deqn{s(0)=1,\qquad \iota(0)=\rho,}
#' \deqn{F_{11}(0)=F_{22}(0)=1,\qquad
#' F_{12}(0)=F_{21}(0)=0.}
#'
#' The first two equations define the deterministic SIR trajectory, while the
#' last four equations define the fundamental matrix
#' \eqn{F(u)=(F_{ij}(u))_{i,j=1}^2}.
#'
#' @seealso \code{\link{calc_K}}, \code{\link[deSolve]{ode}}
#'
#' @export

fwd_ODE = function(beta, gamma, rho, Tmax, dt = 0.1)
{
  time_pts = seq(0, Tmax, by = dt)
  parameters = c(beta, gamma, rho)
  state  = c(S = 1, I = rho, f11 = 1, f12 = 0, f21 = 0, f22 = 1)

  Lorenz = function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      ## rate of change
      dS = -beta * S * I
      dI = beta * S * I - gamma * I
      df11 = beta * I * (f11 - f21)
      df12 = beta*I*(f12 - f22)
      df21 = beta * S * (f11 - f21) + gamma * f21
      df22 = beta * S * (f12 - f22) + gamma * f22
      list(c(dS, dI, df11, df12, df21, df22))
    })
  }
  out <- deSolve::ode(
    y = state,
    times = time_pts,
    func = Lorenz,
    parms = parameters
  )
  return(out)
}


