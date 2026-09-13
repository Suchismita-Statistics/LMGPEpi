#' beta_params
#'
#' Given mean and sd, this function finds the value of the parameters for a Beta distribution.
#'
#' @param mean Mean of the Beta distribution.
#' @param sd Standard deviation of the Beta Distribution.
#'
#' @return The two parameters of the Beta distribution determined by specified mean and sd.
#' @export


beta_params <- function(mean, sd) {
  v <- sd^2
  if (mean <= 0 || mean >= 1) stop("mean must be in (0, 1)")
  if (v >= mean * (1 - mean))
    stop(sprintf("sd too large: variance must be < %.6f for this mean", mean * (1 - mean)))

  nu <- mean * (1 - mean) / v - 1      # a + b
  c(a = mean * nu, b = (1 - mean) * nu)
}
