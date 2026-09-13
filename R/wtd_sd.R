#' wtd_sd
#' Calculates the weighted standard deriviation
#'
#' @param x vector of values
#' @param w vector of sd
#'
#' @return A numeric scalar giving the weighted standard deviation.
#'
#' @details
#' The weighted mean is calculated as
#' \deqn{\bar{x}_w = \sum_i w_i x_i,}
#' and the weighted standard deviation is calculated as
#' \deqn{
#' \sqrt{\frac{\sum_i w_i(x_i-\bar{x}_w)^2}
#' {1-\sum_i w_i^2}}.
#' }
#' This formula assumes that the weights are normalized so that
#' \eqn{\sum_i w_i = 1}.
#' @export

wtd_sd <- function(x, w) {
  m <- sum(w * x)
  sqrt(sum(w * (x - m)^2) / (1 - sum(w^2)))
}
