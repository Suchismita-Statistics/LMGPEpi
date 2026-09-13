#' LMGP_ppc
#'
#' Posterior predictive check for models fitted using LMGP method.
#'
#' @import rstan
#' @import ggplot2
#'
#' @param fit_full The stanfit object from the output of LMGP function.
#' @param data_obs The final data set (after pre-processing), which was used to fitting. A two-column matrix with first column days and second column is corresponding counts.
#' @param N Initial number of susceptible
#' @param raw_data A two-column matrix with first column days and second column is raw counts if it differs from the processed data mentioned in \code{data_obs}.
#' @param w_norm Normalized weights. If \code{NULL}, equal weights will be used.
#' @param quantile_probs A vector of probabilities, will be used to make credible intervals in posterior predictive checks. The default is \code{0.025, 0.5, 0.975}.
#'
#' @export


LMGP_ppc <- function(fit_full,
                     data_obs,
                     N,
                     prob     = TRUE,
                     raw_data = NULL,
                     w_norm   = NULL,
                     quantile_probs    = c(0.025, 0.5, 0.975)) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  stopifnot(requireNamespace("rstan", quietly = TRUE))
  stopifnot(is.logical(prob), length(prob) == 1L, !is.na(prob))

  ## -- 1. Posterior draws of u (M x (Tmax+1)) and, if present, p ---------
  ## Extracted in ONE call so both come from the same cached permutation
  ## of fit_full (rstan fixes one permutation per stanfit object).
  pars_needed <- if (prob)
    c("u", "p")
  else
    "u"
  draws   <- rstan::extract(fit_full, pars = pars_needed)
  u_draws <- draws$u
  M    <- nrow(u_draws)
  Tmax <- ncol(u_draws) - 1
  stopifnot(Tmax == nrow(data_obs))

  if (prob) {
    p_draws <- draws$p
    stopifnot(length(p_draws) == M)
  } else {
    ## No reporting probability in this model: everything is observed,
    ## which is the p == 1 special case of the formulas below.
    p_draws <- rep(1, M)
  }

  ## -- 2. Map each (u, p) draw to an implied count trajectory ------------
  ## Mirrors u_tilde in the Stan model block exactly:
  ##   u_tilde[1]      = p * (1 - sum(u[1:Tmax]))
  ##   u_tilde[2:Tmax] = p * u[1:(Tmax-1)]
  ## (u_tilde[Tmax+1], tied to N - K, is not part of this series.)
  pred_counts <- matrix(NA_real_, nrow = M, ncol = Tmax)
  pred_counts[, 1] <- N * p_draws * (1 - rowSums(u_draws[, 1:Tmax, drop = FALSE]))
  if (Tmax > 1) {
    pred_counts[, 2:Tmax] <- sweep(
      u_draws[, 1:(Tmax - 1), drop = FALSE],
      MARGIN = 1,
      STATS = N * p_draws,
      FUN = "*"
    )
  }

  ## -- 3. Importance weights -----------------------------------------------
  if (!is.null(w_norm) && length(w_norm) == M) {
    wts <- w_norm
  } else {
    warning("w_norm not supplied (or wrong length) -- using uniform weights.")
    wts <- rep(1 / M, M)
  }

  ## -- 4. Weighted quantile helper -------------------------------------------
  wtd_quantile <- function(x, w, quantile_probs) {
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    cw <- cumsum(w) / sum(w)
    sapply(quantile_probs, function(pp)
      x[which(cw >= pp)[1]])
  }

  interval_summary <- t(apply(
    pred_counts,
    2,
    wtd_quantile,
    w = wts,
    quantile_probs = quantile_probs
  ))
  colnames(interval_summary) <- c("lo95", "median", "hi95")

  cum_pred_counts <- t(apply(pred_counts, 1, cumsum))          # M x Tmax
  cum_summary <- t(apply(
    cum_pred_counts,
    2,
    wtd_quantile,
    w = wts,
    quantile_probs = quantile_probs
  ))
  colnames(cum_summary) <- c("lo95", "median", "hi95")

  ## -- 5. Assemble plotting data frames --------------------------------------
  time_idx <- data_obs[, 1]

  ## Optional raw (pre-MA) series, normalized to a 2-col (time, raw) df.
  raw_df <- NULL
  if (!is.null(raw_data)) {
    if (is.null(dim(raw_data))) {
      stopifnot(length(raw_data) == nrow(data_obs))
      raw_df <- data.frame(time = time_idx, raw = as.numeric(raw_data))
    } else {
      stopifnot(nrow(raw_data) == nrow(data_obs))
      raw_df <- data.frame(time = time_idx, raw = as.numeric(raw_data[, 2]))
    }
  }

  df_interval <- data.frame(
    time   = time_idx,
    ma     = data_obs[, 2],
    lo95   = interval_summary[, "lo95"],
    median = interval_summary[, "median"],
    hi95   = interval_summary[, "hi95"]
  )

  df_cum <- data.frame(
    time     = time_idx,
    observed = cumsum(data_obs[, 2]),
    lo95     = cum_summary[, "lo95"],
    median   = cum_summary[, "median"],
    hi95     = cum_summary[, "hi95"]
  )

  ## -- 6. Labels differ only in wording (reported vs. all infections) --------
  y_lab_new <- if (prob)
    "New reported infections"
  else
    "New infections"
  y_lab_cum <- if (prob)
    "Cumulative reported infections"
  else
    "Cumulative infections"

  ## -- 7. Plots ----------------------------------------------------------------
  p_interval <- ggplot2::ggplot(df_interval, ggplot2::aes(x = time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo95, ymax = hi95),
                         fill = "steelblue",
                         alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(y = median),
                       color = "steelblue",
                       linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(y = ma, color = "3-day MA"), size = 2)

  if (!is.null(raw_df)) {
    p_interval <- p_interval +
      ggplot2::geom_point(
        data = raw_df,
        ggplot2::aes(x = time, y = raw, color = "Raw"),
        shape = 17,
        size = 2.5,
        alpha = 1
      )
  }

  p_interval <- p_interval +
    ggplot2::scale_color_manual(name   = NULL,
                                values = c("3-day MA" = "black", "Raw" = "darkorange3")) +
    ggplot2::labs(title = "Predicted New Infections", x = "Time", y = y_lab_new) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "bottom")

  p_cum <- ggplot2::ggplot(df_cum, ggplot2::aes(x = time)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lo95, ymax = hi95),
      fill = "darkorange",
      alpha = 0.25
    ) +
    ggplot2::geom_line(ggplot2::aes(y = median),
                       color = "darkorange",
                       linewidth = 1) +
    ggplot2::geom_point(ggplot2::aes(y = observed),
                        color = "black",
                        size = 2) +
    ggplot2::labs(title = "Predicted Cumulative Infections", x = "Time", y = y_lab_cum) +
    ggplot2::theme_minimal(base_size = 13)

  print(p_interval)
  print(p_cum)

  list(
    interval    = p_interval,
    cumulative  = p_cum,
    pred_counts = pred_counts,
    summary     = interval_summary,
    cum_summary = cum_summary,
    prob        = prob
  )
}
