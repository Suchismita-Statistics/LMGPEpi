#' LMGP
#'
#' This function makes posterior inference using LMGP method about the model parameters - infection rate \eqn{\beta}, recovery rate \eqn{\gamma}, limiting proportion of initial infected and susceptibles \eqn{\rho}. For probability model, it additionally estimates under-reporting probability \eqn{p}.
#'
#' @import rstan
#'
#' @param data_obs A two-column matrix with first column days and second column is corresponding counts.
#' @param Tmax Numeric. Number of observation in the data.
#' @param N Initial number of susceptible.
#' @param prob If TRUE fits a model with underreporting. The default is FALSE.
#' @param dt Numeric. time increment. Default is \code{0.1}.
#' @param n_mc The number of samples to be used for Monte Carlo estimates of simplex probability.
#' @param r_a,r_b The parameters of Beta prior of \eqn{\rho}. Default values are 1.
#' @param p_a,p_b The parameters of Beta prior of \eqn{p}. Default values are 1.
#' @param iteration Total number of posterior samples. The default is \eqn{10^4}.
#' @param nchain Integer specifying the number of Markov chains. The default is 1.
#' @param ncore Same as "cores" in Sampling function of Stan. Default is 1.
#' @param stan_seed Same as "seed" in Sampling function of Stan. Default is 1234.
#' @param initial_value description
#' @param initial_value List of length \code{nchain}, one init list per chain.
#'   When \code{NULL} (the default) the same starting point is used for every
#'   chain: \code{b = 1}, \code{g = 0.1}, \code{r = 0.01} and
#'   \code{u = rep(1 / (Tmax + 1), Tmax + 1)} (the uniform simplex), with
#'   \code{p = 0.5} added when \code{prob = TRUE}. Because all chains then start
#'   from an identical point, between-chain diagnostics such as R-hat are
#'   optimistic; supply dispersed inits when convergence is being assessed.
#' @param single_run_adj Logical, controlling how far the analysis is carried and
#'   hence what is returned. Sampling and the \eqn{Z(\theta)} computation are
#'   identical either way.
#'
#'   \code{FALSE} (the default) stops after \code{Z_hat}, returning the fit object
#'   and the per-draw weights untouched. This is the mode for simulation studies:
#'   run many replicates, then let \code{\link{summary_stan_list_adj}} rebuild the
#'   weights from \code{Z_hat} and pool across them. Use it whenever the posterior
#'   draws are needed afterwards.
#'
#'   \code{TRUE} completes the analysis for a single dataset. The importance
#'   weights are normalised by \code{\link{wnorm_calc}} and, for each parameter,
#'   the weight-corrected posterior mean, the bias-corrected weighted standard
#'   deviation and the weighted 2.5\% and 97.5\% quantiles are returned. Two
#'   efficiency measures are reported separately and should not be merged:
#'   \code{ess_stan} and \code{ess_stan_per_sec} are Stan's bulk-ESS relative to
#'   post-warmup sampling time alone, whereas \code{ess_combined} and
#'   \code{ess_comb_per_sec} scale that bulk-ESS by the weight-degeneracy factor
#'   (Kish ESS divided by the number of draws) and divide by the total time
#'   including the \eqn{Z(\theta)} loop, giving the throughput of the corrected
#'   pipeline as a whole. Note that this branch does \strong{not} return
#'   \code{fit_full} or \code{Z_hat}; the draws are discarded.
#' @param tn Logical. When \code{TRUE}, \eqn{Z(\theta)} is computed with
#'   \code{\link{compute_Z_tn}} rather than \code{\link{compute_Z_mc}}.
#'
#' @return A list. When \code{single_run_adj = FALSE}: \code{fit_full} (the
#'   \code{stanfit}), \code{Z_hat}, \code{total_time} and \code{Tmax}. When
#'   \code{single_run_adj = TRUE}: \code{Tmax}, \code{weighted_mean},
#'   \code{weighted_sd}, \code{q2.5}, \code{q97.5}, \code{ess_stan},
#'   \code{ess_combined}, \code{ess_stan_per_sec} and \code{ess_comb_per_sec},
#'   each over \code{b}, \code{g}, \code{r}, \code{R0} and, when
#'   \code{prob = TRUE}, \code{p}.
#'
#' @section Dependencies:
#' Requires the compiled Stan models \code{sm_eff} and \code{sm_eff_prob}, and the
#' helpers \code{compute_Z_mc}, and \code{wnorm_calc}, to be
#' available in the calling environment.
#'
#' @export

LMGP <- function(data_obs, Tmax, N, prob = FALSE, dt = 0.1, n_mc = 2000, r_a = 1, r_b = 1, p_a = 1, p_b = 1,
                     iteration = 1e4,  ncore = 1, nchain = 1, stan_seed = 1234, initial_value = NULL, single_run_adj = TRUE) {

  obs_len = nrow(data_obs)

  if(is.null(initial_value))
  {
    init_one <- list(b = 1, g = 0.1, r = 0.01, u = rep(1 / ( obs_len  + 1),  obs_len + 1))
    if(prob) init_one$p <- 0.5
    initial_value <- replicate(nchain, init_one, simplify = FALSE)
  }

  K <- sum(data_obs[, 2])

  int_len <- length(seq(dt, 1, dt))
  max_len <- int_len*Tmax+ 1
  Time_mat <- matrix(0, nrow = obs_len, ncol = max_len)
  for (i in 1:obs_len) {
    Time_mat[i, 1:(i * int_len)] <- seq(dt, i, dt)
  }

  A <- matrix(0, obs_len, obs_len)
  diag(A) <- rep(1, obs_len)
  for (i in 1:(obs_len - 1)) {
    A[i, i + 1] <- -1
  }

  obs_t <- data_obs[, 1]
  full_grid <- seq(0, max(obs_t), dt)
  inx <- round(obs_t / dt) + 1
  # halve dt until every observation time is a grid point an even number of
  # steps from 0 (composite Simpson needs an even count of sub-intervals)
  dt_in <- dt
  while (max(abs(full_grid[inx] - obs_t)) > 1e-8 || any(inx %% 2 == 0)) {
    dt <- dt / 2
    if (dt < 1e-6) stop("cannot place data_obs[, 1] on a uniform grid")
    full_grid <- seq(0, max(obs_t), dt)
    inx <- round(obs_t / dt) + 1
  }
  if (dt < dt_in) message(sprintf("dt reduced %g -> %g (%d grid points)", dt_in, dt, length(full_grid)))

  len   <- rep(c(4, 2), length(full_grid) / 2)
  coeff <- c(1, len[-length(len)], 1)

  data_list <- list( N = N, T_max_int = obs_len,
                     #infection_count = if(prob) c(data_obs[, 2], N - K) else data_obs[, 2],
                     infection_count =  c(data_obs[, 2], N - K),
                     t0 = 0, dt = dt,
                     full_len = length(full_grid),
                     simpson_coeff = coeff, full_grid = full_grid[full_grid > 0],
                     index = inx, A = A
  )
  if(prob) data_list <- c(data_list, list(r_a = r_a, r_b = r_b, p_a = p_a, p_b = p_b))

  # if (prob == TRUE) {
  #   stan_file = system.file("stan", "LMGP_with_prob.stan", package = "LMGPEpi")
  # } else {
  #   stan_file = system.file("stan", "LMGP_wo_prob.stan", package = "LMGPEpi")
  # }
  #
  # sm = rstan::stan_model(file = stan_file)

  stan_name <- if (isTRUE(prob)) {
    "LMGP_with_prob.stan"
  } else {
    "LMGP_wo_prob.stan"
  }

  stan_file <- system.file(
    "stan",
    stan_name,
    package = "LMGPEpi"
  )

  if (!nzchar(stan_file) || !file.exists(stan_file)) {
    stop(
      "Could not find Stan file: ", stan_name,
      "\nExpected path: ", stan_file
    )
  }

  message("Using Stan file: ", normalizePath(stan_file))

  sm <- rstan::stan_model(file = stan_file)

  # -------------------------------------------------------------
  # Clock 1: Stan sampling (includes warmup)
  # -------------------------------------------------------------
  t_stan <- system.time({
    fit_full <- rstan::sampling(
      sm, data = data_list, iter = iteration, chain = nchain, cores = ncore,
      init = initial_value, seed = stan_seed
    )
  })["elapsed"]

  # -------------------------------------------------------------
  # Z(theta) per draw -> importance weights
  # -------------------------------------------------------------
  mc_draws <- rstan::extract(fit_full, pars = c("mn", "cov"))
  M <- dim(mc_draws$mn)[1]

  Z_hat <- numeric(M)

  t_correction <- system.time({
    for (m in seq_len(M)) {
        out <- compute_Z_mc(mc_draws$mn[m, ], mc_draws$cov[m, , ], n_mc = n_mc)
      Z_hat[m] <- out$Z
    }
  })["elapsed"]

  weights_calc <- wnorm_calc(Z_hat, M)
  ok <- weights_calc$ok
  w_norm <- weights_calc$w_norm

  total_time       <- t_stan + t_correction

  if(single_run_adj)
  {
    # -------------------------------------------------------------
    # ESS/s -- two numbers, two clocks. Do not merge them.
    # -------------------------------------------------------------

    et <- rstan::get_elapsed_time(fit_full)
    sample_time <-  sum(et[, "sample"])
    draws_arr <- posterior::as_draws_array(fit_full)
    ess_stan  <- posterior::summarise_draws(draws_arr, ess_bulk = posterior::ess_bulk)
    ess_stan$ess_per_sec <- ess_stan$ess_bulk / sample_time

    ess_kish         <- 1 / sum(w_norm^2, na.rm = TRUE)


    ess_combined <- ess_kish*(ess_stan$ess_bulk)/length(w_norm)
    ess_comp_ber_sec <- ess_combined/total_time


    # ---------------------------------------------------------
    # Correction Using Weights
    # ---------------------------------------------------------

    param_names <- if(prob) c("b", "g", "r", "R0", "p") else c("b", "g", "r", "R0")
    ext = rstan::extract(fit_full, pars = param_names)

      keep <- match(param_names, ess_stan$variable)

    wtd_mean <- sapply(param_names, function(a) sum(w_norm * ext[[a]]))
    wtd_sd   <- sapply(param_names, function(a) {
      m <- wtd_mean[[a]]
      sqrt(sum(w_norm * (ext[[a]] - m)^2) / (1 - sum(w_norm^2)))
    })
    q_2.5  <- sapply(param_names, function(a) Hmisc::wtd.quantile(ext[[a]], w_norm, probs = 0.025, normwt = TRUE))
    q_97.5 <- sapply(param_names, function(a) Hmisc::wtd.quantile(ext[[a]], w_norm, probs = 0.975, normwt = TRUE))

    summary_table <- cbind("Weighted mean" = as.numeric(wtd_mean),
                           "Weighted SD"   = as.numeric(wtd_sd),
                           "2.5%"          = as.numeric(q_2.5),
                           "97.5%"         = as.numeric(q_97.5),
                           "ESS_comb/s"    = as.numeric(ess_comp_ber_sec[keep]))
    rownames(summary_table) <- if(prob) c("beta", "gamma", "rho", "R_0", "reporting_prob") else c("beta", "gamma", "rho", "R_0")
    print(round(summary_table, 3))

    list(
      stan_fit = fit_full,
      Z_hat = Z_hat,
      weighted_mean    = wtd_mean,
      weighted_sd      = wtd_sd,
      q2.5             = q_2.5,
      q97.5            = q_97.5,
      ess_stan         = setNames(ess_stan$ess_bulk[keep],      param_names),  # Stan-only bulk-ESS
      ess_combined     = setNames(ess_combined[keep],           param_names),  # MCMC x importance-weight ESS
      ess_stan_per_sec = setNames(ess_stan$ess_per_sec[keep],   param_names),  # Stan-only throughput
      ess_comp_ber_sec = setNames(ess_comp_ber_sec[keep],       param_names)   # full-pipeline throughput
    )
  }else{
    list(
      fit_full = fit_full, Z_hat = Z_hat, total_time = unname(total_time), Tmax = Tmax
    )
  }
}
