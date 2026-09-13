library(LMGPEpi)

# The true parameters are N = 1e4; beta = 1.5; gamma = 1; rho = 0.05; p = 0.6; Tmax = 10


trial_data = noisy_data_sim(N = 1e4, beta = 1.5, gamma = 1, rho = 0.05, p = 0.6, Tmax = 10) #Substitute with your Data
print(trial_data)


#-----------------------------------------------------------------------
# Fitting using true counts
#-----------------------------------------------------------------------

## Substitute according to your model

Tfinal = 10 ## The final observation period of the epidemic
N_obs = 1e4 ## Initial susceptible


## Fitting LMGP Model with true counts without under-reporting

tot_cores = parallel::detectCores()
num_cores = min(tot_cores, 2) ## You can specify the number of cores/chains you need

LMGP_true_cts = LMGP(data_obs = trial_data, Tmax = Tfinal, N = N_obs, prob = FALSE,
                     iteration = 5e3,  ncore = num_cores, nchain = num_cores)


## Check the simplex probability, and the performance of the weights
zhat_stability_replicate(LMGP_true_cts) ## Comments: If the mean is small, probably the
## final time/initial susceptible are incorrect. Try changing the value.



## Check the convergence
traceplot(LMGP_true_cts$stan_fit, par = c("b", "g", "r", "R0"))


wnorm_true_cts = wnorm_calc(LMGP_true_cts$Z_hat, M = length(LMGP_true_cts$Z_hat))

out <- LMGP_ppc(
  prob = FALSE,
  fit_full = LMGP_true_cts$stan_fit,
  data_obs = trial_data,     # cols: time, 3-day MA count
  N = N_obs,
  raw_data = LMGP_true_cts$Incidence,    # optional: pre-MA counts, same length/order
  w_norm = wnorm_true_cts$w_norm
)
print(out$interval)
print(out$cumulative)


#--------------------------------------------------------------------------------------------
## Fitting LMGP Model with under-reported counts using under-reporting model
#--------------------------------------------------------------------------------------------

# For the noisy model, we have to give a strong prior for rho parameter, assuming we have true information for this.
# True rho is 0.05. So, we find the parameters of Beta distribution, which has mean 0.05, and sd 0.001.
beta_prior_params = round(beta_params(0.05, 0.001))

cat("Therefore, r_a = ", beta_prior_params[1], " r_b = ", beta_prior_params[2], "for fitting LMGP adjusted by under-reporting probability.")

LMGP_adj = LMGP(data_obs = trial_data[, c(1, 3)], Tmax = Tfinal, N = N_obs, prob = TRUE,
                r_a = beta_prior_params[1], r_b = beta_prior_params[2], iteration = 5e3)

zhat_stability_replicate(LMGP_adj)
traceplot(LMGP_adj$stan_fit, par = c("b", "g", "r", "R0", "p"))


## PPC

wnorm_adj = wnorm_calc(LMGP_adj$Z_hat, M = length(LMGP_adj$Z_hat))

out <- LMGP_ppc(
  prob = TRUE,
  fit_full = LMGP_adj$stan_fit,
  data_obs = trial_data[, c(1, 3)],
  N = N_obs,
  raw_data = LMGP_adj$Incidence,
  w_norm = wnorm_adj$w_norm
)
print(out$interval)
print(out$cumulative)


#--------------------------------------------------------------------------------------------
## To illustrate the use of summary function
#--------------------------------------------------------------------------------------------

count_data_list = list()

for(i in 1:3)
{
  cat(i, "th simulation is running")
  count_data_list[[i]] = LMGP(data_obs = trial_data, Tmax = Tfinal, N = N_obs, prob = FALSE,
                              iteration = 1e3, single_run_adj = FALSE)
}

summary_stan_list_adj(result_list = count_data_list, true_value = c(1.5, 1, 0.05))
