ebola_drc <- read.csv("ebola_drc_upto23Aug.csv")

print(ebola_drc)

N_obs = 23000
Tfinal <- 101 ## The final Time

rho_mean = 1 / N_obs
rho_sd = 0.001

beta_params(rho_mean, rho_sd)

p_mean = 0.35 ## Using CDC result, under-reporting between 30-40%
p_sd = 0.01

beta_params(p_mean, p_sd)

iter = 1e4
final_data = matrix(c(ebola_drc$Day, ebola_drc$moving_avg_3), ncol = 2)

DRC = LMGP(
  data_obs = final_data,
  Tmax = Tfinal,
  N = N_obs,
  prob = TRUE,
  r_a = 0.002,
  r_b = 43,
  p_a = 795,
  p_b = 1478,
  iteration = iter
)

zhat_stability_replicate(DRC)
traceplot(DRC$stan_fit, par = c("b", "g", "r", "R0", "p"))

wnorm_DRC = wnorm_calc(DRC$Z_hat, M = length(DRC$Z_hat))

out <- LMGP_ppc(
  prob = TRUE,
  fit_full = DRC$stan_fit,
  data_obs = final_data,
  # cols: time, 3-day MA count
  N = N_obs,
  raw_data = ebola_drc_new$Incidence,
  # optional: pre-MA counts, same length/order
  w_norm = wnorm_DRC$w_norm
)
print(out$interval)
print(out$cumulative)
