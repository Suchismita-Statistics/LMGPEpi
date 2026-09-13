functions{
  vector fund_ode(real t, vector y, array[] real parms){
      real b = parms[1];
      real g = parms[2];
      real r = parms[3];
      vector[6] dydt;
      dydt[1] = -b*y[1]*y[2];
      dydt[2] = b*y[1]*y[2] - g*y[2];
      dydt[3] = b*y[2]*y[3] - b*y[2]*y[5];
      dydt[4] = b*y[2]*y[4] - b*y[2]*y[6];
      dydt[5] = b*y[1]*y[3] - b*y[1]*y[5] + g*y[5];
      dydt[6] = b*y[1]*y[4] - b*y[1]*y[6] + g*y[6];
      return dydt;
      }
}
data{
  int<lower = 0> N;                                 // newly infected total
  int<lower = 0> T_max_int;
  array[T_max_int + 1] int<lower = 0> infection_count;
  real t0;
  real dt;
  int<lower = 0> full_len;
  vector[full_len] simpson_coeff;
  array[full_len - 1] real full_grid;
  array[T_max_int] int<lower = 0> index;
  matrix[T_max_int, T_max_int] A;
}
parameters {
  real<lower = 0> b;
  real<lower = 1e-6> g;
  real<lower = 0, upper = 1> r;
  simplex[T_max_int + 1] u;
}
transformed parameters {
   real<lower = 0> R0 = b / g;
  vector[T_max_int] mn;
  matrix[T_max_int, T_max_int] cov;
  {
    array[3] real parms;
    vector[6] ic;
    array[full_len] vector[6] fwd_sol;
    array[T_max_int, 2] real inv_store;
    vector[full_len] s_full;
    vector[full_len] iota_full;
    real N3;
    vector[4] temp;
    real det;
    matrix[full_len, T_max_int] phi11;
    matrix[full_len, T_max_int] phi12;
    matrix[T_max_int, T_max_int] cov_mat;

    parms[1] = b;
    parms[2] = g;
    parms[3] = r;
    ic[1] = 1.0;
    ic[2] = r;
    ic[3] = 1.0;
    ic[4] = 0;
    ic[5] = 0;
    ic[6] = 1.0;
    N3 = 3 * N;

    fwd_sol[1] = ic;
    fwd_sol[2:full_len] = ode_rk45(fund_ode, ic, t0, full_grid, parms);
    s_full    = to_vector(fwd_sol[, 1]);
    iota_full = to_vector(fwd_sol[, 2]);

    for (i in 1:T_max_int) {
      temp = fwd_sol[index[i], 3:6];
      det = temp[1]*temp[4] - temp[2]*temp[3];
      inv_store[i, 1] =  temp[4] / det;
      inv_store[i, 2] = -temp[3] / det;
      phi11[, i] = to_vector(fwd_sol[, 3]) * inv_store[i, 1]
                 + to_vector(fwd_sol[, 4]) * inv_store[i, 2];
      phi12[, i] = to_vector(fwd_sol[, 5]) * inv_store[i, 1]
                 + to_vector(fwd_sol[, 6]) * inv_store[i, 2];
    }

    for (i in 1:T_max_int) {
      int upto = index[i];
      vector[upto] s_i     = s_full[1:upto];
      vector[upto] iota_i  = iota_full[1:upto];
      vector[upto] phi11_i = phi11[1:upto, i];
      vector[upto] phi12_i = phi12[1:upto, i];
      vector[upto] diff_i  = phi11_i - phi12_i;
      vector[upto] simp_coeff = rep_vector(1, upto);
      simp_coeff[1:(upto - 1)] = simpson_coeff[1:(upto - 1)];

      cov_mat[i, i] =
          sum(dt * b * s_i .* iota_i .* (diff_i .* diff_i) .* simp_coeff) / N3
        + sum(dt * g * phi12_i .* phi12_i .* iota_i           .* simp_coeff) / N3;

      if (i < T_max_int) {
        for (j in (i + 1):T_max_int) {
          vector[upto] diff_j  = phi11[1:upto, j] - phi12[1:upto, j];
          vector[upto] phi12_j = phi12[1:upto, j];
          cov_mat[i, j] =
              sum(dt * b * s_i .* iota_i .* (diff_i .* diff_j) .* simp_coeff) / N3
            + sum(dt * g * phi12_i .* phi12_j .* iota_i        .* simp_coeff) / N3;
          cov_mat[j, i] = cov_mat[i, j];
        }
      }
    }

    mn  = A * s_full[index];
    cov = A * cov_mat * A';

    for (i in 1:T_max_int) {
      for (j in (i + 1):T_max_int) {
        real m = 0.5 * (cov[i, j] + cov[j, i]);
        cov[i, j] = m;
        cov[j, i] = m;
      }
    }
  }
}
model {
  target += gamma_lpdf(b | 0.1, 0.1);
  target += gamma_lpdf(g | 0.1, 0.1);
  target += uniform_lpdf(r | 0, 1);

   u[1:T_max_int] ~ multi_normal(mn, cov);

  target += infection_count[1] * log(1 - sum(u[1:T_max_int]));
  for (i in 1:(T_max_int - 1)) {
    target += infection_count[i + 1] * log(u[i]);
  }
  target += infection_count[T_max_int + 1] * log(u[T_max_int]);
}

