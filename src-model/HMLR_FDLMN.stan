// stan file for creating the HMLR_FDLMN model


data {
  int<lower=0> N; // the number of samples
  int<lower=1> L; // number of locations
  int<lower=1> K; // number of clades
  int<lower=1, upper=K> y[N]; // the clades
  int<lower=1, upper=L> ll[N]; // locations
  real x[N]; // the days of the samples
  row_vector [N]  weights; // number of seq per location
}
parameters {
  real<lower=0> bsd; // prior sd for betas
  real bloc[K-1]; // prior means for betas
  vector[K-1] raw_alpha[L]; // alpha's without the reference variant
  vector[K-1] raw_beta[L]; // betas without the reference variant
}
model {
  bsd ~ normal(1, 0.1); // prior for the sd
  bloc ~ normal(0, 0.2); // priors for the betas and alphas

  for (k in 1:(K-1)) {
    raw_beta[:, k] ~ normal(bloc[k], bsd);
    raw_alpha[:, k] ~ normal(0, 6);
  }

  {
    // Temporary variables to avoid repeated calculation
    vector[K] alpha[L];
    vector[K] beta[L];
    for (l in 1:L) {
      alpha[l] = append_row(0, raw_alpha[l]);
      beta[l] = append_row(0, raw_beta[l]);
    }

    // Vectorized likelihood calculation
    vector[N] log_prob;
    for (n in 1:N) {
      log_prob[n] = categorical_logit_lpmf(y[n] | alpha[ll[n]] + beta[ll[n]] * x[n]);
    }
    target += weights*log_prob;
  }
}
