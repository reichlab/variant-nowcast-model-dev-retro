// stan file for creating the HMLR_FDLMN model


data {
  int<lower=0> N; // the number of samples
  int<lower=1> L; // number of locations
  int<lower=1> K; // number of clades
  array[N] int<lower=1, upper=K> y; // the clades
  array[N] int<lower=1, upper=L> ll; // locations
  array[N] real x; // the days of the samples
  row_vector [N]  weights; // number of seq per location
}
parameters {
  real<lower=0> bsd; // prior sd for betas
  array[K-1] real bloc; // prior means for betas
  array[L] vector[K-1] raw_alpha; // alpha's without the reference variant
  array[L] vector[K-1] raw_beta; // betas without the reference variant
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
    array[L] vector[K] alpha;
    array[L] vector[K] beta;
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
