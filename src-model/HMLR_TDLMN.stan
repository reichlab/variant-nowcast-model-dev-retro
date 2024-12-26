data {
  int<lower=0> N; // the number of samples
  int<lower=1> L; // number of locations
  int<lower=1> K; // number of clades
  int<lower=1, upper=K> y[N]; // the clades
  int<lower=1, upper=L> ll[N]; // the locations
  real x[N]; // the days of the samples
  row_vector[N] weights; // number of seq per location
}
parameters {
  vector<lower=0>[K-1] bsd; // prior sd for betas
  vector<lower=0>[K-1] asd; // prior sd for the alphas
  real bloc[K-1]; // prior means for betas
  real aloc[K-1]; // prior means for the alphas
  vector[K-1] alpha_noncentered[L]; // non-centered alpha parameters
  vector[K-1] beta_noncentered[L]; // non-centered beta parameters
}
transformed parameters {
  // Centered alpha and beta parameters

  // Calculated alpha_raw and beta_raw
  vector[K-1] raw_alpha[L];
  vector[K-1] raw_beta[L];

  // Apply non-centered transformation and back-calculate raw values
  for (l in 1:L) {
    for( k in 1:(K-1)){
      raw_alpha[l, k] = aloc[k] + asd[k]*alpha_noncentered[l,k];
      raw_beta[l,k] = bloc[k] + bsd[k]*beta_noncentered[l,k];
    }
  }
}
model {
  // Priors for hyperparameters
  bsd ~ normal(1, 400); // prior for beta standard deviation
  asd ~ normal(1, 400); // prior for alpha standard deviation
  bloc ~ normal(0, 400); // prior for beta mean
  aloc ~ normal(0, 400); // prior for alpha mean

  // Standard normal priors for non-centered parameters
  for (l in 1:L) {
    alpha_noncentered[l] ~ normal(0, 1);
    beta_noncentered[l] ~ normal(0, 1);
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
    target += weights * log_prob;
  }
}
