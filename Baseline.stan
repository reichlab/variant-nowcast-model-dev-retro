data {
  int<lower=0> N; // the number of samples
  int<lower=1> K; // number of clades
  array[N] int<lower=1, upper=K> y; // the clades
  array[N] real x; // the days of the samples
  row_vector [N]  weights; // number of seq per location
}
parameters {
  vector[K-1] raw_alpha; // alpha's without the reference variant
  vector[K-1] raw_beta; // betas without the reference variant
}
model {
    // Temporary variables to avoid repeated calculation
    vector[K] alpha;
    vector[K] beta;
    alpha = append_row(raw_alpha, 0);
    beta = append_row(raw_beta, 0);

    // Vectorized likelihood calculation
    vector[N] log_prob;
    for (n in 1:N) {
      log_prob[n] = categorical_logit_lpmf(y[n] | alpha + beta* x[n]);
    }
    target += weights*log_prob;
}
