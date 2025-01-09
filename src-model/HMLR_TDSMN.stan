functions {
  real partial_sum_lpmf(int[] slice_y, int start, int end, int[] ll, matrix x,
                        row_vector weights, vector[] raw_alpha, matrix[] raw_beta,
                        int K, int B) {
    real acc = 0;
    for (n in 1:size(slice_y)) {
      vector[K] alpha = append_row(0, raw_alpha[ll[start + n - 1]]);
      matrix[K, B] beta;
      for (b in 1:B) {
        beta[:, b] = append_row(0, raw_beta[ll[start + n - 1], :, b]);
      }
      acc += weights[start + n - 1] * categorical_logit_lpmf(slice_y[n] | alpha + beta * x[:, start + n - 1]);
    }
    return acc;
  }
}

data {
  int<lower=0> N; // the number of samples
  int<lower=1> L; // number of locations
  int<lower=1> K; // number of clades
  int<lower=1> B; // the number of df
  array[N] int<lower=1, upper=K> y; // the clades
  array[N] int<lower=1, upper=L> ll; // locations
  matrix[B, N] x; // the spline over the days of the samples
  row_vector[N] weights; // number of seq per location
}

parameters {
  vector<lower=0>[K-1] bsd; // prior sd for betas
  vector<lower=0>[K-1] asd; // prior sd for alphas
  array[K-1] real aloc; // prior means for the alpha's
  array[K-1] real bloc; // prior means for betas

  // Non-centered parameterization
  vector[K-1] alpha_nc[L]; // raw alpha's
  matrix[K-1, B] beta_nc[L]; // raw betas
}

transformed parameters {
  // Apply non-centered parameterization
  vector[K-1] raw_alpha[L];
  matrix[K-1, B] raw_beta[L];
  for (l in 1:L) {
    for (k in 1:(K-1)) {
      raw_alpha[l, k] = aloc[k] + asd[k] * alpha_nc[l, k]; // Reparameterized alpha
      for (b in 1:B) {
        raw_beta[l, k, b] = bloc[k] + bsd[k] * beta_nc[l, k, b]; // Reparameterized beta
      }
    }
  }
}

model {
  // Priors
  bsd ~ normal(1, 400); // prior for the sd
  bloc ~ normal(0, 400); // priors for the betas
  asd ~ normal(1, 400); // prior for the alpha sd
  aloc ~ normal(0, 400); // prior for the alpha means
  
  // Priors on raw parameters (standard normal)
  for (l in 1:L) {
    alpha_nc[l] ~ normal(0, 1); // standard normal for alpha_nc
    for (k in 1:(K-1)) {
      beta_nc[l, k] ~ normal(0, 1); // standard normal for beta_nc
    }
  }
  
  // Likelihood using reduce_sum
  target += reduce_sum(partial_sum_lpmf, y, 50, ll, x, weights, raw_alpha, raw_beta, K, B);
}