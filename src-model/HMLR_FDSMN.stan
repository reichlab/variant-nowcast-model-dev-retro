functions {

  vector map_rect_likelihood(vector phi,vector params, real[] data_r, int[] data_i) {
    int K = data_i[1];  // Number of clades
    int B = data_i[2];  // Degrees of freedom for spline
    int N = data_i[3];  // Number of days
    int y[N, K];
    {
      int idx = 1;
      for (n in 1:N) {
        for (k in 1:K) {
          y[n, k] = data_i[3 + idx];
          idx += 1;
        }
      }
    }

    // Extract shared data
    matrix[B, N] x = to_matrix(data_r[1:(B * N)], B, N);
    vector[K-1] aloc = phi[1:K-1];
    vector[K-1] bloc = phi[K:(2*K-2)];
    real asd = phi[2*K -1];
    real bsd = phi[2*K];

    // Extract location-specific parameters
    vector[K-1] alpha_nc = segment(params, 1, K-1);
    matrix[K-1, B] beta_nc = to_matrix(segment(params, K, (K-1)*B), K-1, B);

    // Compute raw alpha and beta
    vector[K] alpha = append_row(0, aloc + asd * alpha_nc);
    matrix[K, B] beta;
    for (b in 1:B) {
      beta[:, b] = append_row(0, bloc + bsd * beta_nc[:, b]);
    }

    // Calculate log-likelihood
    real acc = 0;
    for (n in 1:N) {
      if (sum(y[n, :]) > 0) {
        vector[K] lambda = exp(alpha + beta * x[:, n]) / sum(exp(alpha + beta * x[:, n]));
        acc += multinomial_lpmf(y[n, :] | lambda);
      }
    }

    return [acc]'; // Return log-likelihood as a vector
  }
}

data {
  int<lower=0> N; // Number of days of samples
  int<lower=1> L; // Number of locations
  int<lower=1> K; // Number of clades
  int<lower=3> B; // Degrees of freedom for the spline
  int<lower=0> y[N, L, K]; // Count data
  matrix[B, N] x; // Spline matrix
}

transformed data{
 // Prepare shared data indices
  int shared_data_i[L, 4 + N * K -1];
  real shared_data_r[L,(B * N)];
  for (l in 1:L) {
  shared_data_r[l] = to_array_1d(x);
  shared_data_i[l,1] = K;
  shared_data_i[l,2] = B;
  shared_data_i[l, 3] = N;
  shared_data_i[l,(4):(4 + (N * K)-1)] = to_array_1d(y[:, l, :]);
  }
}
parameters {
  real<lower=0> bsd; // Prior sd for betas
  real<lower=0> asd; // Prior sd for the alphas
  real bloc[K-1]; // Prior means for betas
  real aloc[K-1]; // Prior means for the alphas
  vector[K-1] alpha_noncentered[L]; // Non-centered alpha parameters
  matrix[K-1, B] beta_noncentered[L]; // Non-centered beta parameters
}

model {
  // Priors
  bsd ~ normal(1, 400);
  asd ~ normal(1, 400);
  bloc ~ normal(0, 400);
  aloc ~ normal(0, 400);

  for (l in 1:L) {
    alpha_noncentered[l] ~ normal(0, 1);
    to_vector(beta_noncentered[l]) ~ normal(0, 1);
  }

    // Vectorize parameters for map_rect
  vector[(K-1) + (K-1)*B] param_shards[L];
  for (l in 1:L) {
    param_shards[l] = append_row(alpha_noncentered[l], to_vector(beta_noncentered[l]));
  }

  // Prepare shared data
    vector[K-1] aloc_vec = to_vector(aloc);
    vector[K-1] bloc_vec = to_vector(bloc);
    vector[2] scalar_vec = [asd, bsd]';
    vector[2*K] phi = append_row(aloc_vec, append_row(bloc_vec, scalar_vec));

  // Call map_rect with the correct arguments
  target += sum(map_rect(map_rect_likelihood,phi, param_shards, shared_data_r, shared_data_i));

}
generated quantities {
  // Declare raw_alpha and raw_beta
 vector[K-1] raw_alpha[L];
  matrix[K-1, B] raw_beta[L];
  for (l in 1:L) {
    for (k in 1:(K-1)) {
      raw_alpha[l, k] = aloc[k] + asd * alpha_noncentered[l, k]; // Reparameterized alpha
      for (b in 1:B) {
        raw_beta[l, k, b] = bloc[k] + bsd * beta_noncentered[l, k, b]; // Reparameterized beta
      }
    }
  }
}
