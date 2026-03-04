functions {
  // function for dirichlet_multinomial_lpmf comes from https://discourse.mc-stan.org/t/transforming-a-multinomial-model-into-a-dirichlet-multinomial/26399/2
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real sum_alpha = sum(alpha);
    return lgamma(sum_alpha) - lgamma(sum(y) + sum_alpha)
           + sum(lgamma(to_vector(y) + alpha)) - sum(lgamma(alpha));
  }
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
    vector[K-1] asd = phi[(2*K -1):(3*K-3)];
    vector[K-1] bsd = phi[(3*K-2):(4*K-4)];
    real kappa = phi[4*K-3];
    matrix[K -1, K-1] Omega_alpha = to_matrix(phi[(4*K-2):(4*K- 3 + (K-1)*(K-1))], (K-1), (K-1));
    matrix[K -1, K-1] Omega_beta = to_matrix(phi[(4*K- 2 + (K-1)*(K-1)):(4*K- 3 + 2*(K-1)*(K-1))], (K-1), (K-1));
    // Extract location-specific parameters
    vector[K-1] alpha_nc = segment(params, 1, K-1);
    matrix[K-1, B] beta_nc = to_matrix(segment(params, K, (K-1)*B), K-1, B);

    // Compute raw alpha and beta and Sigmas
    vector[K-1] raw_alpha;
  matrix[K-1, B] raw_beta;
  matrix[K-1, B] mat;
   matrix[K-1, K-1] Sigma_alpha = diag_pre_multiply(asd, Omega_alpha);
  matrix[K-1, K-1] Sigma_beta = diag_pre_multiply(bsd, Omega_beta);
    raw_alpha = aloc + Sigma_alpha*alpha_nc;
    mat= Sigma_beta*beta_nc;
    for(b in 1:B){
      raw_beta[:, b] = bloc + mat[:, b];
    }
    vector[K] alpha = append_row(0, raw_alpha);
    matrix[K, B] beta;
    for (b in 1:B) {
      beta[:, b] = append_row(0, raw_beta[:, b]);
    }

    // Calculate log-likelihood
    real acc = 0;
    for (n in 1:N) {
      if (sum(y[n, :]) > 0) {
        vector[K] lambda = kappa*exp(alpha + beta * x[:, n]) / sum(exp(alpha + beta * x[:, n]));
        acc += dirichlet_multinomial_lpmf(y[n, :] | lambda);
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
  vector<lower=0>[K-1] bsd; // Prior sd for betas
  vector<lower=0>[K-1] asd; // Prior sd for the alphas
  real<lower=0> kappa; // Scale for the Dirichlet
  cholesky_factor_corr[K -1] Omega_alpha;
  cholesky_factor_corr[K -1] Omega_beta;
  vector[K-1] bloc; // Prior means for betas
  vector[K-1] aloc; // Prior means for the alphas
  vector[K-1] alpha_noncentered[L]; // Non-centered alpha parameters
  matrix[K-1, B] beta_noncentered[L]; // Non-centered beta parameters
}

model {
  // Priors
  Omega_alpha ~ lkj_corr_cholesky(2);
  Omega_beta ~ lkj_corr_cholesky(2);
  bsd ~ normal(1, 400);
  asd ~ normal(1, 400);
  bloc ~ normal(0, 400);
  aloc ~ normal(0, 400);
  kappa ~ normal(1, 5);

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
    vector[K-1] asd_vec = to_vector(asd);
    vector[K-1] bsd_vec = to_vector(bsd);
    vector[1] kappa_vec = [kappa]';
    vector[(K-1)*(K-1)] omega_a_vec = to_vector(Omega_alpha);
    vector[(K-1)*(K-1)] omega_b_vec = to_vector(Omega_beta);
    vector[4*K -3 + 2*(K-1)*(K-1)] phi = append_row(aloc_vec, append_row(bloc_vec, append_row(asd_vec,append_row(bsd_vec, append_row(kappa_vec, append_row(omega_a_vec, omega_b_vec))))));

  target += sum(map_rect(map_rect_likelihood,phi, param_shards, shared_data_r, shared_data_i));

}
generated quantities {
  // Declare raw_alpha and raw_beta
 vector[K-1] raw_alpha[L];
  matrix[K-1, B] raw_beta[L];
  matrix[K-1, K-1] Sigma_alpha = diag_pre_multiply(asd, Omega_alpha);
  matrix[K-1, K-1] Sigma_beta = diag_pre_multiply(bsd, Omega_beta);
  matrix[K-1, B] mat;
  for (l in 1:L) {
      raw_alpha[l, :] = aloc + Sigma_alpha * alpha_noncentered[l,:]; // Reparameterized alpha
      mat = Sigma_beta*beta_noncentered[l,:];
      for (b in 1:B) {
        raw_beta[l, :, b] = bloc + mat[, b]; // Reparameterized beta
      }
    }
  }
