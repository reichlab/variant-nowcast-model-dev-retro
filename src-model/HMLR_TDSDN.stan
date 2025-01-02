functions {
  // function for dirichlet_multinomial_lpmf comes from https://discourse.mc-stan.org/t/transforming-a-multinomial-model-into-a-dirichlet-multinomial/26399/2
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real sum_alpha = sum(alpha);
    return lgamma(sum_alpha) - lgamma(sum(y) + sum_alpha)
           // + lgamma(sum(y)+1) - sum(lgamma(to_vector(y)+1) // constant, may omit
           + sum(lgamma(to_vector(y) + alpha)) - sum(lgamma(alpha));
  }
}

data {
  int<lower=0> N; // the number of days of samples fit to
  int<lower=1> L; // number of locations
  int<lower=1> K; // number of clades
  int<lower=3> B; // the degrees of freedom
  int<lower=0> y[N,L,K]; // the vectors of counts for each day and location
  matrix[B, N] x; // the spline over the days of the samples
}
parameters {
  vector<lower=0>[K-1] bsd; // prior sd for betas
  vector<lower=0>[K-1] asd; // prior sd for the alphas
  real<lower=0> kappa; // the scale for the dirchlet
  real bloc[K-1]; // prior means for betas
  real aloc[K-1]; // prior means for the alphas
  vector[K-1] alpha_noncentered[L]; // non-centered alpha parameters
  matrix[K-1,B] beta_noncentered[L]; // non-centered beta parameters
}
transformed parameters {
  // Centered alpha and beta parameters

  // Calculated alpha_raw and beta_raw
  vector[K-1] raw_alpha[L];
  matrix[K-1, B] raw_beta[L];

  // Apply non-centered transformation and back-calculate raw values
  for (l in 1:L) {
    for( k in 1:(K-1)){
      raw_alpha[l, k] = aloc[k] + asd[k]*alpha_noncentered[l,k];
      for(b in 1:B){
       raw_beta[l,k, b] = bloc[k] + bsd[k]*beta_noncentered[l,k, b];
      }
    }
  }
}
model {
  // Priors for hyperparameters
  bsd ~ normal(1, 400); // prior for beta standard deviation
  asd ~ normal(1, 400); // prior for alpha standard deviation
  bloc ~ normal(0, 400); // prior for beta mean
  aloc ~ normal(0, 400); // prior for alpha mean
  kappa ~ normal(1, 5); // prior for scale

  // Standard normal priors for non-centered parameters
  for (l in 1:L) {
    alpha_noncentered[l] ~ normal(0, 1);
    for (k in 1:(K-1)) {
      beta_noncentered[l, k] ~ normal(0, 1); // standard normal for beta_nc
    }
  }

  {
    // Temporary variables to avoid repeated calculation
    vector[K] alpha[L];
    matrix[K,B] beta[L];
    for (l in 1:L) {
      alpha[l] = append_row(0, raw_alpha[l]);
      for (b in 1:B) {
        beta[l, :, b] = append_row(0, raw_beta[l, :, b]);
      }
    }
    // running the model for each location
    for(l in 1:L){
     for (n in 1:N) {
      if(sum(y[n,l,:]) == 0){
         continue;
       } else{
        y[n, l, :] ~ dirichlet_multinomial(kappa*exp(alpha[l] + beta[l, :, :]* x[:, n])/sum(exp(alpha[l] + beta[l, :, :]* x[:, n])));
       }
    }
    }
  }
}
