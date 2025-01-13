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
  int<lower=0> y[N,L,K]; // the vectors of counts for each day and location
}
parameters {
  real<lower=0> bsd; // prior sd for betas
  real<lower=0> asd; // prior sd for the alphas
  real<lower=0> kappa; // the scale for the dirchlet
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
      raw_alpha[l, k] = aloc[k] + asd*alpha_noncentered[l,k];
      raw_beta[l,k] = bloc[k] + bsd*beta_noncentered[l,k];
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
    // running the model for each location
    for(l in 1:L){
     for (n in 1:N) {
      if(sum(y[n,l,:]) == 0){
         continue;
       } else{
        y[n, l, :] ~ dirichlet_multinomial(kappa*exp(alpha[l] + beta[l]*n)/sum(exp(alpha[l] + beta[l]*n)));
       }
    }
    }
  }
}
