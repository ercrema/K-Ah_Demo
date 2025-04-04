data {
  int<lower=1> N; // Number of observations
  vector[2] coords[N]; // Y coordinates
  int x[N]; // Predictor
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  real gamma0;
  real gamma1;
}

transformed data {
  //matrix[N, N] K = cov_exp_quad(coords, alpha, rho);
  matrix[N, N] K = cov_exp_quad(coords, alpha,rho );
  vector[N] mu = rep_vector(0, N);
  real sq_sigma = square(sigma);
  // diagonal elements
  for (n in 1:N) {
    K[n, n] = K[n, n] + sq_sigma;
  }
}

generated quantities {
  vector[N] lambda;
  vector[N] eta;
  int y[N];
  eta =  multi_normal_rng(mu, K);
  for (i in 1:N)
  {
    lambda[i] = gamma0 + gamma1*x[i] + eta[i];
  }
  y = bernoulli_logit_rng(lambda);
}
