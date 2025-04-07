data {
 int<lower=1> N; // Sample size
 int<lower=1> D; // Number of Dimensions
 array[N] vector[D] x; // Predictor (locations)
 vector[N] y; // response
 }

transformed data {
  vector[N] mu = rep_vector(0, N);
}

parameters {
 real<lower=0> rho;
 real<lower=0> alpha;
 real<lower=0> sigma;
}

model {
  matrix[N, N] L_K;
  matrix[N, N] K = gp_exp_quad_cov(x, alpha, rho);
  real sq_sigma = square(sigma);
  // diagonal elements
  for (n in 1:N) {
    K[n, n] = K[n, n] + sq_sigma;
  }
  L_K = cholesky_decompose(K);
  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  sigma ~ std_normal();
  y ~ multi_normal_cholesky(mu, L_K);
}


