data {
 int<lower=1> N; // Sample size
 int<lower=1> D; // Number of Dimensions
 array[N] vector[D] x; // Predictor (locations)
 array[N] int<lower=0,upper=1> y; // response
 }

transformed data {
  real delta = 1e9;
}

parameters {
 real<lower=0> rho;
 real<lower=0> alpha;
 real gamma0;
 vector[N] eta;
}

model {
  vector[N] f;
  {
    matrix[N, N] L_K;
    matrix[N, N] K = gp_exp_quad_cov(x, alpha, rho);
    // diagonal elements
    for (n in 1:N) {
      K[n, n] = K[n, n] + delta;
    }
    L_K = cholesky_decompose(K);
    f = L_K * eta;
  }
  rho ~ inv_gamma(5, 20);
  alpha ~ std_normal();
  eta ~ std_normal();
  gamma0 ~ std_normal();
  y ~ bernoulli_logit(gamma0+f);
}


