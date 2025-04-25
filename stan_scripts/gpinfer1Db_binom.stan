data {
  int<lower=1> N1; // Observed Locations
  int<lower=1> N2; // Predicted Locations
  int<lower=1> D; // Number of Dimensions
  array[N1] vector[D] x1; // Observed Predictor
  array[N2] vector[D] x2; // Extrapolation Predictror
  array[N1] int<lower=0,upper=1> y1; 
}

transformed data {
  real delta = 1e-9;
  int<lower=1> N = N1 + N2;
  array[N] vector[D] x;
  for (n1 in 1:N1) {
    x[n1] = x1[n1];
  }
  for (n2 in 1:N2) {
    x[N1 + n2] = x2[n2];
  }
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> gamma0;
  vector[N] eta;
}

transformed parameters {
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
}

  model {
    rho ~ inv_gamma(5, 5);
    alpha ~ std_normal();
    gamma0 ~ std_normal();
    eta ~ std_normal();
    y1 ~ bernoulli_logit(gamma0+f[1:N1]);
 //   y1 ~ bernoulli_logit(f[1:N1]);
  }
  generated quantities {
    vector[N2] y2;
    for (n2 in 1:N2) {
      y2[n2] = f[N1 + n2];
    }
  }
