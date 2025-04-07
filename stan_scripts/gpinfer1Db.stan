data {
  int<lower=1> N1;
  array[N1] real x1;
  vector[N1] y1;
  int<lower=1> N2;
  array[N2] real x2;
}
transformed data {
  real delta = 1e-9;
  int<lower=1> N = N1 + N2;
  array[N] real x;
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
  real<lower=0> sigma;
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
    sigma ~ std_normal();
    eta ~ std_normal();
    y1 ~ normal(f[1:N1], sigma);
  }
  generated quantities {
    vector[N2] y2;
    for (n2 in 1:N2) {
      y2[n2] = f[N1 + n2];
    }
  }
