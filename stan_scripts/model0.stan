data {
 int<lower=1> n;
 int<lower=0,upper=1> y[n];
 vector[n] x;
}

parameters {
 real gamma0;
 real gamma1;
}

model {
 vector[n] p;
 p = inv_logit(gamma0+gamma1*x);
 // likelihood
 y ~ bernoulli(p);
 // Priors
 gamma0 ~ normal(0,5);
 gamma1 ~ normal(0,5);
}




