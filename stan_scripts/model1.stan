data {
 int<lower=1> n;
 vector[n] p_y;
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
 for (i in 1:n)
 {
   target += log_mix(p_y[i],
		     bernoulli_lpmf(1|p[i]),
		     bernoulli_lpmf(0|p[i]));
 }
 // Priors
 gamma0 ~ normal(0,5);
 gamma1 ~ normal(0,5);
 }




