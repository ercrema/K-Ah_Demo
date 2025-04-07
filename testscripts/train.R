library(here)
library(rstan)
library(mgcv)

# 1-D Gaussian Response ----
mcycle <- MASS::mcycle
plot(accel~times,data=mcycle)
x <- mcycle$times
y <- mcycle$accel
D  <- 1

# Just fitting
gp <- stan_model(here('stan_scripts','gpinfer1D.stan'))
fit <- sampling(gp,list(N=length(x),D=1,x=matrix(x,ncol=1),y=y),iter=1000,chains=3,cores=3)


gp2 <- stan_model(here('stan_scripts','gpinfer1Db.stan'))

dat <- list(N1 =  length(x),
	    x1 = x,
	    y1 = y,
	    N2 = length(seq(2,60,1)),
	    x2 = seq(2,60,1))

# Basic predictions
fit2 <- sampling(gp2,dat,iter=500,chains=3,cores=3)
posterior  <-  extract(fit2)

plot(x,y)
lines(seq(2,60,1),apply(posterior$y2,2,mean))

# Fast analytical predictions
gp3 <- stan_model(here('stan_scripts','gpinfer1Dc.stan'))
fit3 <- sampling(gp3,dat,iter=1000,chains=3,cores=3)
posterior  <-  extract(fit3)
plot(x,y)
lines(seq(2,60,1),apply(posterior$f2,2,mean))


# 1-D Binomial Response ----
y.p  <- plogis(scale(y))
plot(x,y.p)
y <- rbinom(length(y),size=1,prob=y.p)

# Much slower... (with no prediction!).
gp4 <- stan_model(here('stan_scripts','gpinfer1D_binom.stan'))
fit <- sampling(gp4,list(N=length(x),D=1,x=matrix(x,ncol=1),y=y),iter=1000,chains=3,cores=3)





