library(here)
library(rstan)
library(mgcv)
library(brms)

# 1-D Gaussian Response ----
mcycle <- MASS::mcycle
plot(accel~times,data=mcycle)
x <- mcycle$times
y <- mcycle$accel
D  <- 1

# Just fitting 1D on gaussian response
# very fast
gp <- stan_model(here('stan_scripts','gpinfer1D.stan'))
fit <- sampling(gp,list(N=length(x),D=1,x=matrix(x,ncol=1),y=y),iter=1000,chains=3,cores=3)


# Fitting + predictions
# Slow
gp2 <- stan_model(here('stan_scripts','gpinfer1Db.stan'))
dat <- list(N1 =  length(x),
	    x1 = x,
	    y1 = y,
	    N2 = length(seq(2,60,1)),
	    x2 = seq(2,60,1))

fit2 <- sampling(gp2,dat,iter=300,chains=3,cores=3)
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
# Generate sample data
library(mgcv)
gam_mod <- gam(accel ~ s(times), data = mcycle)
gam.pred <- plogis(scale(predict(gam_mod)))
x <- mcycle$times
set.seed(123)
y  <- rbinom(length(x),size=1,prob=gam.pred)
cdplot(as.factor(y)~x)
d <- data.frame(y=y,x=x)
gamfit <- gam(y~s(x),data=d,family='binomial',method='REML')
plot(gamfit)
# using brms, ok speed
gamfit2 <- brm(y~gp(x),data=d,family=bernoulli(link='logit'),chains=4,cores=4,niter=3000)
gamfit2 <- brm(y~gp(x),data=d,family=bernoulli(link='logit'),chains=4,cores=4,control=list(adapt_delta=0.95,max_treedepth=15),niter=3000)


# Much slower compared to Gaussian...
gp4 <- stan_model(here('stan_scripts','gpinfer1D_binom.stan'))
fit4 <- sampling(gp4,list(N=length(x),D=1,x=matrix(x,ncol=1),y=y),iter=2000,chains=4,cores=4)
posterior  <-  extract(fit4)




# Making predictions (runing but possibly biased)
gp5  <- stan_model(here('stan_scripts','gpinfer1Db_binom.stan'))
dat  <- list(N1=length(x),N2=length(seq(2,60,1)),D=1,x1=matrix(x,ncol=1),x2=matrix(seq(2,60,1),ncol=1),y1=y)
fit5 <- sampling(gp5,dat,iter=2000,chains=3,cores=3)
posterior  <-  extract(fit5)
plot(seq(2,60,1),apply(posterior$gamma0,2,mean)+apply(posterior$f[,c((dat$N1+1):(dat$N1+dat$N2))],2,mean),type='l')
points(x,scale(y))




s
