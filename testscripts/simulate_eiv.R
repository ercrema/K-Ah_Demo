library(here)
library(dplyr)
library(rstan)
library(here)

n.sites <- 100
n.dates <- 250
tephra <- rbinom(n.sites,size=1,prob=0.3)
site.id <- 1:n.sites
d <- data.frame(site.id=site.id,tephra=tephra)
gamma0 <- 0.8
gamma1 <- -2
d$p_before <- plogis(gamma0 + tephra * gamma1)
d.dates <- data.frame(site.id=c(1:n.sites,sample(n.sites,size=n.dates-n.sites,replace=T)))
d.dates <- left_join(d.dates,d)
d.dates$y <- rbinom(n=n.dates,size=1,prob=d.dates$p_before)

model0 <- stan_model(here('stan_scripts','model0.stan'))
fit0 <- sampling(model0,list(n=n.dates,x=d.dates$tephra,y=d.dates$y),iter=1000,chains=3)
post0 <- extract(fit0)
par(mfrow=c(2,1))
hist(post0$gamma0,border=NA)
abline(v=gamma0)
hist(post0$gamma1,border=NA)
abline(v=gamma1)


d.dates$p_y = d.dates$y
i <- sample(nrow(d.dates),size=20)
d.dates$p_y[i] <- d.dates$p_y[i] + runif(20,min=-0.3,max=0.3)
d.dates$p_y[which(d.dates$p_y<0)]=0
d.dates$p_y[which(d.dates$p_y>1)]=1


model1 <- stan_model(here('stan_scripts','model1.stan'))
fit1 <- sampling(model1,list(n=n.dates,x=d.dates$tephra,p_y=d.dates$p_y),iter=1000,chains=3)
post1 <- extract(fit1)
par(mfrow=c(2,1))
hist(post1$gamma0,border=NA)
abline(v=gamma0)
hist(post1$gamma1,border=NA)
abline(v=gamma1)







