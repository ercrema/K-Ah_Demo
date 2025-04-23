library(here)
library(rstan)
library(dplyr)
library(sf)
library(nimble)
library(rnaturalearth)
library(cascsim)
library(ggplot2)
library(brms)

# Load Data ----
load(here('testscripts','simdata.RData'))

# Fit using brms ----
priors <- prior(normal(0, 1), class = "Intercept") + 
	prior(normal(0,1),class='b') + 
	prior(inv_gamma(6, 2000), class = "lscale", coef = "gpXY") +
	prior(normal(0, 0.25), class = "sdgp", coef='gpXY')

fitted.model <- brm(y ~ ash + gp(X,Y), family=bernoulli(link='logit'),data=dates.d,cores=4,chains=4, prior=priors,control=list(adapt_delta=0.99))
post <- as_draws_matrix(fitted.model)



covfun <- function(alpha,rho,d) 
{
	return(alpha^2 * exp(-d^2/(2*rho^2)))
}

plot(1:500,covfun(alpha=1,rho=100,d=1:500),type='l')



# Simulate data using GP model ----
simdat <- list(N = nrow(dates.d),
	    coords = as.matrix(dates.d[,c('X','Y')]/1000),
	    x = dates.d$ash,
	    rho=400,
	    alpha=0.5,
	    sigma=0.00001,
	    gamma0=-0.3,
	    gamma1=1.6)
gpsim <- stan_model(here('stan_scripts','gpsim.stan'))
sim  <- sampling(gpsim,simdat,iter=1,warmup=0,chains=1,algorithm='Fixed_param') |> extract()
hist(as.numeric(sim$eta))
hist(plogis(as.numeric(sim$lambda)))

dates.d$p=as.numeric(plogis(sim$lambda))

ggplot(dates.d,aes(x=X,y=Y,col=p)) +
	geom_point() 

dates.d[5,]
plogis(as.numeric(sim$lambda)[5])
plogis(-0.3 + 1.6 * 1 + as.numeric(sim$eta)[5])
plot(p~as.factor(ash),dates.d)



# Estimate using brms (worked) ----


# Make prediction on new locations
p
























