library(here)
library(rstan)
library(dplyr)
library(sf)
library(nimble)
library(rnaturalearth)
library(cascsim)
library(ggplot2)

# Load Data ----
load(here('testscripts','simdata.RData'))

covfun <- function(alpha,rho,d) 
{
	return(alpha^2 * exp(-d^2/(2*rho^2)))
}

plot(1:500,covfun(alpha=1,rho=100,d=1:500),type='l')






simdat <- list(N = nrow(dates.d),
	    coords = as.matrix(dates.d[,c('X','Y')]/1000),
	    x = dates.d$ash,
	    rho=400,
	    alpha=0.5,
	    sigma=0.00001,
	    gamma0=-0.3,
	    gamma1=1.6)
model2c <- stan_model(here('stan_scripts','model2c_gen.stan'))
sim  <- sampling(model2c,simdat,iter=1,warmup=0,chains=1,algorithm='Fixed_param') |> extract()
hist(as.numeric(sim$eta))
hist(plogis(as.numeric(sim$lambda)))

dates.d$p=as.numeric(plogis(sim$lambda))

ggplot(dates.d,aes(x=X,y=Y,col=p)) +
	geom_point() 

dates.d[5,]
plogis(as.numeric(sim$lambda)[5])
plogis(-0.3 + 1.6 * 1 + as.numeric(sim$eta)[5])
plot(p~as.factor(ash),dates.d)



# Stan Analyses
dat <- list(N = nrow(dates.d),
	    coords = as.matrix(dates.d[,c('X','Y')])/1000,
	    y = dates.d$y,
	    x = dates.d$ash)




model2c <- stan_model(here('stan_scripts','model2c.stan'))

fit  <- sampling(model2c,dat,iter=500,chains=4,cores=4)
params <- extract(fit)

plot(plogis(apply(params$lambda,2,mean)),dates.d$p)



# Nimble Analyses ----
# distmat <- as.matrix(dist(rbind(dates.d[,c('X','Y')],st_coordinates(pred.locations.centroid))))
distmat <- as.matrix(dist(dates.d[,c('X','Y')]))
distmat  <- distmat/1000
dat <- list(y=dates.d$y)
# constants <- list(z=c(dates.d$ash,predictions.d$ash),
# 		  N = nrow(dates.d)+nrow(predictions.d),
# 		  N1 = nrow(dates.d),
# 		  dist_mat = distmat)
constants <- list(z=c(dates.d$ash),
		  N = nrow(dates.d),
		  dist_mat = distmat)
cov_ExpQ <- nimbleFunction(run = function(dists = double(2), rho = double(0), etasq = double(0),sigmasq = double(0)) 
			   {
				   returnType(double(2))
				   n <- dim(dists)[1]
				   result <- matrix(nrow = n, ncol = n, init = FALSE)
				   deltaij <- matrix(nrow = n, ncol = n,init = TRUE)
				   diag(deltaij) <- 1
				   for(i in 1:n)
					   for(j in 1:n)
						   result[i, j] <- etasq*exp(-0.5*(dists[i,j]/rho)^2) + sigmasq*deltaij[i,j]
					   return(result)
				   })
Ccov_ExpQ <- compileNimble(cov_ExpQ)

model1 <- nimbleCode({
	for (i in 1:N)
	{
		# Assign probabilities
		p[i]  <- 1/(1+(exp(-(gamma0 + gamma1 * z[i] + s[i]))))
	        y[j] ~ dbern(p[j])
		
	}
	# Spatial model
	mu_s[1:N] <- 0
	cov_s[1:N, 1:N] <- cov_ExpQ(dist_mat[1:N, 1:N], rho, etasq, 0.000001)
	s[1:N] ~ dmnorm(mu_s[1:N], cov = cov_s[1:N, 1:N])

	#Priors
	gamma0 ~ dnorm(0,5)
	gamma1 ~ dnorm(0,5)
	etasq ~ dexp(20);
	rho ~ T(dgamma(10,(10-1)/150),1,1350); #mode 150
})

model2 <- nimbleCode({
	for (i in 1:N)
	{
		# Assign probabilities
		p[i]  <- 1/(1+(exp(-(gamma0 + gamma1 * z[i] + s[i]))))
		
	}
	for (j in 1:N1)
	{
		# Fi to data
	        y[j] ~ dbern(p[j])
	}

	# Spatial model
	mu_s[1:N] <- 0
	cov_s[1:N, 1:N] <- cov_ExpQ(dist_mat[1:N, 1:N], rho, etasq, 0.000001)
	s[1:N] ~ dmnorm(mu_s[1:N], cov = cov_s[1:N, 1:N])

	#Priors
	gamma0 ~ dnorm(0,5)
	gamma1 ~ dnorm(0,5)
	etasq ~ dexp(20);
	rho ~ T(dgamma(10,(10-1)/150),1,1350); #mode 150
})

inits  <-  list()
inits$gamma0 <- 0
inits$gamma1 <- 0
inits$rho  <- rtgamma(1,shape=10,scale=(10-1)/200,min=1,max=1350)
inits$etasq  <- rexp(1,20)
inits$s  <- rep(0,constants$N)
inits$cov_s <- Ccov_ExpQ(constants$dist_mat, inits$rho, inits$etasq, 0.000001)
inits$s <-  t(chol(inits$cov_s)) %*% rnorm(constants$N)
inits$s <- inits$s[ , 1]  # so can give nimble a vector rather than one-column matrix


model.gp <- nimbleModel(model2,constants = constants,data=dat,inits=inits)
cModel.gp <- compileNimble(model.gp)
conf.gp <- configureMCMC(model.gp)
conf.gp$addMonitors('s')
conf.gp$addMonitors('rho')
conf.gp$addMonitors('etasq')
conf.gp$removeSamplers('s[1:575]')
conf.gp$removeSamplers('gamma0')
# conf.gp$removeSamplers('gamma1')
conf.gp$addSampler(c('gamma0','s[1:575]'), type='AF_slice') 
MCMC.gp <- buildMCMC(conf.gp)
cMCMC.gp <- compileNimble(MCMC.gp)
results <- runMCMC(cMCMC.gp, nchain=1,niter = 200, thin=1, nburnin = 100,samplesAsCodaMCMC = T) 




model.gp <- nimbleModel(model2,constants = constants,data=dat,inits=inits)
cModel.gp <- compileNimble(model.gp)
conf.gp <- configureMCMC(model.gp)
conf.gp$addMonitors('s')
conf.gp$addMonitors('rho')
conf.gp$addMonitors('etasq')
conf.gp$removeSamplers('s[1:575]')
conf.gp$removeSamplers('gamma0')
# conf.gp$removeSamplers('gamma1')
conf.gp$addSampler(c('gamma0','s[1:575]'), type='AF_slice') 
MCMC.gp <- buildMCMC(conf.gp)
cMCMC.gp <- compileNimble(MCMC.gp)
results <- runMCMC(cMCMC.gp, nchain=1,niter = 200, thin=1, nburnin = 100,samplesAsCodaMCMC = T) 






# Prepare Data for Stan ----
distmat <- as.matrix(dist(rbind(dates.d[,c('X','Y')],st_coordinates(pred.locations.centroid))))


dat <- list(n = nrow(dates.d) + nrow(predictions.d),
	    n1 = nrow(dates.d),#Number of Samples
	    n2 = nrow(predictions.d), #Number of prediction locations
	    x  = c(dates.d$ash,predictions.d$ash), #Predictors
	    x2  = predictions.d$ash,  #Predictors (only predicted locations)
	    y  = dates.d$y, #Binary responses
	    dist  = distmat) # Distance Matrix (combined)

model2 <- stan_model(here('stan_scripts','model2.stan'))
fit2 <- sampling(model2,dat,iter=100,chains=1)




