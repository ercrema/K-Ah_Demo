# Model B-i: ICAR model with no covariates
icarmodel.b.i  <- nimbleCode({
	for (i in 1:N)
	{

		theta[i] ~ dDoubleFlat(a=a,b=b,mu=delta0,eta=eta[id.hex[i]])
		c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
		sigmaDate[i] <- (cra.error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
	}

	# ICAR Model Prior
	eta[1:n.hex] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n.hex], tau, zero_mean = 0)
	tau <- 1/sigma1^2
	sigma1 ~ dunif(0,100)
})


run.b.i <- function(seed,d,code,constants,theta.init,nburnin,niter,thin)
{
	library(nimbleCarbon)
	library(here)
	source(here('src','dDoubleFlat.R'))
	#Setup Init
	inits  <- list()
	inits$theta <- theta.init
	inits$sigma1 <- runif(1,0,100)
	inits$eta <- rnorm(constants$n.hex,mean=0,sd=0.001)

	#MCMC
	model <- nimbleModel(code,constants=constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,monitors=c('eta','sigma1','theta'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC,project=cModel)
	samples <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,samplesAsCodaMCMC=T,setSeed=seed)
	return(samples)
}

# Model B-ii: ICAR model with Ash as a covariate
icarmodel.b.ii  <- nimbleCode({
	for (i in 1:N)
	{

		theta[i] ~ dDoubleFlat(a=a,b=b,mu=delta0,eta=eta[id.hex[i]] + beta1 * ash[id.hex[i]])
		c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
		sigmaDate[i] <- (cra.error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
	}

	# ICAR Model Prior
	eta[1:n.hex] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n.hex], tau, zero_mean = 0)
	tau <- 1/sigma1^2
	sigma1 ~ dunif(0,100)

	# Ash Prior
	beta1 ~ dnorm(0,0.001) # effect of ash
})


run.b.ii <- function(seed,d,code,constants,theta.init,nburnin,niter,thin)
{
	library(nimbleCarbon)
	library(here)
	source(here('src','dDoubleFlat.R'))
	#Setup Init
	inits  <- list()
	inits$theta <- theta.init
	inits$sigma1 <- runif(1,0,100)
	inits$eta <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$beta1 <- rnorm(1,0,0.001)

	#MCMC
	model <- nimbleModel(code,constants=constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,monitors=c('eta','sigma1','theta','beta1'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC,project=cModel)
	samples <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,samplesAsCodaMCMC=T,setSeed=seed)
	return(samples)
}

# Model B-iii: ICAR model with Dist as a covariate
icarmodel.b.iii  <- nimbleCode({
	for (i in 1:N)
	{

		theta[i] ~ dDoubleFlat(a=a,b=b,mu=delta0,eta=eta[id.hex[i]] + beta1 * dist[id.hex[i]])
		c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
		sigmaDate[i] <- (cra.error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
	}

	# ICAR Model Prior
	eta[1:n.hex] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n.hex], tau, zero_mean = 0)
	tau <- 1/sigma1^2
	sigma1 ~ dunif(0,100)

	# Dist Prior
	beta1 ~ dnorm(0,0.001) # effect of dist
})


run.b.iii <- function(seed,d,code,constants,theta.init,nburnin,niter,thin)
{
	library(nimbleCarbon)
	library(here)
	source(here('src','dDoubleFlat.R'))
	#Setup Init
	inits  <- list()
	inits$theta <- theta.init
	inits$sigma1 <- runif(1,0,100)
	inits$eta <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$beta1 <- rnorm(1,0,0.001)

	#MCMC
	model <- nimbleModel(code,constants=constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,monitors=c('eta','sigma1','theta','beta1'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC,project=cModel)
	samples <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,samplesAsCodaMCMC=T,setSeed=seed)
	return(samples)
}
