# Model B-i: ICAR model with no covariates
icarmodel.b0  <- nimbleCode({
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


run.b0.icar <- function(seed,d,code,constants,theta.init,nburnin,niter,thin)
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



icarmodel.b1  <- nimbleCode({
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
	beta1 ~ dnorm(mean=0,sd=0.01) # effect of dist
})


run.b1.icar <- function(seed,d,code,constants,theta.init,nburnin,niter,thin)
{
	library(nimbleCarbon)
	library(here)
	source(here('src','dDoubleFlat.R'))
	#Setup Init
	inits  <- list()
	inits$theta <- theta.init
	inits$sigma1 <- runif(1,0,100)
	inits$eta <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$beta1 <- rnorm(1,0,0.0001)

	#MCMC
	model <- nimbleModel(code,constants=constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,monitors=c('eta','sigma1','beta1','theta'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC,project=cModel)
	samples <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,samplesAsCodaMCMC=T,setSeed=seed)
	return(samples)
}



icarmodel.b2  <- nimbleCode({
	for (i in 1:N)
	{

		theta[i] ~ dDoubleFlat(a=a,b=b,mu=delta0,eta=eta[id.hex[i]] + beta1 * logAsh[id.hex[i]] + beta2 * max(0,logAsh[id.hex[i]]-kappa))
		c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
		sigmaDate[i] <- (cra.error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
	}

	# ICAR Model Prior
	eta[1:n.hex] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n.hex], tau, zero_mean = 0)
	tau <- 1/sigma1^2
	sigma1 ~ dunif(0,100)
	beta1 ~ dnorm(mean=0,sd=1) # effect of dist
	beta2 ~ dnorm(mean=0,sd=1) # effect of ash above kappa
	kappa ~ T(dnorm(mean=2,sd=0.1),0,3.21) # change point for ash
})


run.b2.icar <- function(seed,d,code,constants,theta.init,nburnin,niter,thin)
{
	library(nimbleCarbon)
	library(truncnorm)
	library(here)
	source(here('src','dDoubleFlat.R'))
	#Setup Init
	inits  <- list()
	inits$theta <- theta.init
	inits$sigma1 <- runif(1,0,100)
	inits$eta <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$beta1 <- rnorm(1,0,0.001)
	inits$beta2 <- rnorm(1,0,0.001)
	inits$kappa <- rtruncnorm(1,0,3.21,2,0.1)

	#MCMC
	model <- nimbleModel(code,constants=constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,monitors=c('eta','sigma1','beta1','beta2','kappa','theta'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC,project=cModel)
	samples <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,samplesAsCodaMCMC=T,setSeed=seed)
	return(samples)
}

run.b2.icar.sens <- function(seed,d,code,constants,theta.init,nburnin,niter,thin)
{
	library(nimbleCarbon)
	library(truncnorm)
	library(here)
	source(here('src','dDoubleFlat.R'))
	#Setup Init
	inits  <- list()
	inits$theta <- theta.init
	inits$sigma1 <- runif(1,0,100)
	inits$eta <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$beta1 <- rnorm(1,0,0.001)
	inits$beta2 <- rnorm(1,0,0.001)
	inits$kappa <- rtruncnorm(1,0,3.21,2,0.1)

	#MCMC
	model <- nimbleModel(code,constants=constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,monitors=c('eta','sigma1','beta1','beta2','kappa'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC,project=cModel)
	samples <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,samplesAsCodaMCMC=T,setSeed=seed)
	return(samples)
}

icarmodel.b3  <- nimbleCode({
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
	beta1 ~ dnorm(mean=0,sd=0.01) # effect of dist
})


run.b3.icar <- function(seed,d,code,constants,theta.init,nburnin,niter,thin)
{
	library(nimbleCarbon)
	library(here)
	source(here('src','dDoubleFlat.R'))
	#Setup Init
	inits  <- list()
	inits$theta <- theta.init
	inits$sigma1 <- runif(1,0,100)
	inits$eta <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$beta1 <- rnorm(1,0,0.0001)

	#MCMC
	model <- nimbleModel(code,constants=constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,monitors=c('eta','sigma1','beta1','theta'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC,project=cModel)
	samples <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,samplesAsCodaMCMC=T,setSeed=seed)
	return(samples)
}











