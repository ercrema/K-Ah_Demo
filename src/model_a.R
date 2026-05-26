# Model A-i: ICAR model without covariates
icarmodel.a.i  <- nimbleCode({
	for (i in 1:N)
	{
		theta[i] ~ dDoubleExponentialGrowth(a=a,b=b,r1=s1[id.hex[i]],r2=s2[id.hex[i]],mu=delta0)
		c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
		sigmaDate[i] <- (cra.error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
	}

	# ICAR Model Prior
	s1[1:n.hex] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n.hex], tau1, zero_mean =0)
	s2[1:n.hex] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n.hex], tau2, zero_mean =0)
	tau1 <- 1/sigma1^2
	tau2 <- 1/sigma2^2
	sigma1 ~ dunif(0,100)
	sigma2 ~ dunif(0,100)
})

run.a.i <- function(seed,d,code,constants,theta.init,nburnin,niter,thin)
{
	library(nimbleCarbon)
	#Setup Init
	inits  <- list()
	inits$theta <- theta.init
	inits$sigma1 <- runif(1,0,100)
	inits$sigma2 <- runif(1,0,100)
	inits$s1 <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$s2 <- rnorm(constants$n.hex,mean=0,sd=0.001)

	#MCMC
	model <- nimbleModel(code,constants=constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,monitors=c('s1','s2','sigma1','sigma2','theta'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC,project=cModel)
	samples <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,samplesAsCodaMCMC=T,setSeed=seed)
	return(samples)
}


# Model A-ii: ICAR model with Ash as covaqriate for r2
icarmodel.a.ii  <- nimbleCode({
	for (i in 1:N)
	{
		theta[i] ~ dDoubleExponentialGrowth(a=a,b=b,r1=s1[id.hex[i]],r2=s2[id.hex[i]] + beta1 * ash[id.hex[i]],mu=delta0)
		c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
		sigmaDate[i] <- (cra.error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
	}

	# ICAR Model Prior
	s1[1:n.hex] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n.hex], tau1, zero_mean =0)
	s2[1:n.hex] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n.hex], tau2, zero_mean =0)
	tau1 <- 1/sigma1^2
	tau2 <- 1/sigma2^2
	sigma1 ~ dunif(0,100)
	sigma2 ~ dunif(0,100)

	# Ash Prior
	beta1 ~ dnorm(0,0.001) # effect of ash
})

run.a.ii <- function(seed,d,code,constants,theta.init,nburnin,niter,thin)
{
	library(nimbleCarbon)
	#Setup Init
	inits  <- list()
	inits$theta <- theta.init
	inits$sigma1 <- runif(1,0,100)
	inits$sigma2 <- runif(1,0,100)
	inits$s1 <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$s2 <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$beta1 <- rnorm(1,0,0.001)

	#MCMC
	model <- nimbleModel(code,constants=constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,monitors=c('s1','s2','sigma1','sigma2','theta','beta1'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC,project=cModel)
	samples <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,samplesAsCodaMCMC=T,setSeed=seed)
	return(samples)
}

# Model A-iii: ICAR model with dist as covariate for r2
icarmodel.a.iii  <- nimbleCode({
	for (i in 1:N)
	{
		theta[i] ~ dDoubleExponentialGrowth(a=a,b=b,r1=s1[id.hex[i]],r2=s2[id.hex[i]] + beta1 * dist[id.hex[i]],mu=delta0)
		c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
		sigmaDate[i] <- (cra.error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
	}

	# ICAR Model Prior
	s1[1:n.hex] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n.hex], tau1, zero_mean =0)
	s2[1:n.hex] ~ dcar_normal(adj[1:L], weights[1:L], num[1:n.hex], tau2, zero_mean =0)
	tau1 <- 1/sigma1^2
	tau2 <- 1/sigma2^2
	sigma1 ~ dunif(0,100)
	sigma2 ~ dunif(0,100)

	# Ash Prior
	beta1 ~ dnorm(0,0.001) # effect of dist
})

run.a.iii <- function(seed,d,code,constants,theta.init,nburnin,niter,thin)
{
	library(nimbleCarbon)
	#Setup Init
	inits  <- list()
	inits$theta <- theta.init
	inits$sigma1 <- runif(1,0,100)
	inits$sigma2 <- runif(1,0,100)
	inits$s1 <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$s2 <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$beta1 <- rnorm(1,0,0.001)

	#MCMC
	model <- nimbleModel(code,constants=constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,monitors=c('s1','s2','sigma1','sigma2','theta','beta1'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC,project=cModel)
	samples <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,samplesAsCodaMCMC=T,setSeed=seed)
	return(samples)
}



