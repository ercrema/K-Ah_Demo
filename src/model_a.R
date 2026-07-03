# ICAR model without covariates
icarmodel.a0  <- nimbleCode({
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

run.a0.icar <- function(seed,d,code,constants,theta.init,nburnin,niter,thin)
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


# ICAR model with distance as covariate
icarmodel.a1  <- nimbleCode({
	for (i in 1:N)
	{
		theta[i] ~ dDoubleExponentialGrowth(a=a,b=b,r1=s1[id.hex[i]],r2= s2[id.hex[i]] + beta1 * dist[id.hex[i]], mu=delta0)
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
	beta1 ~ dnorm(mean=0,sd=0.01) # effect of dist
})

run.a1.icar <- function(seed,d,code,constants,theta.init,nburnin,niter,thin)
{
	library(nimbleCarbon)
	#Setup Init
	inits  <- list()
	inits$theta <- theta.init
	inits$sigma1 <- runif(1,0,100)
	inits$sigma2 <- runif(1,0,100)
	inits$s1 <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$s2 <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$beta1 <- rnorm(1,0,0.0001)

	#MCMC
	model <- nimbleModel(code,constants=constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,monitors=c('s1','s2','sigma1','sigma2','beta1','theta'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC,project=cModel)
	samples <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,samplesAsCodaMCMC=T,setSeed=seed)
	return(samples)
}


# ICAR model with change point for ash as covariate
icarmodel.a2  <- nimbleCode({
	for (i in 1:N)
	{
		theta[i] ~ dDoubleExponentialGrowth(a=a,b=b,r1=s1[id.hex[i]],r2= s2[id.hex[i]] + beta1 * logAsh[id.hex[i]] + beta2 * max(0,logAsh[id.hex[i]]-kappa), mu=delta0)
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
	beta1 ~ dnorm(mean=0,sd=0.1) 
	beta2 ~ dnorm(mean=0,sd=0.1) 
	kappa ~ T(dnorm(mean=2,sd=0.1),0,3.21) # change point for ash
})

run.a2.icar <- function(seed,d,code,constants,theta.init,nburnin,niter,thin)
{
	library(nimbleCarbon)
	library(truncnorm)
	#Setup Init
	inits  <- list()
	inits$theta <- theta.init
	inits$sigma1 <- runif(1,0,100)
	inits$sigma2 <- runif(1,0,100)
	inits$s1 <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$s2 <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$beta1 <- rnorm(1,0,0.0001)
	inits$beta2 <- rnorm(1,0,0.0001)
	inits$kappa <- rtruncnorm(1,0,3.21,2,0.1)

	#MCMC
	model <- nimbleModel(code,constants=constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,monitors=c('s1','s2','sigma1','sigma2','beta1','beta2','kappa','theta'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC,project=cModel)
	samples <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,samplesAsCodaMCMC=T,setSeed=seed)
	return(samples)
}


# ICAR model with ash as covariate
icarmodel.a3  <- nimbleCode({
	for (i in 1:N)
	{
		theta[i] ~ dDoubleExponentialGrowth(a=a,b=b,r1=s1[id.hex[i]],r2= s2[id.hex[i]] + beta1 * ash[id.hex[i]], mu=delta0)
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
	beta1 ~ dnorm(mean=0,sd=0.01) # effect of dist
})

run.a3.icar <- function(seed,d,code,constants,theta.init,nburnin,niter,thin)
{
	library(nimbleCarbon)
	#Setup Init
	inits  <- list()
	inits$theta <- theta.init
	inits$sigma1 <- runif(1,0,100)
	inits$sigma2 <- runif(1,0,100)
	inits$s1 <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$s2 <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$beta1 <- rnorm(1,0,0.0001)

	#MCMC
	model <- nimbleModel(code,constants=constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,monitors=c('s1','s2','sigma1','sigma2','beta1','theta'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC,project=cModel)
	samples <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,samplesAsCodaMCMC=T,setSeed=seed)
	return(samples)
}
