# Load Library
library(here)
library(nimbleCarbon)
library(parallel)
load(here('01_dataprep_out.RData'))

# Model
run <- function(seed,d,constants,theta.init,nburnin,niter,thin)
{
	library(nimbleCarbon)
	#icar model
	icarmodel  <- nimbleCode({
		for (i in 1:N)
		{
# 			theta[i] ~ dDoubleExponentialGrowth(a=a,b=b,r1=s1[id.hex[i]],r2=s2[id.hex[i]],mu=delta0)
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
 		beta1 ~ dnorm(0,0.0001) # effect of ash
	})

	#Setup Init
	inits  <- list()
	inits$theta <- theta.init
	inits$sigma1 <- runif(1,0,100)
	inits$sigma2 <- runif(1,0,100)
	inits$s1 <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$s2 <- rnorm(constants$n.hex,mean=0,sd=0.001)
	inits$beta1 <- rnorm(1,0,0.0001)

	#MCMC
	model <- nimbleModel(icarmodel,constants=constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model,monitors=c('s1','s2','sigma1','sigma2','beta1'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC,project=cModel)
	samples <- runMCMC(cMCMC,niter=niter,thin=thin,nburnin=nburnin,samplesAsCodaMCMC=T,setSeed=seed)
	return(samples)
}

# Run MCMC ----
ncores <- 4
cl <- makeCluster(ncores)
seeds <- c(1,2,3,4)
niter <- 800000
nburnin <- 400000
thin <- 80

out  <- parLapply(cl = cl, X = seeds, fun = run, d = dat.icar, constants = constants.icar, theta.init = theta.init.icar, niter = niter, nburnin = nburnin, thin = thin)  
stopCluster(cl)

# Post-process and save ----
samples  <- coda::mcmc.list(out)
rhats.a <- coda::gelman.diag(samples)
if(any(rhats.a[[1]][,1]>1.01)){rhats.a[[1]][which(rhats.a[[1]][,1]>1.01),1]}
ess.a  <- coda::effectiveSize(samples)
posterior.a <- do.call(rbind.data.frame,samples)
save(posterior.a,rhats.a,ess.a,file=here('results','02a_fitmodel_out.RData'))
