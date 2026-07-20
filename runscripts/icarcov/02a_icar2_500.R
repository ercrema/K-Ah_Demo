# Load Library ----
library(here)
library(nimbleCarbon)
library(parallel)

# Load Data ----
load(here('data','01_dataprep_out_d500.RData'))

# Load Model ----
source(here('src','model_a.R'))


# MCMC Settings ----
ncores <- 4
seeds <- c(1,2,3,4)
niter <- 1000000
nburnin <- 500000
thin <- 100

# Fit Models ---- 
cl <- makeCluster(ncores)
out.a2.500  <- parLapply(cl = cl, X = seeds, fun = run.a2.icar.sens, code = icarmodel.a2, d = dat.icar.d500, constants = constants.icar.d500, theta.init = theta.init.icar.d500, niter = niter, nburnin = nburnin, thin = thin)  
stopCluster(cl)

# Model Diagnostics ----
#rhat
rhats.a2.500 <- coda::gelman.diag(coda::mcmc.list(out.a2.500))
if(any(rhats.a2.500[[1]][,1]>1.01)){rhats.a2.500[[1]][which(rhats.a2.500[[1]][,1]>1.01),1]}
#ess
ess.a2.500 <- coda::effectiveSize(coda::mcmc.list(out.a2.500))
posterior.a2.500 <- do.call(rbind.data.frame,out.a2.500)

# Save ----
save(posterior.a2.500,rhats.a2.500,ess.a2.500,file=here('results','02a_fitmodel2_500_out.RData'))
