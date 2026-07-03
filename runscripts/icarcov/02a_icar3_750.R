# Load Library ----
library(here)
library(nimbleCarbon)
library(parallel)

# Load Data ----
load(here('data','01_dataprep_out_d750.RData'))

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
out.a3.750  <- parLapply(cl = cl, X = seeds, fun = run.a3.icar, code = icarmodel.a3, d = dat.icar.d750, constants = constants.icar.d750, theta.init = theta.init.icar.d750, niter = niter, nburnin = nburnin, thin = thin)  
stopCluster(cl)

# Model Diagnostics ----
#rhat
rhats.a3.750 <- coda::gelman.diag(coda::mcmc.list(out.a3.750))
if(any(rhats.a3.750[[1]][,1]>1.01)){rhats.a3.750[[1]][which(rhats.a3.750[[1]][,1]>1.01),1]}
#ess
ess.a3.750 <- coda::effectiveSize(coda::mcmc.list(out.a3.750))
#waic
m.a3.750 <- nimbleModel(code=icarmodel.a3, data = dat.icar.d750, constants = constants.icar.d750)
Cm.a3.750 <- compileNimble(m.a3.750)
samples.a3.750 <- do.call(rbind.data.frame,out.a3.750)
waic.a3.750 <- calculateWAIC(samples.a3.750,m.a3.750)
posterior.a3.750 <- samples.a3.750[,!grepl('theta',colnames(samples.a3.750))]

# Save ----
save(posterior.a3.750,rhats.a3.750,ess.a3.750,waic.a3.750,file=here('results','02a_fitmodel3_750_out.RData'))
