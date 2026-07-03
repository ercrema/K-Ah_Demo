# Load Library ----
library(here)
library(nimbleCarbon)
library(parallel)

# Load Data ----
load(here('data','01_dataprep_out_d750.RData'))

# Load Model ----
source(here('src','model_b.R'))


# MCMC Settings ----
ncores <- 4
seeds <- c(1,2,3,4)
niter <- 1000000
nburnin <- 500000
thin <- 100

# Fit Models ---- 
cl <- makeCluster(ncores)
out.b0.750  <- parLapply(cl = cl, X = seeds, fun = run.b0.icar, code = icarmodel.b0, d = dat.icar.d750, constants = constants.icar.d750, theta.init = theta.init.icar.d750, niter = niter, nburnin = nburnin, thin = thin)  
stopCluster(cl)

# Model Diagnostics ----
#rhat
rhats.b0.750 <- coda::gelman.diag(coda::mcmc.list(out.b0.750))
if(any(rhats.b0.750[[1]][,1]>1.01)){rhats.b0.750[[1]][which(rhats.b0.750[[1]][,1]>1.01),1]}
#ess
ess.b0.750 <- coda::effectiveSize(coda::mcmc.list(out.b0.750))
#waic
source(here('src','dDoubleFlat.R'))
m.b0.750 <- nimbleModel(code=icarmodel.b0, data = dat.icar.d750, constants = constants.icar.d750)
Cm.b0.750 <- compileNimble(m.b0.750)
samples.b0.750 <- do.call(rbind.data.frame,out.b0.750)
waic.b0.750 <- calculateWAIC(samples.b0.750,m.b0.750)
posterior.b0.750 <- samples.b0.750[,!grepl('theta',colnames(samples.b0.750))]

# Save ----
save(posterior.b0.750,rhats.b0.750,ess.b0.750,waic.b0.750,file=here('results','02b_fitmodel0_750_out.RData'))
