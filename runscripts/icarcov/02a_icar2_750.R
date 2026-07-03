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
out.a2.750  <- parLapply(cl = cl, X = seeds, fun = run.a2.icar, code = icarmodel.a2, d = dat.icar.d750, constants = constants.icar.d750, theta.init = theta.init.icar.d750, niter = niter, nburnin = nburnin, thin = thin)  
stopCluster(cl)

# Model Diagnostics ----
#rhat
rhats.a2.750 <- coda::gelman.diag(coda::mcmc.list(out.a2.750))
if(any(rhats.a2.750[[1]][,1]>1.01)){rhats.a2.750[[1]][which(rhats.a2.750[[1]][,1]>1.01),1]}
#ess
ess.a2.750 <- coda::effectiveSize(coda::mcmc.list(out.a2.750))
#waic
m.a2.750 <- nimbleModel(code=icarmodel.a2, data = dat.icar.d750, constants = constants.icar.d750)
Cm.a2.750 <- compileNimble(m.a2.750)
samples.a2.750 <- do.call(rbind.data.frame,out.a2.750)
waic.a2.750 <- calculateWAIC(samples.a2.750,m.a2.750)
posterior.a2.750 <- samples.a2.750[,!grepl('theta',colnames(samples.a2.750))]

# Save ----
save(posterior.a2.750,rhats.a2.750,ess.a2.750,waic.a2.750,file=here('results','02a_fitmodel2_750_out.RData'))
