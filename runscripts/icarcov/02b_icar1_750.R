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
out.b1.750  <- parLapply(cl = cl, X = seeds, fun = run.b1.icar, code = icarmodel.b1, d = dat.icar.d750, constants = constants.icar.d750, theta.init = theta.init.icar.d750, niter = niter, nburnin = nburnin, thin = thin)  
stopCluster(cl)

# Model Diagnostics ----
#rhat
rhats.b1.750 <- coda::gelman.diag(coda::mcmc.list(out.b1.750))
if(any(rhats.b1.750[[1]][,1]>1.01)){rhats.b1.750[[1]][which(rhats.b1.750[[1]][,1]>1.01),1]}
#ess
ess.b1.750 <- coda::effectiveSize(coda::mcmc.list(out.b1.750))
#waic
source(here('src','dDoubleFlat.R'))
m.b1.750 <- nimbleModel(code=icarmodel.b1, data = dat.icar.d750, constants = constants.icar.d750)
Cm.b1.750 <- compileNimble(m.b1.750)
samples.b1.750 <- do.call(rbind.data.frame,out.b1.750)
waic.b1.750 <- calculateWAIC(samples.b1.750,m.b1.750)
posterior.b1.750 <- samples.b1.750[,!grepl('theta',colnames(samples.b1.750))]

# Save ----
save(posterior.b1.750,rhats.b1.750,ess.b1.750,waic.b1.750,file=here('results','02b_fitmodel1_750_out.RData'))
