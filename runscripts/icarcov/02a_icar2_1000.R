# Load Library ----
library(here)
library(nimbleCarbon)
library(parallel)

# Load Data ----
load(here('data','01_dataprep_out_d1000.RData'))

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
out.a2.1000  <- parLapply(cl = cl, X = seeds, fun = run.a2.icar, code = icarmodel.a2, d = dat.icar.d1000, constants = constants.icar.d1000, theta.init = theta.init.icar.d1000, niter = niter, nburnin = nburnin, thin = thin)  
stopCluster(cl)

# Model Diagnostics ----
#rhat
rhats.a2.1000 <- coda::gelman.diag(coda::mcmc.list(out.a2.1000))
if(any(rhats.a2.1000[[1]][,1]>1.01)){rhats.a2.1000[[1]][which(rhats.a2.1000[[1]][,1]>1.01),1]}
#ess
ess.a2.1000 <- coda::effectiveSize(coda::mcmc.list(out.a2.1000))
#waic
m.a2.1000 <- nimbleModel(code=model.a2.1000, data = dat.icar.d1000, constants = constants.icar.d1000)
Cm.a2.1000 <- compileNimble(m.a2.1000)
samples.a2.1000 <- do.call(rbind.data.frame,out.a2.1000)
waic.a2.1000 <- calculateWAIC(samples.a2.1000,m.a2.1000)
posterior.a2.1000 <- samples.a2.1000[,!grepl('theta',colnames(samples.a2.1000))]

# Save ----
save(posterior.a2.1000,rhats.a2.1000,ess.a2.1000,waic.a2.1000,file=here('results','02a_fitmodel2_1000_out.RData'))
