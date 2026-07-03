# Load Library ----
library(here)
library(nimbleCarbon)
library(parallel)

# Load Data ----
load(here('data','01_dataprep_out_d1000.RData'))

# Load Model ----
source(here('src','model_b2.R'))


# MCMC Settings ----
ncores <- 4
seeds <- c(1,2,3,4)
niter <- 1000000
nburnin <- 500000
thin <- 100

# Fit Models ---- 
cl <- makeCluster(ncores)
out.b2.1000  <- parLapply(cl = cl, X = seeds, fun = run.b2.icar, code = icarmodel.b2, d = dat.icar.d1000, constants = constants.icar.d1000, theta.init = theta.init.icar.d1000, niter = niter, nburnin = nburnin, thin = thin)  
stopCluster(cl)

# Model Diagnostics ----
#rhat
rhats.b2.1000 <- coda::gelman.diag(coda::mcmc.list(out.b2.1000),multivariate = FALSE)
if(any(rhats.b2.1000[[1]][,1]>1.01)){rhats.b2.1000[[1]][which(rhats.b2.1000[[1]][,1]>1.01),1]}
#ess
ess.b2.1000 <- coda::effectiveSize(coda::mcmc.list(out.b2.1000))
#waic
source(here('src','dDoubleFlat.R'))
m.b2.1000 <- nimbleModel(code=icarmodel.b2, data = dat.icar.d1000, constants = constants.icar.d1000)
Cm.b2.1000 <- compileNimble(m.b2.1000)
samples.b2.1000 <- do.call(rbind.data.frame,out.b2.1000)
waic.b2.1000 <- calculateWAIC(samples.b2.1000,m.b2.1000)
posterior.b2.1000 <- samples.b2.1000[,!grepl('theta',colnames(samples.b2.1000))]

# Save ----
save(posterior.b2.1000,rhats.b2.1000,ess.b2.1000,waic.b2.1000,file=here('results','02b_fitmodel2_1000_out.RData'))
