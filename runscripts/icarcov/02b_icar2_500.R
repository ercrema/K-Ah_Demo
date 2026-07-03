# Load Library ----
library(here)
library(nimbleCarbon)
library(parallel)

# Load Data ----
load(here('data','01_dataprep_out_d500.RData'))

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
out.b2.500  <- parLapply(cl = cl, X = seeds, fun = run.b2.icar, code = icarmodel.b2, d = dat.icar.d500, constants = constants.icar.d500, theta.init = theta.init.icar.d500, niter = niter, nburnin = nburnin, thin = thin)  
stopCluster(cl)

# Model Diagnostics ----
#rhat
rhats.b2.500 <- coda::gelman.diag(coda::mcmc.list(out.b2.500))
if(any(rhats.b2.500[[1]][,1]>1.01)){rhats.b2.500[[1]][which(rhats.b2.500[[1]][,1]>1.01),1]}
#ess
ess.b2.500 <- coda::effectiveSize(coda::mcmc.list(out.b2.500))
#waic
source(here('src','dDoubleFlat.R'))
m.b2.500 <- nimbleModel(code=icarmodel.b2, data = dat.icar.d500, constants = constants.icar.d500)
Cm.b2.500 <- compileNimble(m.b2.500)
samples.b2.500 <- do.call(rbind.data.frame,out.b2.500)
waic.b2.500 <- calculateWAIC(samples.b2.500,m.b2.500)
posterior.b2.500 <- samples.b2.500[,!grepl('theta',colnames(samples.b2.500))]

# Save ----
save(posterior.b2.500,rhats.b2.500,ess.b2.500,waic.b2.500,file=here('results','02b_fitmodel2_500_out.RData'))
