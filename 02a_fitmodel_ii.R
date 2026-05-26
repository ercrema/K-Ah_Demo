# Load Library ----
library(here)
library(nimbleCarbon)
library(parallel)

# Load Data ----
load(here('data','01_dataprep_out_d500.RData'))
load(here('data','01_dataprep_out_d750.RData'))
load(here('data','01_dataprep_out_d1000.RData'))

# Load Model ----
source(here('src','model_a.R'))


# MCMC Settings ----
ncores <- 4
seeds <- c(1,2,3,4)
#niter <- 800000
niter <- 2000
#nburnin <- 400000
nburnin <- 1000
#thin <- 80
thin <- 1

# Fit Models ---- 
#500
cl <- makeCluster(ncores)
out.ii.d500  <- parLapply(cl = cl, X = seeds, fun = run.a.ii, code = icarmodel.a.ii, d = dat.icar.d500, constants = constants.icar.d500, theta.init = theta.init.icar.d500, niter = niter, nburnin = nburnin, thin = thin)  
stopCluster(cl)
#750
cl <- makeCluster(ncores)
out.ii.d750  <- parLapply(cl = cl, X = seeds, fun = run.a.ii, code = icarmodel.a.ii, d = dat.icar.d750, constants = constants.icar.d750, theta.init = theta.init.icar.d750, niter = niter, nburnin = nburnin, thin = thin)  
stopCluster(cl)
#1000
cl <- makeCluster(ncores)
out.ii.d1000  <- parLapply(cl = cl, X = seeds, fun = run.a.ii, code = icarmodel.a.ii, d = dat.icar.d1000, constants = constants.icar.d1000, theta.init = theta.init.icar.d1000, niter = niter, nburnin = nburnin, thin = thin)  
stopCluster(cl)

# Model Diagnostics ----
#500
#rhat
rhats.a.ii.d500 <- coda::gelman.diag(coda::mcmc.list(out.ii.d500))
if(any(rhats.a.ii.d500[[1]][,1]>1.01)){rhats.a.ii.d500[[1]][which(rhats.a.ii.d500[[1]][,1]>1.01),1]}
#ess
ess.a.ii.d500 <- coda::effectiveSize(coda::mcmc.list(out.ii.d500))
#waic
m.a.ii.d500 <- nimbleModel(code=icarmodel.a.ii, data = dat.icar.d500, constants = constants.icar.d500)
Cm.a.ii.d500 <- compileNimble(m.a.ii.d500)
samples.a.ii.d500 <- do.call(rbind.data.frame,out.ii.d500)
waic.a.ii.d500 <- calculateWAIC(samples.a.ii.d500,m.a.ii.d500)
posterior.a.ii.d500 <- samples.a.ii.d500[,!grepl('theta',colnames(samples.a.ii.d500))]

#750
#rhat
rhats.a.ii.d750 <- coda::gelman.diag(coda::mcmc.list(out.ii.d750))
if(any(rhats.a.ii.d750[[1]][,1]>1.01)){rhats.a.ii.d750[[1]][which(rhats.a.ii.d750[[1]][,1]>1.01),1]}
#ess
ess.a.ii.d750 <- coda::effectiveSize(coda::mcmc.list(out.ii.d750))
#waic
m.a.ii.d750 <- nimbleModel(code=icarmodel.a.ii, data = dat.icar.d750, constants = constants.icar.d750)
Cm.a.ii.d750 <- compileNimble(m.a.ii.d750)
samples.a.ii.d750 <- do.call(rbind.data.frame,out.ii.d750)
waic.a.ii.d750 <- calculateWAIC(samples.a.ii.d750,m.a.ii.d750)
posterior.a.ii.d750 <- samples.a.ii.d750[,!grepl('theta',colnames(samples.a.ii.d750))]

#1000
#rhat
rhats.a.ii.d1000 <- coda::gelman.diag(coda::mcmc.list(out.ii.d1000))
if(any(rhats.a.ii.d1000[[1]][,1]>1.01)){rhats.a.ii.d1000[[1]][which(rhats.a.ii.d1000[[1]][,1]>1.01),1]}
#ess
ess.a.ii.d1000 <- coda::effectiveSize(coda::mcmc.list(out.ii.d1000))
#waic
m.a.ii.d1000 <- nimbleModel(code=icarmodel.a.ii, data = dat.icar.d1000, constants = constants.icar.d1000)
Cm.a.ii.d1000 <- compileNimble(m.a.ii.d1000)
samples.a.ii.d1000 <- do.call(rbind.data.frame,out.ii.d1000)
waic.a.ii.d1000 <- calculateWAIC(samples.a.ii.d1000,m.a.ii.d1000)
posterior.a.ii.d1000 <- samples.a.ii.d1000[,!grepl('theta',colnames(samples.a.ii.d1000))]

# Save ----
save(posterior.a.ii.d500,rhats.a.ii.d500,ess.a.ii.d500,waic.a.ii.d500,file=here('results','02a_fitmodel_ii_d500_out.RData'))
save(posterior.a.ii.d750,rhats.a.ii.d750,ess.a.ii.d750,waic.a.ii.d750,file=here('results','02a_fitmodel_ii_d750_out.RData'))
save(posterior.a.ii.d1000,rhats.a.ii.d1000,ess.a.ii.d1000,waic.a.ii.d1000,file=here('results','02a_fitmodel_ii_d1000_out.RData'))
