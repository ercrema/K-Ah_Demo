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
out.iii.d500  <- parLapply(cl = cl, X = seeds, fun = run.a.iii, code = icarmodel.a.iii, d = dat.icar.d500, constants = constants.icar.d500, theta.init = theta.init.icar.d500, niter = niter, nburnin = nburnin, thin = thin)  
stopCluster(cl)
#750
cl <- makeCluster(ncores)
out.iii.d750  <- parLapply(cl = cl, X = seeds, fun = run.a.iii, code = icarmodel.a.iii, d = dat.icar.d750, constants = constants.icar.d750, theta.init = theta.init.icar.d750, niter = niter, nburnin = nburnin, thin = thin)  
stopCluster(cl)
#1000
cl <- makeCluster(ncores)
out.iii.d1000  <- parLapply(cl = cl, X = seeds, fun = run.a.iii, code = icarmodel.a.iii, d = dat.icar.d1000, constants = constants.icar.d1000, theta.init = theta.init.icar.d1000, niter = niter, nburnin = nburnin, thin = thin)  
stopCluster(cl)

# Model Diagnostics ----
#500
#rhat
rhats.a.iii.d500 <- coda::gelman.diag(coda::mcmc.list(out.iii.d500))
if(any(rhats.a.iii.d500[[1]][,1]>1.01)){rhats.a.iii.d500[[1]][which(rhats.a.iii.d500[[1]][,1]>1.01),1]}
#ess
ess.a.iii.d500 <- coda::effectiveSize(coda::mcmc.list(out.iii.d500))
#waic
m.a.iii.d500 <- nimbleModel(code=icarmodel.a.iii, data = dat.icar.d500, constants = constants.icar.d500)
Cm.a.iii.d500 <- compileNimble(m.a.iii.d500)
samples.a.iii.d500 <- do.call(rbind.data.frame,out.iii.d500)
waic.a.iii.d500 <- calculateWAIC(samples.a.iii.d500,m.a.iii.d500)
posterior.a.iii.d500 <- samples.a.iii.d500[,!grepl('theta',colnames(samples.a.iii.d500))]

#750
#rhat
rhats.a.iii.d750 <- coda::gelman.diag(coda::mcmc.list(out.iii.d750))
if(any(rhats.a.iii.d750[[1]][,1]>1.01)){rhats.a.iii.d750[[1]][which(rhats.a.iii.d750[[1]][,1]>1.01),1]}
#ess
ess.a.iii.d750 <- coda::effectiveSize(coda::mcmc.list(out.i.d750))
#waic
m.a.iii.d750 <- nimbleModel(code=icarmodel.a.iii, data = dat.icar.d750, constants = constants.icar.d750)
Cm.a.iii.d750 <- compileNimble(m.a.iii.d750)
samples.a.iii.d750 <- do.call(rbind.data.frame,out.iii.d750)
waic.a.iii.d750 <- calculateWAIC(samples.a.iii.d750,m.a.iii.d750)
posterior.a.iii.d750 <- samples.a.iii.d750[,!grepl('theta',colnames(samples.a.iii.d750))]

#1000
#rhat
rhats.a.iii.d1000 <- coda::gelman.diag(coda::mcmc.list(out.iii.d1000))
if(any(rhats.a.iii.d1000[[1]][,1]>1.01)){rhats.a.iii.d1000[[1]][which(rhats.a.iii.d1000[[1]][,1]>1.01),1]}
#ess
ess.a.iii.d1000 <- coda::effectiveSize(coda::mcmc.list(out.iii.d1000))
#waic
m.a.iii.d1000 <- nimbleModel(code=icarmodel.a.iii, data = dat.icar.d1000, constants = constants.icar.d1000)
Cm.a.iii.d1000 <- compileNimble(m.a.iii.d1000)
samples.a.iii.d1000 <- do.call(rbind.data.frame,out.iii.d1000)
waic.a.iii.d1000 <- calculateWAIC(samples.a.iii.d1000,m.a.iii.d1000)
posterior.a.iii.d1000 <- samples.a.iii.d1000[,!grepl('theta',colnames(samples.a.iii.d1000))]

# Save ----
save(posterior.a.iii.d500,rhats.a.iii.d500,ess.a.iii.d500,waic.a.iii.d500,file=here('results','02a_fitmodel_iii_d500_out.RData'))
save(posterior.a.iii.d750,rhats.a.iii.d750,ess.a.iii.d750,waic.a.iii.d750,file=here('results','02a_fitmodel_iii_d750_out.RData'))
save(posterior.a.iii.d1000,rhats.a.iii.d1000,ess.a.iii.d1000,waic.a.iii.d1000,file=here('results','02a_fitmodel_iii_d1000_out.RData'))
