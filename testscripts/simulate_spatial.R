library(here)
library(rstan)
library(dplyr)
library(sf)
library(nimble)
library(rnaturalearth)
library(cascsim)
# Code from https://aheblog.com/2016/12/07/geostatistical-modelling-with-r-and-stan/

source(here('nimble_scripts','gpSIM.R'))

# settings
n.locations <- 200
n.samples <- 300
perc_misclass <- 0.2


# generate random locations in the main islands of japan
japan <- ne_countries(country = "Japan", returnclass = "sf",scale=10) |> st_transform(6684) |> st_geometry() |> st_cast('POLYGON') |> st_as_sf()
japan$area <- st_area(japan)
japan <- arrange(japan,desc(area)) |> slice_head(n = 4) |> st_geometry()
locations <- st_sample(japan,size=n.locations,type='random')
pred.locations <- st_make_grid(japan,cellsize=50000,square=FALSE)
ii <- apply(st_intersects(japan,pred.locations,sparse=F),2,any)
pred.locations.centroid  <- st_centroid(pred.locations[ii])
pred.locations <- st_intersection(japan,pred.locations)

sites.d <- data.frame(site.id=1:n.locations)
sites.d <- cbind.data.frame(sites.d,st_coordinates(locations))

# Generatec an ellipses representing volcanic eruption
inside_ellipse <- function(x, y, h, k, a, b, theta = 0) {
  # Convert theta to radians
  theta <- theta * pi / 180
  # Rotate points
  x_prime <- (x - h) * cos(theta) + (y - k) * sin(theta)
  y_prime <- -(x - h) * sin(theta) + (y - k) * cos(theta)
  # Check the ellipse condition
  inside <- (x_prime^2 / a^2) + (y_prime^2 / b^2) <= 1
  return(inside)
}
# Simulate ash-fall
j <- inside_ellipse(x=sites.d$X,y=sites.d$Y,h=646610,k=652395,a=1000*1000,b=150*1000,theta=35)
sites.d$ash <- 0
sites.d$ash[j] <- 1

predictions.d <- as.data.frame(st_coordinates(pred.locations.centroid))
j <- inside_ellipse(x=predictions.d$X,y=predictions.d$Y,h=646610,k=652395,a=1000*1000,b=150*1000,theta=35)
predictions.d$ash <- 0
predictions.d$ash[j] <- 1

# Simulate response
res <- gpSim(x=sites.d$X,y=sites.d$Y,z=sites.d$ash,gamma0=-0.2,gamma1=1.2)

# Combine to original
sites.d <- left_join(sites.d,res,by=c('site.id'='ID'))

# Simulate binary dates 
dates.d <- data.frame(site.id=c(1:n.locations,sample(n.locations,size=n.samples-n.locations,replace=T)))
dates.d <- left_join(dates.d,sites.d)
dates.d$y <- rbinom(n=nrow(dates.d),size=1,prob=dates.d$p)
dates.d$p_y<- dates.d$y
i <- sample(nrow(dates.d),size=round(n.samples*perc_misclass))
dates.d$p_y[i] <- dates.d$p_y[i] + runif(round(n.samples*perc_misclass),min=-0.3,max=0.3)
dates.d$p_y[which(dates.d$p_y<0)]=0
dates.d$p_y[which(dates.d$p_y>1)]=1

# Nimble Analyses ----
# distmat <- as.matrix(dist(rbind(dates.d[,c('X','Y')],st_coordinates(pred.locations.centroid))))
distmat <- as.matrix(dist(dates.d[,c('X','Y')]))
distmat  <- distmat/1000
dat <- list(y=dates.d$y)
# constants <- list(z=c(dates.d$ash,predictions.d$ash),
# 		  N = nrow(dates.d)+nrow(predictions.d),
# 		  N1 = nrow(dates.d),
# 		  dist_mat = distmat)
constants <- list(z=c(dates.d$ash),
		  N = nrow(dates.d),
		  dist_mat = distmat)
cov_ExpQ <- nimbleFunction(run = function(dists = double(2), rho = double(0), etasq = double(0),sigmasq = double(0)) 
			   {
				   returnType(double(2))
				   n <- dim(dists)[1]
				   result <- matrix(nrow = n, ncol = n, init = FALSE)
				   deltaij <- matrix(nrow = n, ncol = n,init = TRUE)
				   diag(deltaij) <- 1
				   for(i in 1:n)
					   for(j in 1:n)
						   result[i, j] <- etasq*exp(-0.5*(dists[i,j]/rho)^2) + sigmasq*deltaij[i,j]
					   return(result)
				   })
Ccov_ExpQ <- compileNimble(cov_ExpQ)

model1 <- nimbleCode({
	for (i in 1:N)
	{
		# Assign probabilities
		p[i]  <- 1/(1+(exp(-(gamma0 + gamma1 * z[i] + s[i]))))
	        y[j] ~ dbern(p[j])
		
	}
	# Spatial model
	mu_s[1:N] <- 0
	cov_s[1:N, 1:N] <- cov_ExpQ(dist_mat[1:N, 1:N], rho, etasq, 0.000001)
	s[1:N] ~ dmnorm(mu_s[1:N], cov = cov_s[1:N, 1:N])

	#Priors
	gamma0 ~ dnorm(0,5)
	gamma1 ~ dnorm(0,5)
	etasq ~ dexp(20);
	rho ~ T(dgamma(10,(10-1)/150),1,1350); #mode 150
})

model2 <- nimbleCode({
	for (i in 1:N)
	{
		# Assign probabilities
		p[i]  <- 1/(1+(exp(-(gamma0 + gamma1 * z[i] + s[i]))))
		
	}
	for (j in 1:N1)
	{
		# Fi to data
	        y[j] ~ dbern(p[j])
	}

	# Spatial model
	mu_s[1:N] <- 0
	cov_s[1:N, 1:N] <- cov_ExpQ(dist_mat[1:N, 1:N], rho, etasq, 0.000001)
	s[1:N] ~ dmnorm(mu_s[1:N], cov = cov_s[1:N, 1:N])

	#Priors
	gamma0 ~ dnorm(0,5)
	gamma1 ~ dnorm(0,5)
	etasq ~ dexp(20);
	rho ~ T(dgamma(10,(10-1)/150),1,1350); #mode 150
})

inits  <-  list()
inits$gamma0 <- 0
inits$gamma1 <- 0
inits$rho  <- rtgamma(1,shape=10,scale=(10-1)/200,min=1,max=1350)
inits$etasq  <- rexp(1,20)
inits$s  <- rep(0,constants$N)
inits$cov_s <- Ccov_ExpQ(constants$dist_mat, inits$rho, inits$etasq, 0.000001)
inits$s <-  t(chol(inits$cov_s)) %*% rnorm(constants$N)
inits$s <- inits$s[ , 1]  # so can give nimble a vector rather than one-column matrix


model.gp <- nimbleModel(model2,constants = constants,data=dat,inits=inits)
cModel.gp <- compileNimble(model.gp)
conf.gp <- configureMCMC(model.gp)
conf.gp$addMonitors('s')
conf.gp$addMonitors('rho')
conf.gp$addMonitors('etasq')
conf.gp$removeSamplers('s[1:575]')
conf.gp$removeSamplers('gamma0')
# conf.gp$removeSamplers('gamma1')
conf.gp$addSampler(c('gamma0','s[1:575]'), type='AF_slice') 
MCMC.gp <- buildMCMC(conf.gp)
cMCMC.gp <- compileNimble(MCMC.gp)
results <- runMCMC(cMCMC.gp, nchain=1,niter = 200, thin=1, nburnin = 100,samplesAsCodaMCMC = T) 




model.gp <- nimbleModel(model2,constants = constants,data=dat,inits=inits)
cModel.gp <- compileNimble(model.gp)
conf.gp <- configureMCMC(model.gp)
conf.gp$addMonitors('s')
conf.gp$addMonitors('rho')
conf.gp$addMonitors('etasq')
conf.gp$removeSamplers('s[1:575]')
conf.gp$removeSamplers('gamma0')
# conf.gp$removeSamplers('gamma1')
conf.gp$addSampler(c('gamma0','s[1:575]'), type='AF_slice') 
MCMC.gp <- buildMCMC(conf.gp)
cMCMC.gp <- compileNimble(MCMC.gp)
results <- runMCMC(cMCMC.gp, nchain=1,niter = 200, thin=1, nburnin = 100,samplesAsCodaMCMC = T) 






# Prepare Data for Stan ----
distmat <- as.matrix(dist(rbind(dates.d[,c('X','Y')],st_coordinates(pred.locations.centroid))))


dat <- list(n = nrow(dates.d) + nrow(predictions.d),
	    n1 = nrow(dates.d),#Number of Samples
	    n2 = nrow(predictions.d), #Number of prediction locations
	    x  = c(dates.d$ash,predictions.d$ash), #Predictors
	    x2  = predictions.d$ash,  #Predictors (only predicted locations)
	    y  = dates.d$y, #Binary responses
	    dist  = distmat) # Distance Matrix (combined)

model2 <- stan_model(here('stan_scripts','model2.stan'))
fit2 <- sampling(model2,dat,iter=2000,chains=3)




