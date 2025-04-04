library(here)
library(rstan)
library(dplyr)
library(sf)
library(rnaturalearth)

# Compile GP simulation in Stan
gpsim <- stan_model(here('stan_scripts','gpsim.stan'))

# sample settings
set.seed(123)
n.locations <- 200
n.samples <- 300
perc_misclass <- 0.2

# parameters
gamma0  <- -0.3
gamma1  <- 1.6
rho  <- 500
alpha  <- 0.5
sigma  <- 0.00001

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
simdat <- list(N = nrow(sites.d),
	    coords = as.matrix(sites.d[,c('X','Y')]/1000),#convert in KM
	    x = sites.d$ash,
	    rho=rho,
	    alpha=alpha,
	    sigma=sigma,
	    gamma0=gamma0,
	    gamma1=gamma1)

sim.out  <- sampling(gpsim,simdat,iter=1,warmup=0,chains=1,algorithm='Fixed_param') |> extract()
sites.d$p <- plogis(as.numeric(sim.out$lambda))

# Simulate binary dates 
dates.d <- data.frame(site.id=c(1:n.locations,sample(n.locations,size=n.samples-n.locations,replace=T)))
dates.d <- left_join(dates.d,sites.d)
dates.d$y <- rbinom(n=nrow(dates.d),size=1,prob=dates.d$p)
dates.d$p_y<- dates.d$y
i <- sample(nrow(dates.d),size=round(n.samples*perc_misclass))
dates.d$p_y[i] <- dates.d$p_y[i] + runif(round(n.samples*perc_misclass),min=-0.3,max=0.3)
dates.d$p_y[which(dates.d$p_y<0)]=0
dates.d$p_y[which(dates.d$p_y>1)]=1


# Save Everything
save(dates.d,sites.d,predictions.d,pred.locations,file=here('testscripts','simdata.RData'))

