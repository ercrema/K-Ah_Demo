# Load Library ----
library(here)
library(rcarbon)
library(nimbleCarbon)
library(rnaturalearth)
library(spdep)
library(sf)
library(terra)
library(dplyr)
data(intcal20)

# Filter Parameters ----
# From https://doi.org/10.1016/j.quascirev.2013.01.026
# 7165–7303 cal BP 
delta0 <- 7250
delta <- c(500,750,1000)
lower <- delta0+delta
upper <- delta0-delta
binsize <- 100
hexgridsize  <- 125000 #125km
buffer <- 1000 #1km coastal buffer
ash.threshold  <- 200 #200 mm
max.buffer  <- 750*1000 #750km buffer for eastern edge
k.ah.coord  <- c(130.308,30.789)

# Read ashfall data ----
tephra  <- rast(here('data','K_Ah_SUM.tif'))
tephra.4326 <- project(tephra,'EPSG:4326')

# Create Distance Raster File from Caldera ----
caldera.pts  <- vect(matrix(k.ah.coord,ncol=2),type='points',crs='EPSG:4326')
distance.rast <- distance(tephra.4326,caldera.pts) / 1000 #convert to km

# Prepare hexgrid ----
japan <- ne_countries(country = "Japan", returnclass = "sf",scale=10) |> st_transform(6684) |> st_geometry() |> st_cast('POLYGON') |> st_as_sf()
japan$area <- st_area(japan)
japan <- arrange(japan,desc(area)) |> slice_head(n = 4) |> st_geometry()
japan  <-  japan |> st_make_valid() |> st_union()
japan.buffer  <- st_buffer(japan,buffer)
hexgrid <- st_make_grid(japan,cellsize=hexgridsize,square=FALSE)
hexgrid <- st_filter(st_as_sf(hexgrid),st_as_sf(japan.buffer),.predicate=st_intersects)
hexgrid <- hexgrid |> rename(geometry=x) |> mutate(hexid = row_number(),centroid = st_centroid(hexgrid))

# make a clipped version for plotting purposes
hexgrid_plot <- st_intersection(hexgrid,japan)
# extract ashfall data and distance from caldera
hexgrid.latlong.vect <- st_transform(hexgrid,4326) |> vect()
avg.tephra  <- extract(tephra.4326,hexgrid.latlong.vect,fun=mean,na.rm=TRUE)
avg.dist  <- extract(distance.rast,hexgrid.latlong.vect,fun=mean,na.rm=TRUE)
hexgrid$ash <- avg.tephra[,2]
hexgrid$dist <- avg.dist[,2]

# extract neighbourhood data
nb_areas <- poly2nb(sf::as_Spatial(hexgrid), queen=FALSE, row.names = hexgrid$area_ID) #neighboring areas using sp library 
nbInfo <- nb2WB(nb_areas) #transform into iCAR inputs: adjacent matrix, weights, number of neighbors (for WinBUGS)

# Read 14C dates ----
# Read Rekihaku Radiocarbon Database, English Curated Version, v.1.2.0 from https://github.com/ercrema/japan_c14db
d <- readRDS(url('https://github.com/ercrema/japan_c14db/raw/refs/heads/master/output/c14db_1.2.0.Rds'))

# Preprocessing & Filtering ----
# Exclude sites with no coordinates
d <- subset(d,!is.na(Latitude)&!is.na(Longitude)&Prefecture!='Okinawa')
# Aggregate sites based on proximity/dbscan
source(url('https://github.com/ercrema/yayoi_demo/raw/refs/heads/main/src/dbscanID.R'))
d$SiteID  <- dbscanID(sitename=d$SiteNameEn,longitude = d$Longitude,latitude = d$Latitude,eps=100) #Assign temporary SiteID
# Compute average coordinate per site cluster
d$cLongitude <- d$cLatitude <- NA

for (i in 1:max(d$SiteID))
{
	tmp <- subset(d,SiteID==i)
	j <- which(d$SiteID==i)
	tmp2 <- select(tmp,SiteNameEn,Longitude,Latitude) |> unique()
	d$cLongitude[j] <- mean(tmp2$Longitude)
	d$cLatitude[j] <- mean(tmp2$Latitude)
}

# Extract Easting and Northing
sf_points_d <- st_as_sf(d, coords = c("cLongitude", "cLatitude"), crs = 4326)
sf_points_d_trans <- st_transform(sf_points_d,crs=6684)
coords <- st_coordinates(sf_points_d_trans)
d$Easting  <- coords[,'X']
d$Northing <- coords[,'Y']

# Handle Dates
d$C14Age = d$UnroundedCRA
i = which(is.na(d$C14Age))
d$C14Age[i] = d$CRA[i]
d$C14Error = d$UnroundedCRAError
i = which(is.na(d$C14Error))
d$C14Error[i] = d$CRAError[i]
d  <- subset(d,!is.na(C14Age) & !is.na(C14Error))
# Consider only terrestrial dates
d  <- subset(d,Material=='Terrestrial')
# Filter to wide range, and calibrate, and consider only dates within delta range
d  <- subset(d,C14Age < (uncalibrate(max(lower))$ccCRA+500) & C14Age > (uncalibrate(min(upper))$ccCRA-500))
cald <- calibrate(d$C14Age,d$C14Error,calMatrix=T)

# Define three sets of dates within 500,750,1000 years of delta
d500.i  <- which.CalDates(cald,BP<lower[1]&BP>upper[1],p=0.5)
d750.i  <- which.CalDates(cald,BP<lower[2]&BP>upper[2],p=0.5)
d1000.i <- which.CalDates(cald,BP<lower[3]&BP>upper[3],p=0.5)


# delta 500 ----
d.d500 <- d[d500.i,]
cald.d500 <- cald[d500.i]
# Random thinning
bins.d500 <- binPrep(ages=cald.d500,h=binsize,sites=d.d500$SiteID)
j.d500 <- thinDates(ages=d.d500$C14Age,errors=d.d500$C14Error,bins=bins.d500,size=1,thresh=1,seed=123)
d.d500 <- d[j.d500,]

# Filter to Sites located within japan_buffer
i.d500 <- which(st_within(st_as_sf(d.d500,coords=c('Easting','Northing'),crs=6684),japan.buffer,sparse=FALSE))
d.d500 <- d.d500[i.d500,]

# Rename SiteID 
d.d500$SiteID <- as.numeric(as.factor(d.d500$SiteID))

# Separate Site Table and Date Table 
sites.df.d500 <- select(d.d500,SiteID,Easting,Northing) |> unique() |> arrange(SiteID)
dates.df.d500 <- select(d.d500,SiteID,C14Age,C14Error)
dates.df.d500$ID <- 1:nrow(dates.df.d500)
# Extract hexid from sites
sites.df.d500$hexid <- st_join(st_as_sf(sites.df.d500,coords=c('Easting','Northing'),crs=6684),hexgrid)$hexid
# Join Table to dates.df
dates.df.d500 <- left_join(dates.df.d500,sites.df.d500)

# Define Lists for nimble analyses
# ICAR Models 
dat.icar.d500 <- list()
dat.icar.d500$cra <- dates.df.d500$C14Age
dat.icar.d500$cra.error <- dates.df.d500$C14Error

constants.icar.d500 <- list()
constants.icar.d500$N  <- nrow(dates.df.d500) #Number of Dates
constants.icar.d500$id.hex <- dates.df.d500$hexid #hexid of each sample
constants.icar.d500$n.hex  <- max(hexgrid$hexid) #number of hex grid
#adjacency info for icar
constants.icar.d500$adj  <- nbInfo$adj
constants.icar.d500$weights <- nbInfo$weights
constants.icar.d500$num <- nbInfo$num
constants.icar.d500$L <- length(nbInfo$adj)
# Predictors
constants.icar.d500$ash <- hexgrid$ash
constants.icar.d500$dist <- hexgrid$dist
#calibration vectors
constants.icar.d500$calBP <- intcal20$CalBP
constants.icar.d500$C14BP <- intcal20$C14Age
constants.icar.d500$C14err <- intcal20$C14Age.sigma
#fixed parameters
constants.icar.d500$a  <- lower[1]
constants.icar.d500$b  <- upper[1]
constants.icar.d500$delta0 <- delta0

mcal.d500 <- medCal(calibrate(dat.icar.d500$cra,dat.icar.d500$cra.error))
mcal.d500 <- ifelse(mcal.d500>=lower[1],lower[1]-1,mcal.d500)
mcal.d500 <- ifelse(mcal.d500<=upper[1],upper[1]+1,mcal.d500)
theta.init.icar.d500 <- mcal.d500

# Save everying into an R image file 
save(dat.icar.d500,constants.icar.d500,theta.init.icar.d500,sites.df.d500,dates.df.d500,hexgrid,hexgrid_plot,file=here('01_dataprep_out_d500.RData'))



# delta 750 ----
d.d750 <- d[d750.i,]
cald.d750 <- cald[d750.i]
# Random thinning
bins.d750 <- binPrep(ages=cald.d750,h=binsize,sites=d.d750$SiteID)
j.d750 <- thinDates(ages=d.d750$C14Age,errors=d.d750$C14Error,bins=bins.d750,size=1,thresh=1,seed=123)
d.d750 <- d[j.d750,]

# Filter to Sites located within japan_buffer
i.d750 <- which(st_within(st_as_sf(d.d750,coords=c('Easting','Northing'),crs=6684),japan.buffer,sparse=FALSE))
d.d750 <- d.d750[i.d750,]

# Rename SiteID 
d.d750$SiteID <- as.numeric(as.factor(d.d750$SiteID))

# Separate Site Table and Date Table 
sites.df.d750 <- select(d.d750,SiteID,Easting,Northing) |> unique() |> arrange(SiteID)
dates.df.d750 <- select(d.d750,SiteID,C14Age,C14Error)
dates.df.d750$ID <- 1:nrow(dates.df.d750)
# Extract hexid from sites
sites.df.d750$hexid <- st_join(st_as_sf(sites.df.d750,coords=c('Easting','Northing'),crs=6684),hexgrid)$hexid
# Join Table to dates.df
dates.df.d750 <- left_join(dates.df.d750,sites.df.d750)

# Define Lists for nimble analyses
# ICAR Models 
dat.icar.d750 <- list()
dat.icar.d750$cra <- dates.df.d750$C14Age
dat.icar.d750$cra.error <- dates.df.d750$C14Error

constants.icar.d750 <- list()
constants.icar.d750$N  <- nrow(dates.df.d750) #Number of Dates
constants.icar.d750$id.hex <- dates.df.d750$hexid #hexid of each sample
constants.icar.d750$n.hex  <- max(hexgrid$hexid) #number of hex grid
#adjacency info for icar
constants.icar.d750$adj  <- nbInfo$adj
constants.icar.d750$weights <- nbInfo$weights
constants.icar.d750$num <- nbInfo$num
constants.icar.d750$L <- length(nbInfo$adj)
# Predictors
constants.icar.d750$ash <- hexgrid$ash
constants.icar.d750$dist <- hexgrid$dist
#calibration vectors
constants.icar.d750$calBP <- intcal20$CalBP
constants.icar.d750$C14BP <- intcal20$C14Age
constants.icar.d750$C14err <- intcal20$C14Age.sigma
#fixed parameters
constants.icar.d750$a  <- lower[2]
constants.icar.d750$b  <- upper[2]
constants.icar.d750$delta0 <- delta0

mcal.d750 <- medCal(calibrate(dat.icar.d750$cra,dat.icar.d750$cra.error))
mcal.d750 <- ifelse(mcal.d750>=lower[2],lower[2]-1,mcal.d750)
mcal.d750 <- ifelse(mcal.d750<=upper[2],upper[2]+1,mcal.d750)
theta.init.icar.d750 <- mcal.d750

# Save everying into an R image file 
save(dat.icar.d750,constants.icar.d750,theta.init.icar.d750,sites.df.d750,dates.df.d750,hexgrid,hexgrid_plot,file=here('01_dataprep_out_d750.RData'))


# delta 1000 ----
d.d1000 <- d[d1000.i,]
cald.d1000 <- cald[d1000.i]
# Random thinning
bins.d1000 <- binPrep(ages=cald.d1000,h=binsize,sites=d.d1000$SiteID)
j.d1000 <- thinDates(ages=d.d1000$C14Age,errors=d.d1000$C14Error,bins=bins.d1000,size=1,thresh=1,seed=123)
d.d1000 <- d[j.d1000,]

# Filter to Sites located within japan_buffer
i.d1000 <- which(st_within(st_as_sf(d.d1000,coords=c('Easting','Northing'),crs=6684),japan.buffer,sparse=FALSE))
d.d1000 <- d.d1000[i.d1000,]

# Rename SiteID 
d.d1000$SiteID <- as.numeric(as.factor(d.d1000$SiteID))

# Separate Site Table and Date Table 
sites.df.d1000 <- select(d.d1000,SiteID,Easting,Northing) |> unique() |> arrange(SiteID)
dates.df.d1000 <- select(d.d1000,SiteID,C14Age,C14Error)
dates.df.d1000$ID <- 1:nrow(dates.df.d1000)
# Extract hexid from sites
sites.df.d1000$hexid <- st_join(st_as_sf(sites.df.d1000,coords=c('Easting','Northing'),crs=6684),hexgrid)$hexid
# Join Table to dates.df
dates.df.d1000 <- left_join(dates.df.d1000,sites.df.d1000)

# Define Lists for nimble analyses
# ICAR Models 
dat.icar.d1000 <- list()
dat.icar.d1000$cra <- dates.df.d1000$C14Age
dat.icar.d1000$cra.error <- dates.df.d1000$C14Error

constants.icar.d1000 <- list()
constants.icar.d1000$N  <- nrow(dates.df.d1000) #Number of Dates
constants.icar.d1000$id.hex <- dates.df.d1000$hexid #hexid of each sample
constants.icar.d1000$n.hex  <- max(hexgrid$hexid) #number of hex grid
#adjacency info for icar
constants.icar.d1000$adj  <- nbInfo$adj
constants.icar.d1000$weights <- nbInfo$weights
constants.icar.d1000$num <- nbInfo$num
constants.icar.d1000$L <- length(nbInfo$adj)
# Predictors
constants.icar.d1000$ash <- hexgrid$ash
constants.icar.d1000$dist <- hexgrid$dist
#calibration vectors
constants.icar.d1000$calBP <- intcal20$CalBP
constants.icar.d1000$C14BP <- intcal20$C14Age
constants.icar.d1000$C14err <- intcal20$C14Age.sigma
#fixed parameters
constants.icar.d1000$a  <- lower[3]
constants.icar.d1000$b  <- upper[3]
constants.icar.d1000$delta0 <- delta0

mcal.d1000 <- medCal(calibrate(dat.icar.d1000$cra,dat.icar.d1000$cra.error))
mcal.d1000 <- ifelse(mcal.d1000>=lower[3],lower[3]-1,mcal.d1000)
mcal.d1000 <- ifelse(mcal.d1000<=upper[3],upper[3]+1,mcal.d1000)
theta.init.icar.d1000 <- mcal.d1000

# Save everying into an R image file 
save(dat.icar.d1000,constants.icar.d1000,theta.init.icar.d1000,sites.df.d1000,dates.df.d1000,hexgrid,hexgrid_plot,file=here('01_dataprep_out_d1000.RData'))


