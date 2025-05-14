# Load Library ----
library(here)
library(rcarbon)
library(nimbleCarbon)
library(rnaturalearth)
library(spdep)
library(sf)
library(terra)
library(dplyr)

# Filter Parameters ----
# From https://doi.org/10.1016/j.quascirev.2013.01.026
# 7165–7303 cal BP 
delta0 <- 7250
delta <- 750
lower <- delta0+delta
upper <- delta0-delta
binsize <- 100
hexgridsize  <- 125000 #125km
buffer <- 1000 #1km coastal buffer
ash.threshold  <- 200 #200 mm
max.buffer  <- 750*1000 #750km buffer for eastern edge
k.ah.coord  <- c(130.308,30.789)
# east.long.threshold  <- 1150000 #Use distance from Eruption point instead?

# Read ashfall data ----
tephra  <- rast(here('data','K_Ah_SUM.tif'))
tephra.4326 <- project(tephra,'EPSG:4326')

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
# extract ashfall data
hexgrid.latlong.vect <- st_transform(hexgrid,4326) |> vect()
avg.tephra  <- extract(tephra.4326,hexgrid.latlong.vect,fun=mean,na.rm=TRUE)
hexgrid$ash <- avg.tephra[,2]

# extract neighbourhood data
nb_areas <- poly2nb(sf::as_Spatial(hexgrid), queen=FALSE, row.names = hexgrid$area_ID) #neighboring areas using sp library 
nbInfo <- nb2WB(nb_areas) #transform into iCAR inputs: adjacent matrix, weights, number of neighbors (for WinBUGS)


## Create a subregion based on ash fall threshold
# Extract contour for pre-defined ashfall
cont <- as.contour(tephra,levels=ash.threshold) |> st_as_sf() |> st_transform(6684) |> st_cast('POLYGON')
k.ah.location  <- st_sfc(st_point(c(k.ah.coord)),crs=4326) |> st_transform(6684)
k.ah.buffer <- st_buffer(k.ah.location,dist=max.buffer) |> st_cast('MULTIPOLYGON')

# Create binary value
japan.sub <- st_intersection(japan.buffer,k.ah.buffer)
overlap <- st_intersection(japan.sub,cont) |> st_as_sf()
overlap$z <- 1
dif <- st_difference(japan.sub,st_union(cont)) |> st_cast('MULTIPOLYGON') |> st_as_sf()
dif$z <- 0
japan.sub.sf <- rbind(overlap,dif)



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
d  <- subset(d,C14Age < (uncalibrate(lower)$ccCRA+500) & C14Age > (uncalibrate(upper)$ccCRA-500))
cald <- calibrate(d$C14Age,d$C14Error,calMatrix=T)
i <- which.CalDates(cald,BP<lower&BP>upper,p=0.5)
d <- d[i,]
cald <- cald[i]
# Random thinning
bins <- binPrep(ages=cald,h=binsize,sites=d$SiteID)
j <- thinDates(ages=d$C14Age,errors=d$C14Error,bins=bins,size=1,thresh=1,seed=123)
d <- d[j,]

# Filter to Sites located within japan_buffer
i <- which(st_within(st_as_sf(d,coords=c('Easting','Northing'),crs=6684),japan.buffer,sparse=FALSE))
d <- d[i,]

# Rename SiteID 
d$SiteID <- as.numeric(as.factor(d$SiteID))

# Separate Site Table and Date Table ----
sites.df <- select(d,SiteID,Easting,Northing) |> unique() |> arrange(SiteID)
dates.df <- select(d,SiteID,C14Age,C14Error)
dates.df$ID <- 1:nrow(dates.df)
# Extract hexid from sites
sites.df$hexid <- st_join(st_as_sf(sites.df,coords=c('Easting','Northing'),crs=6684),hexgrid)$hexid
# Extract z from sites
sites.df$z <- st_join(st_as_sf(sites.df,coords=c('Easting','Northing'),crs=6684),japan.sub.sf)$z
# Join Table to dates.df
dates.df <- left_join(dates.df,sites.df)

# Define Lists for nimble analyses ----

# ICAR Models 
data(intcal20)
dat.icar <- list()
dat.icar$cra <- dates.df$C14Age
dat.icar$cra.error <- dates.df$C14Error

constants.icar <- list()
constants.icar$N  <- nrow(dates.df) #Number of Dates
constants.icar$id.hex <- dates.df$hexid #hexid of each sample
constants.icar$n.hex  <- max(hexgrid$hexid) #number of hex grid
#adjacency info for icar
constants.icar$adj  <- nbInfo$adj
constants.icar$weights <- nbInfo$weights
constants.icar$num <- nbInfo$num
constants.icar$L <- length(nbInfo$adj)
# Ash data
constants.icar$ash <- hexgrid$ash
#calibration vectors
constants.icar$calBP <- intcal20$CalBP
constants.icar$C14BP <- intcal20$C14Age
constants.icar$C14err <- intcal20$C14Age.sigma
#fixed parameters
constants.icar$a  <- lower
constants.icar$b  <- upper
constants.icar$delta0 <- delta0

mcal <- medCal(calibrate(dat.icar$cra,dat.icar$cra.error))
mcal <- ifelse(mcal>=lower,lower-1,mcal)
mcal <- ifelse(mcal<=upper,upper+1,mcal)
theta.init.icar <- mcal

# Two region comparison 
dates.df2 <- subset(dates.df,!is.na(z))

dat.two <- list()
dat.two$cra <- dates.df2$C14Age
dat.two$cra.error <- dates.df2$C14Error

constants.two <- list()
constants.two$N  <- nrow(dates.df2) #Number of Dates
constants.two$z <- dates.df2$z
#calibration vectors
constants.two$calBP <- intcal20$CalBP
constants.two$C14BP <- intcal20$C14Age
constants.two$C14err <- intcal20$C14Age.sigma
#fixed parameters
constants.two$a  <- lower
constants.two$b  <- upper
constants.two$delta0 <- delta0

mcal <- medCal(calibrate(dat.two$cra,dat.two$cra.error))
mcal <- ifelse(mcal>=lower,lower-1,mcal)
mcal <- ifelse(mcal<=upper,upper+1,mcal)
theta.init.two <- mcal

# Save everying into an R image file ----
save(dat.icar,dat.two,constants.icar,constants.two,theta.init.icar,theta.init.two,sites.df,dates.df,hexgrid,hexgrid_plot,file=here('01_dataprep_out.RData'))
