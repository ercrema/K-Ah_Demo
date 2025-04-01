# Load Library ----
library(here)
library(rcarbon)

# Filter Parameters ----
# From https://doi.org/10.1016/j.quascirev.2013.01.026
# 7165–7303 cal BP 
delta0 <- 7250
delta <- 750
lower <- delta0+delta
upper <- delta0-delta
binsize <- 100

# Read 14C dates ----
# Read Rekihaku Radiocarbon Database, English Curated Version, v.1.2.0 from https://github.com/ercrema/japan_c14db
d <- readRDS(url('https://github.com/ercrema/japan_c14db/raw/refs/heads/master/output/c14db_1.2.0.Rds'))

# Preprocessing & Filtering ----
# Exclude sites with no coordinates
d <- subset(d,!is.na(Latitude)&!is.na(Longitude))
# Aggregate sites based on proximity/dbscan
source(url('https://github.com/ercrema/yayoi_demo/raw/refs/heads/main/src/dbscanID.R'))
d$SiteID  <- dbscanID(sitename=d$SiteNameEn,longitude = d$Longitude,latitude = d$Latitude,eps=100)
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
cald <- cald[j]
# Compute probability mass before delta0
d$p_before <- apply(cald$calmatrix[1:which(rownames(cald$calmatrix)==as.character(delta0)),],2,sum)
# Extract Center point for each SiteID
SiteIDUniques <- unique(d$SiteID)
d$Latitude2 <- NA
d$Longitude2 <- NA
for (k in SiteIDUniques)
{
	tmp.i <- which(d$SiteID==k)
	lats  <- mean(unique(d$Latitude[tmp.i]))
	longs <- mean(unique(d$Longitude[tmp.i]))
	d$Latitude2[tmp.i] <- lats
	d$Longitude2[tmp.i] <- longs
}

# Extract Ash Levels from Raster File
d$ash <- NA
# <to do>

#






