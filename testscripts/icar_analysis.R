library(nimbleCarbon)
library(sf)
library(rnaturalearth)
library(dplyr)
library(spdep)

# generate random locations in the main islands of japan
japan <- ne_countries(country = "Japan", returnclass = "sf",scale=10) |> st_transform(6684) |> st_geometry() |> st_cast('POLYGON') |> st_as_sf()
japan$area <- st_area(japan)
japan <- arrange(japan,desc(area)) |> slice_head(n = 4) |> st_geometry()
japan  <-  japan |> st_make_valid() |> st_union()
pred.locations <- st_make_grid(japan,cellsize=100000,square=FALSE)
pred.locations <- st_filter(st_as_sf(pred.locations),st_as_sf(japan),.predicate=st_intersects)
pred.locations <- pred.locations |> rename(geometry=x) |> mutate(id = row_number(),centroid = st_centroid(pred.locations))
nb_areas <- poly2nb(sf::as_Spatial(pred.locations), queen=FALSE, row.names = pred.locations$area_ID) #neighboring areas using sp library 
nbInfo <- nb2WB(nb_areas) #transform into iCAR inputs: adjacent matrix, weights, number of neighbors (for WinBUGS)
adj <- nbInfo$adj
weights <- nbInfo$weights
num <- nbInfo$num
L <- length(nbInfo$adj)




