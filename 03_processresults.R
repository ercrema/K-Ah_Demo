# Load Libraries ----
library(here)
library(ggplot2)
library(sf)
library(coda)
library(dplyr)
library(rnaturalearth)
library(tidybayes)
library(latex2exp)
library(terra)
library(viridis)
library(gridExtra)
library(tidyr)
library(ggridges)

# Load Utility Functions ----
source(here("src", "ppscripts.R"))

# Load Data and Results ----
# Full 14C dates
d <- readRDS(url('https://github.com/ercrema/japan_c14db/raw/refs/heads/master/output/c14db_1.2.0.Rds'))
# Data
load(here('data','01_dataprep_out_d500.RData'))
load(here('data','01_dataprep_out_d750.RData'))
load(here('data','01_dataprep_out_d1000.RData'))
# Results (main - 750yrs)
load(here('results','02a_fitmodel0_750_out.RData'))
load(here('results','02a_fitmodel1_750_out.RData'))
load(here('results','02a_fitmodel2_750_out.RData'))
load(here('results','02a_fitmodel3_750_out.RData'))

load(here('results','02b_fitmodel0_750_out.RData'))
load(here('results','02b_fitmodel1_750_out.RData'))
load(here('results','02b_fitmodel2_750_out.RData'))
load(here('results','02b_fitmodel3_750_out.RData'))

# Results (sensitivity analysis)
load(here('results','02a_fitmodel2_500_out.RData'))
load(here('results','02a_fitmodel2_1000_out.RData'))
load(here('results','02b_fitmodel2_500_out.RData'))
load(here('results','02b_fitmodel2_1000_out.RData'))

# Sample Summaries ----

sample.sizes <- data.frame(Sample=c('original','d500','d750','d1000'),
			   N=c(nrow(d),constants.icar.d500$N,constants.icar.d750$N,constants.icar.d1000$N))
write.csv(sample.sizes,here('figures_and_tables','sample_sizes.csv'),row.names=FALSE)	

# WAIC-based model comparison ----

waic.table<- data.frame(Model=c('0','1','2','3'),
			WAIC.a=c(waic.a0.750$WAIC,
				 waic.a1.750$WAIC,
				 waic.a2.750$WAIC,
				 waic.a3.750$WAIC),
			dWAIC1.a=deltaWAIC(c(waic.a0.750$WAIC,
					    waic.a1.750$WAIC,
					    waic.a2.750$WAIC,
					    waic.a3.750$WAIC)),	
			dWAIC2.a=c(waic.a0.750$WAIC,
					    waic.a1.750$WAIC,
					    waic.a2.750$WAIC,
					    waic.a3.750$WAIC)-waic.a2.750$WAIC,	
			WAICweight.a=WAICweights(c(waic.a0.750$WAIC,
						   waic.a1.750$WAIC,
						   waic.a2.750$WAIC,
						   waic.a3.750$WAIC)),
			WAIC.b=c(waic.b0.750$WAIC,
				 waic.b1.750$WAIC,
				 waic.b2.750$WAIC,
				 waic.b3.750$WAIC),
			dWAIC1.b=deltaWAIC(c(waic.b0.750$WAIC,
					    waic.b1.750$WAIC,
					    waic.b2.750$WAIC,
					    waic.b3.750$WAIC)),
			dWAIC2.b=c(waic.b0.750$WAIC,
					    waic.b1.750$WAIC,
					    waic.b2.750$WAIC,
					    waic.b3.750$WAIC)-waic.b2.750$WAIC,
			WAICweight.b=WAICweights(c(waic.b0.750$WAIC,
						   waic.b1.750$WAIC,
						   waic.b2.750$WAIC,
						   waic.b3.750$WAIC)))

write.csv(waic.table,here('figures_and_tables','waic.csv'),row.names=FALSE)

# Post-processing for Model A ----

# Extract posteriors of r1 and r2 into matrices
r1.0.d750.matrix  <- posterior.a0.750[,grep('s1',colnames(posterior.a0.750))]
r2.0.d750.matrix  <- posterior.a0.750[,grep('s2',colnames(posterior.a0.750))]

r1.1.d750.matrix  <- posterior.a1.750[,grep('s1',colnames(posterior.a1.750))]
r2.1.d750.matrix  <- posterior.a1.750[,grep('s2',colnames(posterior.a1.750))] + posterior.a1.750[,'beta1'] %*% t(constants.icar.d750$dist)

r1.2.d750.matrix  <- posterior.a2.750[,grep('s1',colnames(posterior.a2.750))]
r2.2.d750.matrix  <- posterior.a2.750[,grep('s2',colnames(posterior.a2.750))] + posterior.a2.750[,'beta1'] %*% t(constants.icar.d750$logAsh) + posterior.a2.750[,'beta2'] * matrix(pmax(0,outer(constants.icar.d750$logAsh,posterior.a2.750[,'kappa'],'-')),ncol=61,byrow=T) #ncol is number hex

r1.3.d750.matrix  <- posterior.a3.750[,grep('s1',colnames(posterior.a3.750))]
r2.3.d750.matrix  <- posterior.a3.750[,grep('s2',colnames(posterior.a3.750))] + posterior.a3.750[,'beta1'] %*% t(constants.icar.d750$ash)


# Additional r1 and r2 for sensitivity analysis

r1.2.d500.matrix  <- posterior.a2.500[,grep('s1',colnames(posterior.a2.500))]
r2.2.d500.matrix  <- posterior.a2.500[,grep('s2',colnames(posterior.a2.500))] + posterior.a2.500[,'beta1'] %*% t(constants.icar.d500$logAsh) + posterior.a2.500[,'beta2'] * matrix(pmax(0,outer(constants.icar.d500$logAsh,posterior.a2.500[,'kappa'],'-')),ncol=61,byrow=T) #ncol is number hex

r1.2.d1000.matrix  <- posterior.a2.1000[,grep('s1',colnames(posterior.a2.1000))]
r2.2.d1000.matrix  <- posterior.a2.1000[,grep('s2',colnames(posterior.a2.1000))] + posterior.a2.1000[,'beta1'] %*% t(constants.icar.d1000$logAsh) + posterior.a2.1000[,'beta2'] * matrix(pmax(0,outer(constants.icar.d1000$logAsh,posterior.a2.1000[,'kappa'],'-')),ncol=61,byrow=T) #ncol is number hex


post.model.a0.d750 <- post_summary_a(r1=r1.0.d750.matrix,r2=r2.0.d750.matrix,n=constants.icar.d750$n.hex,p=0.9)
post.model.a1.d750 <- post_summary_a(r1=r1.1.d750.matrix,r2=r2.1.d750.matrix,n=constants.icar.d750$n.hex,p=0.9)
post.model.a2.d750 <- post_summary_a(r1=r1.2.d750.matrix,r2=r2.2.d750.matrix,n=constants.icar.d750$n.hex,p=0.9)
post.model.a3.d750 <- post_summary_a(r1=r1.3.d750.matrix,r2=r2.3.d750.matrix,n=constants.icar.d750$n.hex,p=0.9)


post.model.a2.d500 <- post_summary_a(r1=r1.2.d500.matrix,r2=r2.2.d500.matrix,n=constants.icar.d500$n.hex,p=0.9)
post.model.a2.d1000 <- post_summary_a(r1=r1.2.d1000.matrix,r2=r2.2.d1000.matrix,n=constants.icar.d1000$n.hex,p=0.9)

# Post-processing for Model B ----
# Extract posteriors of eta into matrices
eta.0.d750.matrix <- posterior.b0.750[,grep('^eta',colnames(posterior.b0.750))] 
eta.1.d750.matrix <- posterior.b1.750[,grep('^eta',colnames(posterior.b1.750))] + posterior.b1.750[,'beta1'] %*% t(constants.icar.d750$dist)
eta.2.d750.matrix <- posterior.b2.750[,grep('^eta',colnames(posterior.b2.750))] + posterior.b2.750[,'beta1'] %*% t(constants.icar.d750$logAsh) + posterior.b2.750[,'beta2'] * matrix(pmax(0,outer(constants.icar.d750$logAsh,posterior.b2.750[,'kappa'],'-')),ncol=61,byrow=T) #ncol is number hex
eta.3.d750.matrix <- posterior.b3.750[,grep('^eta',colnames(posterior.b3.750))] + posterior.b3.750[,'beta1'] %*% t(constants.icar.d750$ash)

# Additional eta for sensitivity analysis
eta.2.d500.matrix <- posterior.b2.500[,grep('^eta',colnames(posterior.b2.500))] + posterior.b2.500[,'beta1'] %*% t(constants.icar.d500$logAsh) + posterior.b2.500[,'beta2'] * matrix(pmax(0,outer(constants.icar.d500$logAsh,posterior.b2.500[,'kappa'],'-')),ncol=61,byrow=T) #ncol is number hex
eta.2.d1000.matrix <- posterior.b2.1000[,grep('^eta',colnames(posterior.b2.1000))] + posterior.b2.1000[,'beta1'] %*% t(constants.icar.d1000$logAsh) + posterior.b2.1000[,'beta2'] * matrix(pmax(0,outer(constants.icar.d1000$logAsh,posterior.b2.1000[,'kappa'],'-')),ncol=61,byrow=T) #ncol is number hex


# Aggregate summary stats
post.model.b0.d750 <- post_summary_b(eta=eta.0.d750.matrix,n=constants.icar.d750$n.hex,p=0.9)	
post.model.b1.d750 <- post_summary_b(eta=eta.1.d750.matrix,n=constants.icar.d750$n.hex,p=0.9)	
post.model.b2.d750 <- post_summary_b(eta=eta.2.d750.matrix,n=constants.icar.d750$n.hex,p=0.9)	
post.model.b3.d750 <- post_summary_b(eta=eta.3.d750.matrix,n=constants.icar.d750$n.hex,p=0.9)	


post.model.b2.d500 <- post_summary_b(eta=eta.2.d500.matrix,n=constants.icar.d500$n.hex,p=0.9)	
post.model.b2.d1000 <- post_summary_b(eta=eta.2.d1000.matrix,n=constants.icar.d1000$n.hex,p=0.9)	

# Extract posteriors from key regions for ICAR models
keys <- c(1,4,14,16,21,49,50)
key.letters  <- letters[1:length(keys)]
key_hexgrid <- st_geometry(hexgrid_plot[keys,]) |> st_cast('MULTILINESTRING')
r1r2keys  <- cbind(r1.2.d750.matrix[,paste0('s1[',keys,']')]*100,r2.2.d750.matrix[,paste0('s2[',keys,']')]*100)  
r1r2keys <- gather(r1r2keys)
r1r2keys$param <- 'r1'
r1r2keys$param[grep('s2',r1r2keys$key)] <- 'r2'
r1r2keys$hex <- key.letters[match(as.integer(sub(".*\\[(\\d+)\\].*", "\\1", r1r2keys$key)),keys)]

etakeys  <- eta.2.d750.matrix[,paste0('eta[',keys,']')] 
etakeys <- gather(etakeys)
etakeys$hex <- key.letters[match(as.integer(sub(".*\\[(\\d+)\\].*", "\\1", etakeys$key)),keys)]

# Map Figure Preparation ----

# combine to plot grid
icar.a.hexgrid.d500 <- left_join(hexgrid_plot,post.model.a2.d500[[1]])
icar.a.hexgrid.d750 <- left_join(hexgrid_plot,post.model.a2.d750[[1]])
icar.a.hexgrid.d1000 <- left_join(hexgrid_plot,post.model.a2.d1000[[1]])

icar.b.hexgrid.d500 <- left_join(hexgrid_plot,post.model.b2.d500[[1]])
icar.b.hexgrid.d750 <- left_join(hexgrid_plot,post.model.b2.d750[[1]])
icar.b.hexgrid.d1000 <- left_join(hexgrid_plot,post.model.b2.d1000[[1]])

# extract sample sizes
mt500 <- st_intersects(hexgrid_plot,st_as_sf(sites.df.d500,coords=c('Easting','Northing'),crs=6684))
mt750 <- st_intersects(hexgrid_plot,st_as_sf(sites.df.d750,coords=c('Easting','Northing'),crs=6684))
mt1000 <- st_intersects(hexgrid_plot,st_as_sf(sites.df.d1000,coords=c('Easting','Northing'),crs=6684))
icar.a.hexgrid.d500$obs  <- icar.b.hexgrid.d500$obs <- lengths(mt500) > 0
icar.a.hexgrid.d750$obs  <- icar.b.hexgrid.d750$obs <- lengths(mt750) > 0
icar.a.hexgrid.d1000$obs  <- icar.b.hexgrid.d1000$obs <- lengths(mt1000) > 0

#convert to lat/long
icar.a.hexgrid.d500 <- st_transform(icar.a.hexgrid.d500,crs=4326)
icar.a.hexgrid.d750 <- st_transform(icar.a.hexgrid.d750,crs=4326)
icar.a.hexgrid.d1000 <- st_transform(icar.a.hexgrid.d1000,crs=4326)	
icar.b.hexgrid.d500 <- st_transform(icar.b.hexgrid.d500,crs=4326)
icar.b.hexgrid.d750 <- st_transform(icar.b.hexgrid.d750,crs=4326)
icar.b.hexgrid.d1000 <- st_transform(icar.b.hexgrid.d1000,crs=4326)	

sites.d500 <- st_as_sf(sites.df.d500,coords=c('Easting','Northing'),crs=6684) |> st_transform(crs=4326)
sites.d750 <- st_as_sf(sites.df.d750,coords=c('Easting','Northing'),crs=6684) |> st_transform(crs=4326)
sites.d1000 <- st_as_sf(sites.df.d1000,coords=c('Easting','Northing'),crs=6684) |> st_transform(crs=4326)

# Download background map
win  <- ne_countries(scale=10,returnclass='sf') |> st_combine() #download background map

# Read ashfall data
tephra  <- rast(here('data','K_Ah_SUM.tif'))
tephra <- project(tephra,"EPSG:4326")
tephra <- as.data.frame(tephra,xy=T)
tephra[tephra==0] <- NA
colnames(tephra)[3] <- 'value'
hexgrid_plot$ash <- hexgrid$ash

# Make Figures ----

# Sample Distribution Map and Eruption Model  ----
sample_map <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour='darkgrey') +
	geom_sf(data=hexgrid_plot,aes(fill=ash),color='white',alpha=0.7)+
	scale_fill_viridis(option='turbo',trans='log10') +
	geom_sf(data=key_hexgrid,color='black') +
	geom_sf(data=st_geometry(sites.d750),size=0.5,col='grey20',alpha=0.8) +
#       geom_contour(data=tephra,aes(x=x,y=y,z=value),linetype='dashed',breaks=200,color='black') +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('point',x=130.308,y=30.789,fill='red',size=2,shape=24) +
	annotate('text',x=133.5,y=30.9,label='Kikai Caldera',size=3) +
	annotate('text',x=129.5,y=32.5,label='a') +
	annotate('text',x=129.5,y=31,label='b') +
	annotate('text',x=132,y=36,label='c') +
	annotate('segment',x=132.3,xend=134,y=35.8,yend=35,colour='black') +
	annotate('text',x=134.7,y=36.3,label='d') +
	annotate('text',x=137,y=34,label='e') +
	annotate('text',x=142,y=37.5,label='f') +
	annotate('text',x=142.5,y=39.5,label='g') +
	labs(x='Longitude',y='Latitude',fill='Average Deposit \nThickness (mm)') +
	theme(legend.position='inside',legend.position.inside=c(0.2,0.75),legend.background=element_rect(fill=alpha('white',0.5)),legend.key.size=unit(0.2,'in'),legend.text=element_text(size=5.5),legend.title=element_text(size=7))

ashfall_model_raw <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour='darkgrey') +
#	geom_contour_filled(data=tephra,aes(x=x,y=y,z=value),alpha=0.5,breaks=c(0,150,300,500,1000,1500,2000,2500,3000)) +
	geom_contour_filled(data=tephra,aes(x=x,y=y,z=value),breaks=c(0,5,25,150,1000,2000,3500)) +
	scale_fill_viridis_d(
			     option = "plasma",  # or "viridis", "inferno", etc.
			     direction = 1,      # 1 = low to high
			     alpha = 0.5,        # semitransparency
			     labels = c(" < 5mm (Thin)", "5 - 25 mm (Modetate)", "25 - 150mm (Thick)", "150 - 1000mm (Very Thick)", "1000 - 2000 mm ", "> 2000 mm "),
			     name = "Tephra Thickness (in mm)"
	) + 	 
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('point',x=130.308,y=30.789,fill='red',size=2,shape=24) +
	annotate('text',x=133.5,y=30.9,label='Kikai Caldera',size=3) +
	annotate('text',x=131,y=33,label='Kyushu',size=2.5) +
	annotate('text',x=134,y=33.7,label='Shikoku',size=2.5) +
	annotate('text',x=139,y=36.5,label='Honshu',size=2.5) +
	annotate('text',x=143,y=43.5,label='Hokkaido',size=2.5) +
	labs(fill='Tephra thickness (in mm)',x='Longitude',y='Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.75),legend.background=element_rect(fill=alpha('white',0.5)),legend.key.size=unit(0.2,'in'),legend.text=element_text(size=5.5),legend.title=element_text(size=7))

ashandhex <- grid.arrange(sample_map,ashfall_model_raw,ncol=2,padding=unit(0.2,'cm'))
ggsave(here('figures_and_tables','samplemap_and_ashfall.pdf'),plot=ashandhex,width=7,height=4)


# Sample Distribution Map for d500 and d750  ----
sample_map500 <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour='darkgrey') +
	geom_sf(data=hexgrid_plot,fill='lightgrey',color='black',alpha=0.5)+
	geom_sf(data=st_geometry(sites.d500),size=0.5,col='grey20',alpha=0.8) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',x=130,y=45,label='a',size=6) +
	labs(x='Longitude',y='Latitude')

sample_map1000 <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour='darkgrey') +
	geom_sf(data=hexgrid_plot,fill='lightgrey',color='black',alpha=0.5)+
	geom_sf(data=st_geometry(sites.d1000),size=0.5,col='grey20',alpha=0.8) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',x=130,y=45,label='b',size=6) +
	labs(x='Longitude',y='Latitude')

samplemaps <- grid.arrange(sample_map500,sample_map1000,ncol=2,padding=unit(0.2,'cm'))
ggsave(here('figures_and_tables','samplemap_500_and_1000.pdf'),plot=samplemaps,width=7,height=4)



# r1 posterior model A d750----


r1plot.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d750,aes(fill=r1.p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=post.model.a2.d750$minmax.lohi) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX('Posterior Mean $r_1$',output='character'),x=132,y=45,parse=T) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())

r1plot.lo <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d750,aes(fill=r1.p_lo),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=post.model.a2.d750$minmax.lohi) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(lower)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 

r1plot.hi <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d750,aes(fill=r1.p_hi),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=post.model.a2.d750$minmax.lohi) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(lower)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 

r1plot <- grid.arrange(r1plot.main,r1plot.lo,r1plot.hi,layout_matrix = rbind(c(1,2),c(1,3)),widths=c(2,1),heights=c(1,1),padding=unit(0.2,'cm'))

ggsave(here('figures_and_tables','posterior_r1_d750.pdf'),plot=r1plot,width=7,height=5)

# r2 posterior model A d750 ----

r2plot.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d750,aes(fill=r2.p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=post.model.a2.d750$minmax.lohi) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX("Posterior Mean $r_2$",output='character'),x=132,y=45,parse=T) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())

r2plot.lo <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d750,aes(fill=r2.p_lo),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=post.model.a2.d750$minmax.lohi) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(lower)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 

r2plot.hi <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d750,aes(fill=r2.p_hi),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=post.model.a2.d750$minmax.lohi) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(higher)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 

r2plot <- grid.arrange(r2plot.main,r2plot.lo,r2plot.hi,layout_matrix = rbind(c(1,2),c(1,3)),widths=c(2,1),heights=c(1,1),padding=unit(0.2,'cm'))

ggsave(here('figures_and_tables','posterior_r2_d750.pdf'),plot=r2plot,width=7,height=5)


# Posterior eta model B d750----
etaplot.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.b.hexgrid.d750,aes(fill=p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.b.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name=TeX('$\\eta$'),limits=post.model.b2.d750$minmax.lohi) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX('Posterior Mean $\\eta$',output='character'),x=132,y=45,parse=TRUE) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())

etaplot.lo <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.b.hexgrid.d750,aes(fill=p_lo90),color='darkgrey') +
	geom_sf(data=subset(icar.b.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name=TeX('$\\eta$'),limits=post.model.b2.d750$minmax.lohi) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(lower)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 

etaplot.hi <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.b.hexgrid.d750,aes(fill=p_hi90),color='darkgrey') +
	geom_sf(data=subset(icar.b.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name=TeX('$\\eta$'),limits=post.model.b2.d750$minmax.lohi) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(lower)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 

etaplot <- grid.arrange(etaplot.main,etaplot.lo,etaplot.hi,layout_matrix = rbind(c(1,2),c(1,3)),widths=c(2,1),heights=c(1,1),padding=unit(0.2,'cm'))
ggsave(here('figures_and_tables','posterior_eta_d750.pdf'),plot=etaplot,width=7,height=5)


# Mean Posterior Sensitivity  ----
rangeminmax.mean.a  <- range(c(post.model.a.i.d500$minmax.mean,post.model.a.i.d750$minmax.mean,post.model.a.i.d1000$minmax.mean))
rangeminmax.mean.b  <- range(c(post.model.b.i.d500$minmax.mean,post.model.b.i.d750$minmax.mean,post.model.b.i.d1000$minmax.mean))
r1plot.d500 <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d500,aes(fill=r1.p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d500,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=rangeminmax.mean.a) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX("Posterior Mean $r_1$",output='character'),x=132,y=45,parse=TRUE) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank()) +
	ggtitle('500yrs')
r1plot.d750 <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d750,aes(fill=r1.p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=rangeminmax.mean.a) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX("Posterior Mean $r_1$",output='character'),x=132,y=45,parse=TRUE) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank()) +
	ggtitle('750yrs')
r1plot.d1000 <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d1000,aes(fill=r1.p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d1000,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=rangeminmax.mean.a) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX("Posterior Mean $r_1$",output='character'),x=132,y=45,parse=TRUE) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank()) +
	ggtitle('1000yrs')

r2plot.d500 <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d500,aes(fill=r2.p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d500,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=rangeminmax.mean.a) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX("Posterior Mean $r_2$",output='character'),x=132,y=45,parse=TRUE) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank()) +
	ggtitle('500yrs')
r2plot.d750 <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d750,aes(fill=r2.p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=rangeminmax.mean.a) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX("Posterior Mean $r_2$",output='character'),x=132,y=45,parse=TRUE) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank()) +
	ggtitle('750yrs')
r2plot.d1000 <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d1000,aes(fill=r2.p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d1000,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=rangeminmax.mean.a) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX("Posterior Mean $r_2$",output='character'),x=132,y=45,parse=TRUE) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank()) +
	ggtitle('1000yrs')

etaplot.d500 <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.b.hexgrid.d500,aes(fill=p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.b.hexgrid.d500,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name=TeX('$\\eta$'),limits=rangeminmax.mean.b) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX('Posterior Mean $\\eta$',output='character'),x=132,y=45,parse=TRUE) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank()) +
	ggtitle('500yrs')

etaplot.d750 <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.b.hexgrid.d750,aes(fill=p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.b.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name=TeX('$\\eta$'),limits=rangeminmax.mean.b) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX('Posterior Mean $\\eta$',output='character'),x=132,y=45,parse=TRUE) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank()) +
	ggtitle('750yrs')

etaplot.d1000 <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.b.hexgrid.d1000,aes(fill=p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.b.hexgrid.d1000,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name=TeX('$\\eta$'),limits=rangeminmax.mean.b) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX('Posterior Mean $\\eta$',output='character'),x=132,y=45,parse=TRUE) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank()) +
	ggtitle('1000yrs')


sensitivityplot <- grid.arrange(r1plot.d500,r1plot.d750,r1plot.d1000,r2plot.d500,r2plot.d750,r2plot.d1000,etaplot.d500,etaplot.d750,etaplot.d1000,layout_matrix = rbind(c(1,2,3),c(4,5,6),c(7,8,9)),widths=c(1,1,1),heights=c(1,1,1),padding=unit(0.2,'cm'))

ggsave(here('figures_and_tables','mean_posterior_sensitivity.pdf'),plot=sensitivityplot,width=14,height=14)


# Mean Posterior Plots for r1, r2, r2-r1, and eta for 750 ----
r1plot.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d750,aes(fill=r1.p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=post.model.a2.d750$minmax.mean) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX('Posterior Mean $r_1$',output='character'),x=134,y=45,parse=T) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())

r2plot.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d750,aes(fill=r2.p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=post.model.a2.d750$minmax.mean) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX("Posterior Mean $r_2$",output='character'),x=134,y=45,parse=T) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())

etaplot.main <- ggplot(st_geometry(win)) + 
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.b.hexgrid.d750,aes(fill=p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.b.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name=TeX('$\\eta$'),limits=post.model.b2.d750$minmax.mean) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX('Posterior Mean $\\eta$',output='character'),x=134,y=45,parse=TRUE) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.15,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())

r2r1plot.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d750,aes(fill=dr2r1.p_mean),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='darkgreen',mid='white',high='orange',midpoint=0,name='Annual Growth Rate (%)') +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX("Posterior Mean $r_2$ - $r_1$",output='character'),x=134,y=45,parse=TRUE) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())

meanpost  <- grid.arrange(r1plot.main,r2plot.main,r2r1plot.main,etaplot.main,layout_matrix = rbind(c(1,2),c(3,4)),widths=c(1,1),heights=c(1,1),padding=unit(0.2,'cm'))
ggsave(here('figures_and_tables','mean_posterior.pdf'),plot=meanpost,width=7,height=7)

# Decrease Probability Plots for r1, r2, r2-r1, and eta for 750 ----

r1decplot.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d750,aes(fill=r1.p_decrease),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	#scale_fill_gradient2(low='darkgreen',mid='white',high='darkorange',midpoint=0.5,name='Probability',limits=c(0,1)) +
	scale_fill_gradient(low='white',high='steelblue',name='Probability',limits=c(0,1)) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX("P[$r_1 < 0]$",output='character'),parse=TRUE,x=132,y=45) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.2,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())

r2decplot.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d750,aes(fill=r2.p_decrease),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	#scale_fill_gradient2(low='darkgreen',mid='white',high='darkorange',midpoint=0.5,name='Probability',limits=c(0,1)) +
	scale_fill_gradient(low='white',high='steelblue',name='Probability',limits=c(0,1)) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX("P[$r_2 < 0]$",output='character'),parse=TRUE,x=132,y=45) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.2,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())


r2r1decplot.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.a.hexgrid.d750,aes(fill=dr2r1.p_decrease),color='darkgrey') +
	geom_sf(data=subset(icar.a.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	#scale_fill_gradient2(low='darkgreen',mid='white',high='darkorange',midpoint=0.5,name='Probability',limits=c(0,1)) +
	scale_fill_gradient(low='white',high='steelblue',name='Probability',limits=c(0,1)) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX("P[$r_2-r_1<0]$",output='character'),parse=TRUE,x=132,y=45) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.2,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())

etadecplot.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=icar.b.hexgrid.d750,aes(fill=p_decrease),color='darkgrey') +
	geom_sf(data=subset(icar.b.hexgrid.d750,obs==FALSE),fill='darkgrey',color='darkgrey') +
	#scale_fill_gradient2(low='darkgreen',mid='white',high='darkorange',midpoint=0.5,name='Probability',limits=c(0,1)) +
	scale_fill_gradient(low='white',high='steelblue',name='Probability',limits=c(0,1)) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX("P[$\\eta<0]$",output='character'),parse=TRUE,x=132,y=45) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.2,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())

decprob  <- grid.arrange(r1decplot.main,r2decplot.main,r2r1decplot.main,etadecplot.main,layout_matrix = rbind(c(1,2),c(3,4)),widths=c(1,1),heights=c(1,1),padding=unit(0.2,'cm'))
ggsave(here('figures_and_tables','decrease_probability.pdf'),plot=decprob,width=7,height=7)

# HexSpecific Plot ----
r1r2keyplot <- ggplot(r1r2keys,aes(x=value,col=param,y=hex)) +
	stat_pointinterval(position=position_dodge(0.3),.width=c(.5,.9),point_interval='median_hdi') +
	scale_color_manual(values=c('r1'='darkorange','r2'='darkgreen'),labels=c('r1'=TeX('$r_1$'),'r2'=TeX('$r_2$'))) +
	theme_minimal() +
	xlim(-0.5,0.5) +
	labs(x='Annual Growth Rate',y='Hexagon',col='Parameters') +
	geom_vline(xintercept=0,linetype='dashed') +
	theme(legend.position = 'inside', legend.position.inside=c(0.2,0.8))


etakeyplot <- ggplot(etakeys,aes(x=value,y=hex)) + 
	stat_pointinterval(position=position_dodge(0.3),.width=c(.5,.9),point_interval='median_hdi') +
	theme_minimal() +
	geom_vline(xintercept=0,linetype='dashed') +
	scale_x_continuous(name=TeX('$\\eta$'),sec.axis=sec_axis(transform=~1/(1+exp(-.)),breaks=round(plogis(c(-3,-2,-1,0,1)),2),name='Relative proportion of dates after the eruption')) +
	labs(x=TeX('$\\eta$'),y='Hexagon')


hex_combined <- grid.arrange(r1r2keyplot,etakeyplot,ncol=1)
ggsave(here('figures_and_tables','hex_focus_plot.pdf'),plot=hex_combined,width=4.75,height=7)

# Posterior parameters vs predictor ----
icar.a.hexgrid.d750$ash <- constants.icar.d750$ash
icar.a.hexgrid.d750$keys <- NA
icar.a.hexgrid.d750$keys[keys] <- key.letters

icar.b.hexgrid.d750$ash <- constants.icar.d750$ash
icar.b.hexgrid.d750$keys <- NA
icar.b.hexgrid.d750$keys[keys] <- key.letters

ashplot1 <- ggplot(data=icar.a.hexgrid.d750,aes(x=ash+1,y=r2.p_mean)) +
	geom_errorbar(aes(ymin=r2.p_lo,ymax=r2.p_hi),alpha=0.8,color='lightgrey') +
	geom_point() +
	geom_vline(xintercept=150,linetype='dashed') +
	geom_hline(yintercept=0,linetype='dotted') +
	scale_x_log10() +
	labs(x='Average Thickness of Tephra Fall Deposit',y='Posterior Estimate of Annual Growth Rate (%)') +
	theme_minimal()

ashplot2 <- ggplot(data=icar.b.hexgrid.d750,aes(x=ash+1,y=p_mean)) +
	geom_errorbar(aes(ymin=p_lo90,ymax=p_hi90),alpha=0.8,color='lightgrey') +
	geom_point() +
	geom_vline(xintercept=150,linetype='dashed') +
	geom_hline(yintercept=0,linetype='dotted') +
	scale_x_log10() +
	scale_y_continuous(name=TeX('$\\eta$'),sec.axis=sec_axis(trans=~1/(1+exp(-.)),breaks=round(plogis(c(-3,-2,-1,0,1)),2),name='Relative proportion of dates after the eruption')) +
	labs(x='Average Thickness of Tephra Fall Deposit') +
	theme_minimal()

ashplot <- grid.arrange(ashplot1,ashplot2,ncol=1)

ggsave(here('figures_and_tables','post_scatter_ashfall.pdf'),plot=ashplot,width=4.75,height=7)


# Posterior Covariate Parameters  ----
params.a.beta1 <- data.frame(key='beta1',value=posterior.a2.750[,'beta1'])
params.a.beta2 <- data.frame(key='beta2',value=posterior.a2.750[,'beta2'])
params.a.beta12 <- data.frame(key='beta1+beta2',value=posterior.a2.750[,'beta1']+posterior.a2.750[,'beta2'])
params.a.kappa <- data.frame(key='kappa',value=10^posterior.a2.750[,'kappa'])

beta1.a <- ggplot(aes(x = value),data=params.a.beta1) + 
	stat_halfeye(slab_fill='steelblue',slab_alpha=0.5,scale=0.8,.width=c(0.9),point_interval="median_hdi",size=1.5,linewidth=1) +
	labs(x=TeX("$\\beta_1$"),y="Posterior Density",title=TeX('Model A:$\\beta_1$')) +
	geom_vline(aes(xintercept=0),linetype='dashed') +
	theme_minimal() +
	theme(plot.title=element_text(size=8),axis.title=element_text(size=6),axis.text=element_text(size=4)) 

beta2.a <- ggplot(aes(x = value),data=params.a.beta2) + 
	stat_halfeye(slab_fill='steelblue',slab_alpha=0.5,scale=0.8,.width=c(0.9),point_interval="median_hdi",size=1.5,linewidth=1) +
	labs(x=TeX("$\\beta_2$"),y="Posterior Density",title=TeX('Model A:$\\beta_2$')) +
	geom_vline(aes(xintercept=0),linetype='dashed') +
	theme_minimal() +
	theme(plot.title=element_text(size=8),axis.title=element_text(size=6),axis.text=element_text(size=4)) 

beta12.a <- ggplot(aes(x = value),data=params.a.beta12) + 
	stat_halfeye(slab_fill='steelblue',slab_alpha=0.5,scale=0.8,.width=c(0.9),point_interval="median_hdi",size=1.5,linewidth=1) +
	labs(x=TeX("$\\beta_1+beta_2$"),y="Posterior Density",title=TeX('Model A:$\\beta_1+beta_2$')) +
	geom_vline(aes(xintercept=0),linetype='dashed') +
	theme_minimal() +
	theme(plot.title=element_text(size=8),axis.title=element_text(size=6),axis.text=element_text(size=4)) 

kappa.a <- ggplot(aes(x = value),data=params.a.kappa) + 
	stat_halfeye(slab_fill='steelblue',slab_alpha=0.5,scale=0.8,.width=c(0.9),point_interval="median_hdi",size=1.5,linewidth=1) +
	labs(x=TeX("$\\kappa$"),y="Posterior Density",title=TeX('Model A:$\\kappa$')) +
	theme_minimal() +
	theme(plot.title=element_text(size=8),axis.title=element_text(size=6),axis.text=element_text(size=4)) 


params.a <- grid.arrange(beta1.a,beta2.a,beta12.a,kappa.a,ncol=2,nrow=2)
ggsave(here('figures_and_tables','params_posterior_a.pdf'),plot=params.a,width=7,height=5)



params.b.beta1 <- data.frame(key='beta1',value=posterior.b2.750[,'beta1'])
params.b.beta2 <- data.frame(key='beta2',value=posterior.b2.750[,'beta2'])
params.b.beta12 <- data.frame(key='beta1+beta2',value=posterior.b2.750[,'beta1']+posterior.b2.750[,'beta2'])
params.b.kappa <- data.frame(key='kappa',value=10^posterior.b2.750[,'kappa'])


beta1.b <- ggplot(aes(x = value),data=params.b.beta1) + 
	stat_halfeye(slab_fill='steelblue',slab_alpha=0.5,scale=0.8,.width=c(0.9),point_interval="median_hdi",size=1.5,linewidth=1) +
	labs(x=TeX("$\\beta_1$"),y="Posterior Density",title=TeX('Model A:$\\beta_1$')) +
	geom_vline(aes(xintercept=0),linetype='dashed') +
	theme_minimal() +
	theme(plot.title=element_text(size=8),axis.title=element_text(size=6),axis.text=element_text(size=4)) 

beta2.b <- ggplot(aes(x = value),data=params.b.beta2) + 
	stat_halfeye(slab_fill='steelblue',slab_alpha=0.5,scale=0.8,.width=c(0.9),point_interval="median_hdi",size=1.5,linewidth=1) +
	labs(x=TeX("$\\beta_2$"),y="Posterior Density",title=TeX('Model A:$\\beta_2$')) +
	geom_vline(aes(xintercept=0),linetype='dashed') +
	theme_minimal() +
	theme(plot.title=element_text(size=8),axis.title=element_text(size=6),axis.text=element_text(size=4)) 

beta12.b <- ggplot(aes(x = value),data=params.b.beta12) + 
	stat_halfeye(slab_fill='steelblue',slab_alpha=0.5,scale=0.8,.width=c(0.9),point_interval="median_hdi",size=1.5,linewidth=1) +
	labs(x=TeX("$\\beta_1+beta_2$"),y="Posterior Density",title=TeX('Model A:$\\beta_1+beta_2$')) +
	geom_vline(aes(xintercept=0),linetype='dashed') +
	theme_minimal() +
	theme(plot.title=element_text(size=8),axis.title=element_text(size=6),axis.text=element_text(size=4)) 

kappa.b <- ggplot(aes(x = value),data=params.b.kappa) + 
	stat_halfeye(slab_fill='steelblue',slab_alpha=0.5,scale=0.8,.width=c(0.9),point_interval="median_hdi",size=1.5,linewidth=1) +
	labs(x=TeX("$\\kappa$"),y="Posterior Density",title=TeX('Model A:$\\kappa$')) +
	theme_minimal() +
	theme(plot.title=element_text(size=8),axis.title=element_text(size=6),axis.text=element_text(size=4)) 


params.b <- grid.arrange(beta1.b,beta2.b,beta12.b,kappa.b,ncol=2,nrow=2)
ggsave(here('figures_and_tables','params_posterior_b.pdf'),plot=params.b,width=7,height=5)


## Posterior Summaries ----

model.a.posterior <- data.frame(params = colnames(posterior.a2.750),
				mean = apply(posterior.a2.750,2,mean),
				lo90 = apply(posterior.a2.750,2,function(x){coda::HPDinterval(mcmc(x),0.90)[1]}),
				hi90 = apply(posterior.a2.750,2,function(x){coda::HPDinterval(mcmc(x),0.90)[2]}),
				ess = ess.a2.750[!grepl('theta',names(ess.a2.750))],
				rhat = rhats.a2.750[[1]][,1][!grepl('theta',names(rhats.a2.750[[1]][,1]))])


write.csv(model.a.posterior,here('figures_and_tables','posteriors_model_a.csv'))


model.b.posterior <- data.frame(params = colnames(posterior.b2.750),
				mean = apply(posterior.b2.750,2,mean),
				lo90 = apply(posterior.b2.750,2,function(x){coda::HPDinterval(mcmc(x),0.90)[1]}),
				hi90 = apply(posterior.b2.750,2,function(x){coda::HPDinterval(mcmc(x),0.90)[2]}),
				ess = ess.b2.750[!grepl('theta',names(ess.b2.750))],
				rhat = rhats.b2.750[[1]][,1][!grepl('theta',names(rhats.b2.750[[1]][,1]))])


write.csv(model.b.posterior,here('figures_and_tables','posteriors_model_b.csv'))






model.b.posterior <- data.frame(params = colnames(posterior.b),
				mean = apply(posterior.b,2,mean),
				lo90 = apply(posterior.b,2,function(x){coda::HPDinterval(mcmc(x),0.90)[1]}),
				hi90 = apply(posterior.b,2,function(x){coda::HPDinterval(mcmc(x),0.90)[2]}),
				ess = ess.b,
				rhat = rhats.b[[1]][,1])
write.csv(model.b.posterior,here('figures_and_tables','posteriors_model_b.csv'))

