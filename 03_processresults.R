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

load(here('01_dataprep_out.RData'))
load(here('results','02a_fitmodel_out.RData'))
load(here('results','02b_fitmodel_out.RData'))

# Post-processing Posteriors ----
# Model A 
# Extract posterior summary stats
post.r1 <- post.r2 <- post.delta  <- data.frame(hexid=1:constants.icar$n.hex)
r1.matrix <- posterior.a[,grep('s1',colnames(posterior.a))] 
r2.matrix <- posterior.a[,grep('s2',colnames(posterior.a))] +  posterior.a[,'beta1'] %*% t(constants.icar$ash) #check if appropriately calculated

post.r1$p_mean <- apply(r1.matrix,2,mean) * 100
post.r1$p_lo90 <- apply(r1.matrix,2,function(x){HPDinterval(mcmc(x),prob=0.9)[1]}) * 100
post.r1$p_hi90 <- apply(r1.matrix,2,function(x){HPDinterval(mcmc(x),prob=0.9)[2]}) * 100
post.r1$prob_decrease  <- apply(r1.matrix,2,function(x){sum(x<0)/length(x)})
post.r2$p_mean <- apply(r2.matrix,2,mean) * 100
post.r2$p_lo90 <- apply(r2.matrix,2,function(x){HPDinterval(mcmc(x),prob=0.9)[1]}) * 100
post.r2$p_hi90 <- apply(r2.matrix,2,function(x){HPDinterval(mcmc(x),prob=0.9)[2]}) * 100
post.r2$prob_decrease  <- apply(r2.matrix,2,function(x){sum(x<0)/length(x)})

minmax.r12  <- range(c(post.r1$p_lo90,post.r2$p_lo90,post.r1$p_hi90,post.r2$p_hi90)) |> round(1)
minmax.r12 <- c(-1.3,0.5)

# Model B
post.eta <- data.frame(hexid=1:constants.icar$n.hex)
eta.matrix <- posterior.b[,grep('^eta',colnames(posterior.b))] +  posterior.b[,'beta1'] %*% t(constants.icar$ash)
post.eta$p_mean <- apply(eta.matrix,2,mean) 
post.eta$p_lo90 <- apply(eta.matrix,2,function(x){HPDinterval(mcmc(x),prob=0.9)[1]}) 
post.eta$p_hi90 <- apply(eta.matrix,2,function(x){HPDinterval(mcmc(x),prob=0.9)[2]}) 
minmax.eta  <- range(c(post.eta$p_lo90,post.eta$p_hi90))


# Map Figure Preparation ----

# keyregion ID
keys <- c(1,4,14,16,21,49,50)
key.letters  <- letters[1:length(keys)]
key_hexgrid <- st_geometry(hexgrid_plot[keys,]) |> st_cast('MULTILINESTRING')
r1r2keys  <- cbind(r1.matrix[,paste0('s1[',keys,']')]*100,r2.matrix[,paste0('s2[',keys,']')]*100)  
r1r2keys <- gather(r1r2keys)
r1r2keys$param <- 'r1'
r1r2keys$param[grep('s2',r1r2keys$key)] <- 'r2'
r1r2keys$hex <- key.letters[match(as.integer(sub(".*\\[(\\d+)\\].*", "\\1", r1r2keys$key)),keys)]

etakeys  <- posterior.b[,paste0('eta[',keys,']')] + posterior.b[,'beta1'] %*% t(constants.icar$ash[keys])
etakeys <- gather(etakeys)
etakeys$hex <- key.letters[match(as.integer(sub(".*\\[(\\d+)\\].*", "\\1", etakeys$key)),keys)]

# combine to plot grid
r1_hexgrid <- left_join(hexgrid_plot,post.r1)
r2_hexgrid <- left_join(hexgrid_plot,post.r2)
delta_hexgrid <- left_join(hexgrid_plot,post.delta)
eta_hexgrid  <- left_join(hexgrid_plot,post.eta)

mt <- st_intersects(hexgrid_plot,st_as_sf(sites.df,coords=c('Easting','Northing'),crs=6684))
eta_hexgrid$obs <- delta_hexgrid$obs <- r1_hexgrid$obs <- r2_hexgrid$obs <- lengths(mt) > 0

# convert to lat/long
r1_hexgrid <- st_transform(r1_hexgrid,crs=4326)
r2_hexgrid <- st_transform(r2_hexgrid,crs=4326)
delta_hexgrid <- st_transform(delta_hexgrid,crs=4326)
eta_hexgrid <- st_transform(eta_hexgrid,crs=4326)

sites <- st_as_sf(sites.df,coords=c('Easting','Northing'),crs=6684) |> st_transform(crs=4326)
win  <- ne_countries(scale=10,returnclass='sf') |> st_combine() #download background map

# Read ashfall data
tephra  <- rast(here('data','K_Ah_SUM.tif'))
tephra <- project(tephra,"EPSG:4326")
tephra <- as.data.frame(tephra,xy=T)
tephra[tephra==0] <- NA
colnames(tephra)[3] <- 'value'
hexgrid_plot$ash <- hexgrid$ash

# Make Figures ----

# Sample Distribution Map ----
sample_map <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour='darkgrey') +
	geom_sf(data=hexgrid_plot,aes(fill=ash),color='white',alpha=0.7)+
	scale_fill_viridis(option='turbo',trans='log10') +
	geom_sf(data=key_hexgrid,color='black') +
	geom_sf(data=st_geometry(sites),size=0.5,col='grey20',alpha=0.8) +
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
	labs(x='Longitude',y='Laitude',fill='Average Deposit \nThickness (mm)') +
	theme(legend.position='inside',legend.position.inside=c(0.2,0.75),legend.background=element_rect(fill=alpha('white',0.5)),legend.key.size=unit(0.2,'in'),legend.text=element_text(size=5.5),legend.title=element_text(size=7))

ggsave(here('figures_and_tables','sample_map.pdf'),plot=sample_map,width=4,height=4)

# K-Akahoya Ashfall model ----
ashfall_model_raw <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour='darkgrey') +
	geom_contour_filled(data=tephra,aes(x=x,y=y,z=value),alpha=0.5,breaks=c(0,150,300,500,1000,1500,2000,2500,3000)) +
	scale_fill_viridis_d(option='turbo') +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('point',x=130.308,y=30.789,fill='red',size=2,shape=24) +
	annotate('text',x=133.5,y=30.9,label='Kikai Caldera',size=3) +
	labs(fill='Tephra thickness (in mm)') +
	theme(legend.position='inside',legend.position.inside=c(0.25,0.75),legend.background=element_rect(fill=alpha('white',0.5)),legend.key.size=unit(0.15,'in'),legend.text=element_text(size=5),legend.title=element_text(size=6))

ggsave(here('figures_and_tables','ashfall_model_raw.pdf'),plot=ashfall_model_raw,width=4,height=4)

# r1 posterior model A ----

r1plot.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=r1_hexgrid,aes(fill=p_mean),color='darkgrey') +
	geom_sf(data=subset(r1_hexgrid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.r12) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX("Posterior Mean $r_1$"),x=132,y=45) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())


r1plot.lo <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=r1_hexgrid,aes(fill=p_lo90),color='darkgrey') +
	geom_sf(data=subset(r1_hexgrid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.r12) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(lower)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 

r1plot.hi <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=r1_hexgrid,aes(fill=p_hi90),color='darkgrey') +
	geom_sf(data=subset(r1_hexgrid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.r12) +
	# 	geom_sf(data=st_geometry(sites),size=0.8) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(upper)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 


r1plot <- grid.arrange(r1plot.main,r1plot.lo,r1plot.hi,layout_matrix = rbind(c(1,2),c(1,3)),widths=c(2,1),heights=c(1,1),padding=unit(0.2,'cm'))

ggsave(here('figures_and_tables','posterior_r1.pdf'),plot=r1plot,width=7,height=5)

# r2 posterior model A ----

r2plot.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=r2_hexgrid,aes(fill=p_mean),color='darkgrey') +
	geom_sf(data=subset(r2_hexgrid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.r12) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX('Posterior Mean $r_2$'),x=132,y=45) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())


r2plot.lo <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=r2_hexgrid,aes(fill=p_lo90),color='darkgrey') +
	geom_sf(data=subset(r2_hexgrid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.r12) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(lower)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 

r2plot.hi <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=r2_hexgrid,aes(fill=p_hi90),color='darkgrey') +
	geom_sf(data=subset(r2_hexgrid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.r12) +
	# 	geom_sf(data=st_geometry(sites),size=0.8) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(upper)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 


r2plot <- grid.arrange(r2plot.main,r2plot.lo,r2plot.hi,layout_matrix = rbind(c(1,2),c(1,3)),widths=c(2,1),heights=c(1,1),padding=unit(0.2,'cm'))
ggsave(here('figures_and_tables','posterior_r2.pdf'),plot=r2plot,width=7,height=5)

# Probability of Decrease model A----

decrProb.left <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=r1_hexgrid,aes(fill=prob_decrease),color='darkgrey') +
	geom_sf(data=subset(r1_hexgrid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='darkgreen',mid='white',high='darkorange',midpoint=0.5,name='Probability Negative Growth',limits=c(0,1)) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='Before Eruption',x=132,y=45) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(1, 1, 1, 1),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())

decrProb.right <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=r2_hexgrid,aes(fill=prob_decrease),color='darkgrey') +
	geom_sf(data=subset(r2_hexgrid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='darkgreen',mid='white',high='darkorange',midpoint=0.5,name='Probability Negative Growth',limits=c(0,1)) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='After Eruption',x=132,y=45) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='none',plot.margin = margin(1, 1, 1, 1),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())

decrProb <- grid.arrange(decrProb.left,decrProb.right,ncol=2,padding=unit(0.2,'cm'))
ggsave(here('figures_and_tables','decrease_probability.pdf'),plot=decrProb,width=7,height=5)

# Posterior eta model B ----
etaplot.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=eta_hexgrid,aes(fill=p_mean),color='darkgrey') +
	geom_sf(data=subset(eta_hexgrid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name=TeX('$\\eta$'),limits=minmax.eta) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX('Posterior Mean $\\eta$'),x=132,y=45) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())


etaplot.lo <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=eta_hexgrid,aes(fill=p_lo90),color='darkgrey') +
	geom_sf(data=subset(eta_hexgrid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.eta) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(lower)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 

etaplot.hi <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=eta_hexgrid,aes(fill=p_hi90),color='darkgrey') +
	geom_sf(data=subset(eta_hexgrid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.eta) +
	# 	geom_sf(data=st_geometry(sites),size=0.8) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(upper)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 


etaplot <- grid.arrange(etaplot.main,etaplot.lo,etaplot.hi,layout_matrix = rbind(c(1,2),c(1,3)),widths=c(2,1),heights=c(1,1),padding=unit(0.2,'cm'))
ggsave(here('figures_and_tables','posterior_eta.pdf'),plot=etaplot,width=7,height=5)


# Posterior Means of the two models
r1mean <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=r1_hexgrid,aes(fill=p_mean),color='darkgrey') +
	geom_sf(data=subset(r1_hexgrid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.r12) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX("Posterior Mean $r_1$"),x=133,y=45,size=5) +
	labs(x='',y='') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=8),legend.title=element_text(size=10),legend.key.size=unit(0.25,'in'),legend.background=element_blank(),axis.title=element_text(size=6),plot.title=element_text(size=6),axis.text=element_text(size=8))

r2mean <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=r2_hexgrid,aes(fill=p_mean),color='darkgrey') +
	geom_sf(data=subset(r2_hexgrid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.r12) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	labs(x='',y='') +
	annotate('text',label=TeX("Posterior Mean $r_2$"),x=133,y=45,size=5) +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=8),legend.title=element_text(size=10),legend.key.size=unit(0.25,'in'),legend.background=element_blank(),axis.title=element_text(size=6),plot.title=element_text(size=6),axis.text=element_text(size=8))

etamean <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=eta_hexgrid,aes(fill=p_mean),color='darkgrey') +
	geom_sf(data=subset(eta_hexgrid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name=TeX('$\\eta$'),limits=minmax.eta) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	labs(x='',y='') +
	annotate('text',label=TeX('Posterior Mean $\\eta$'),x=133,y=45,size=5) +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=8),legend.title=element_text(size=10),legend.key.size=unit(0.25,'in'),legend.background=element_blank(),axis.title=element_text(size=6),plot.title=element_text(size=6),axis.text=element_text(size=8))

posterior.means.combined  <- grid.arrange(r1mean,r2mean,etamean,ncol=3)
ggsave(here('figures_and_tables','mean_posterior.pdf'),plot=posterior.means.combined,width=7.25,height=5)

# Posterior beta  ----
post.beta1.r12.plot <- data.frame(key='beta1',value=posterior.a[,'beta1'])
beta.r12 <- ggplot(aes(x = value),data=post.beta1.r12.plot) + 
	stat_halfeye(slab_fill='steelblue',slab_alpha=0.5,scale=0.8,.width=c(0.9),point_interval="median_hdi",size=1.5,linewidth=1) +
	labs(x=TeX("$\\beta$"),y="Posterior Density",title='Model a') +
	geom_vline(aes(xintercept=0),linetype='dashed') +
	scale_x_continuous(breaks=c(0,-0.00001,-0.00002),labels=c('0',TeX('$1 \\times 10^{-5}$'),TeX('$2 \\times 10^{-5}$'))) +
	theme_minimal() +
	theme(plot.title=element_text(size=5),axis.title=element_text(size=5),axis.text=element_text(size=4)) 

post.beta1.eta.plot <- data.frame(key='beta1',value=posterior.b[,'beta1'])
beta.eta <- ggplot(aes(x = value),data=post.beta1.eta.plot) + 
	stat_halfeye(slab_fill='steelblue',slab_alpha=0.5,scale=0.8,.width=c(0.9),point_interval="median_hdi",size=1.5,linewidth=1) +
	labs(x=TeX("$\\beta$"),y="Posterior Density",title='Model b') +
	geom_vline(aes(xintercept=0),linetype='dashed') +
	theme_minimal() +
	theme(plot.title=element_text(size=5),axis.title=element_text(size=5),axis.text=element_text(size=4)) 

beta.combined  <- grid.arrange(beta.r12,beta.eta,ncol=1)

ggsave(here('figures_and_tables','beta_posterior.pdf'),plot=beta.combined,width=2.25,height=3)

# HexSpecific Plot ----
r1r2keyplot <- ggplot(r1r2keys,aes(x=value,col=param,y=hex)) +
	stat_pointinterval(position=position_dodge(0.3),.width=c(.5,.75,.9),point_interval='median_hdi') +
	scale_color_manual(values=c('r1'='darkorange','r2'='darkgreen'),labels=c('r1'=TeX('$r_1$'),'r2'=TeX('$r_2$'))) +
	theme_minimal() +
	xlim(-1,0.5) +
	labs(x='Annual Growth Rate',y='Hexagon',col='Parameters') +
	geom_vline(xintercept=0,linetype='dashed') +
	theme(legend.position = 'inside', legend.position.inside=c(0.2,0.8))


etakeyplot <- ggplot(etakeys,aes(x=value,y=hex)) + 
	stat_pointinterval(position=position_dodge(0.3),.width=c(.5,.75,.9),point_interval='median_hdi') +
	theme_minimal() +
	geom_vline(xintercept=0,linetype='dashed') +
	scale_x_continuous(name=TeX('$\\eta$'),sec.axis=sec_axis(transform=~1/(1+exp(-.)),breaks=round(plogis(c(-3,-2,-1,0,1)),2),name='Relative proportion of dates after the eruption')) +
	labs(x=TeX('$\\eta$'),y='Hexagon')


hex_combined <- grid.arrange(r1r2keyplot,etakeyplot,ncol=1)
ggsave(here('figures_and_tables','hex_focus_plot.pdf'),plot=hex_combined,width=4.75,height=7)


# Posterior parameters vs predictors ----
post.r2$ash <- constants.icar$ash
post.r2$keys <- NA
post.r2$keys[keys] <- key.letters
post.eta$ash <- constants.icar$ash
post.eta$keys <- NA
post.eta$keys[keys] <- key.letters

ashplot1 <- ggplot(data=post.r2,aes(x=ash,y=p_mean)) +
	geom_errorbar(aes(ymin=p_lo90,ymax=p_hi90),alpha=0.5) +
	geom_point() +
	geom_vline(xintercept=150,linetype='dashed') +
	geom_hline(yintercept=0,linetype='dotted') +
	scale_x_log10() +
	labs(x='Average Thickness of Tephra Fall Deposit',y='Posterior Estimate of Annual Growth Rate (%)') +
	theme_minimal()

ashplot2 <- ggplot(data=post.eta,aes(x=ash,y=p_mean)) +
	geom_errorbar(aes(ymin=p_lo90,ymax=p_hi90),alpha=0.5) +
	geom_point() +
	geom_vline(xintercept=150,linetype='dashed') +
	geom_hline(yintercept=0,linetype='dotted') +
	scale_x_log10() +
	scale_y_continuous(name=TeX('$\\eta$'),sec.axis=sec_axis(trans=~1/(1+exp(-.)),breaks=round(plogis(c(-3,-2,-1,0,1)),2),name='Relative proportion of dates after the eruption')) +
	labs(x='Average Thickness of Tephra Fall Deposit') +
	theme_minimal()

ashplot <- grid.arrange(ashplot1,ashplot2,ncol=1)

ggsave(here('figures_and_tables','post_scatter_ashfall.pdf'),plot=ashplot,width=4.75,height=7)

## Posterior Summaries ----

model.a.posterior <- data.frame(params = colnames(posterior.a),
				mean = apply(posterior.a,2,mean),
				lo90 = apply(posterior.a,2,function(x){coda::HPDinterval(mcmc(x),0.90)[1]}),
				hi90 = apply(posterior.a,2,function(x){coda::HPDinterval(mcmc(x),0.90)[2]}),
				ess = ess.a,
				rhat = rhats.a[[1]][,1])

write.csv(model.a.posterior,here('figures_and_tables','posteriors_model_a.csv'))

model.b.posterior <- data.frame(params = colnames(posterior.b),
				mean = apply(posterior.b,2,mean),
				lo90 = apply(posterior.b,2,function(x){coda::HPDinterval(mcmc(x),0.90)[1]}),
				hi90 = apply(posterior.b,2,function(x){coda::HPDinterval(mcmc(x),0.90)[2]}),
				ess = ess.b,
				rhat = rhats.b[[1]][,1])
write.csv(model.b.posterior,here('figures_and_tables','posteriors_model_b.csv'))

