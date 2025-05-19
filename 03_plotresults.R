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

load(here('01_dataprep_out.RData'))
load(here('results','02a_fitmodel_out.RData'))
load(here('results','02b_fitmodel_out.RData'))
load(here('results','02c_fitmodel_out.RData'))
load(here('results','02d_fitmodel_out.RData'))

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

# Model C
post.r <- data.frame(par=c('r1','r1','r2','r2'),z=c(0,1,0,1)) 
r.matrix <- cbind(posterior.c[,'beta1'], posterior.c[,'beta1'] + posterior.c[,'eta1'],posterior.c[,'beta2'],posterior.c[,'beta2']+posterior.c[,'eta2'])
post.r$p_mean <- apply(r.matrix,2,mean) 
post.r$p_lo90 <- apply(r.matrix,2,function(x){HPDinterval(mcmc(x),prob=0.9)[1]}) 
post.r$p_hi90 <- apply(r.matrix,2,function(x){HPDinterval(mcmc(x),prob=0.9)[2]}) 

# Model D
post.beta <- data.frame(par=c('beta0','beta1'))
post.beta$p_mean <- c(mean(posterior.d[,'beta0']),mean(posterior.d[,'beta1']))
post.beta$p_lo90 <- c(HPDinterval(mcmc(posterior.d[,'beta0']),prob=0.9)[1],HPDinterval(mcmc(posterior.d[,'beta1']),prob=0.9)[1])
post.beta$p_hi90 <- c(HPDinterval(mcmc(posterior.d[,'beta0']),prob=0.9)[2],HPDinterval(mcmc(posterior.d[,'beta1']),prob=0.9)[2])
											    


# Map Figure Preparation ----

# keyregion ID
keys <- c(1,4,14,16,49,50)
key.letters  <- letters[length(keys)]
key_hexgird <- st_geometry(hexgrid_plot[keys,]) |> st_cast('MULTILINESTRING')

# combine to plot grid
r1_hexgrid <- left_join(hexgrid_plot,post.r1)
r2_hexgrid <- left_join(hexgrid_plot,post.r2)
delta_hexgrid <- left_join(hexgrid_plot,post.delta)
eta_hexgrid  <- left_join(hexgrid_plot,post.eta)

# binary field to determine whether a given grid has or does not have any sites
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
# sample_map <- ggplot(st_geometry(win)) +
# 	geom_sf(fill='darkgrey',colour='darkgrey') +
# 	geom_sf(data=st_geometry(sites),size=0.8) +
# 	geom_sf(data=st_geometry(hexgrid_plot),fill=NA,colour='lightgrey')+
# 	geom_contour(data=tephra,aes(x=x,y=y,z=value,linetype=factor(..level..)),breaks=c(100,200,500),color='black') +
# scale_linetype_manual(values=c('dotted','dashed','solid')) +
# 	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
# 	annotate('point',x=130.308,y=30.789,color='red',size=3,shape=17) +
# 	annotate('text',x=132,y=31,label='Kikai \n Caldera') +
# 	labs(x='Longitude',y='Laitude',linetype='Deposit Thickness (mm)') +
# 	theme(legend.position='inside',legend.position.inside=c(0.3,0.6))


sample_map <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour='darkgrey') +
	geom_sf(data=hexgrid_plot,aes(fill=ash),color='white',alpha=0.7)+
	scale_fill_viridis(option='turbo',trans='log10') +
	geom_sf(data=key_hexgird,color='black') +
	geom_sf(data=st_geometry(sites),size=0.5,col='grey20',alpha=0.8) +
#       geom_contour(data=tephra,aes(x=x,y=y,z=value),linetype='dashed',breaks=200,color='black') +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('point',x=130.308,y=30.789,fill='red',size=2,shape=24) +
	annotate('text',x=133.5,y=30.9,label='Kikai Caldera',size=3) +
	annotate('text',x=129.5,y=32.5,label='a') +
	annotate('text',x=129.5,y=31,label='b') +
	annotate('text',x=132,y=36,label='c') +
	annotate('text',x=134.7,y=36.3,label='d') +
	annotate('text',x=142,y=37.5,label='e') +
	annotate('text',x=142.5,y=39.5,label='f') +
	labs(x='Longitude',y='Laitude',fill='Average Deposit \nThickness (mm)') +
	theme(legend.position='inside',legend.position.inside=c(0.2,0.75),legend.background=element_rect(fill=alpha('white',0.5)),legend.key.size=unit(0.2,'in'),legend.text=element_text(size=5.5),legend.title=element_text(size=7))

ggsave(here('figures','sample_map.pdf'),plot=sample_map,width=4,height=4)


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

ggsave(here('figures','r1_icar_post.pdf'),plot=r1plot,width=7,height=5)



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
ggsave(here('figures','r2_icar_post.pdf'),plot=r2plot,width=7,height=5)

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
ggsave(here('figures','decrProb.pdf'),plot=decrProb,width=7,height=5)


# Posterior beta1 model A ----
post.beta1.modela.plot <- data.frame(key='beta1',value=posterior.a[,'beta1'])

#TODO : fix below
beta.r12 <- ggplot(aes(x = value),data=post.beta1.modela.plot) + 
	stat_halfeye(slab_fill='steelblue',slab_alpha=0.5,scale=0.8,.width=c(0.5,0.9),point_interval="median_hdi") +
	labs(x=TeX("$\\beta$"),y="Posterior Density") +
	geom_vline(aes(xintercept=0),linetype='dashed') +
	theme_minimal() 

ggsave(here('figures','beta_r12.pdf'),plot=beta.r12,width=5,height=4)

# Posterior eta model B ----
etaplot.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=eta_hexgrid,aes(fill=p_mean),color='darkgrey') +
	geom_sf(data=subset(eta_hexgrid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name=TeX('$\\eta$'),limits=minmax.eta) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label=TeX('Posterior Mean $\\eta$'),x=132,y=45) +
#TODO Could be done in latex2exp
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
ggsave(here('figures','eta_posterior.pdf'),plot=etaplot,width=7,height=5)


# HexSpecific Plot ----

ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour='darkgrey') +
	geom_sf(data=hexgrid_plot,alpha=0.7)+
	geom_sf_label(data=hexgrid_plot,aes(label=hexid)) +
#       geom_contour(data=tephra,aes(x=x,y=y,z=value),linetype='dashed',breaks=200,color='black') +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8))

