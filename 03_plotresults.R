library(here)
library(ggplot2)
library(sf)
library(coda)
library(dplyr)
library(rnaturalearth)
library(gridExtra)
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

post.delta.raw <- (r2.matrix - r1.matrix) * 100
post.delta$p_mean <- apply(post.delta.raw,2,mean) 
post.delta$p_lo90 <- apply(post.delta.raw,2,function(x){HPDinterval(mcmc(x),prob=0.9)[1]}) 
post.delta$p_hi90 <- apply(post.delta.raw,2,function(x){HPDinterval(mcmc(x),prob=0.9)[2]}) 

minmax.r12.a  <- range(c(post.r1$p_lo90,post.r2$p_lo90,post.r1$p_hi90,post.r2$p_hi90))
minmax.delta.a <- range(c(post.delta$p_lo,post.delta$p_hi))

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
post.beta$p_lo90 <- HPDinterval(mcmc(posterior.d[,'beta0']),prob=0.9)[1]
post.beta$p_hi90 <- HPDinterval(mcmc(posterior.d[,'beta1']),prob=0.9)[2] 





# combine to plot grid
pre_hex_grid <- left_join(hexgrid_plot,post.r1)
post_hex_grid <- left_join(hexgrid_plot,post.r2)
delta_hex_grid <- left_join(hexgrid_plot,post.delta)
# binary field to determine whether a given grid has or does not have any sites
mt <- st_intersects(hexgrid_plot,st_as_sf(sites.df,coords=c('Easting','Northing'),crs=6684))
delta_hex_grid$obs <- pre_hex_grid$obs <- post_hex_grid$obs <- lengths(mt) > 0

# convert to lat/long
pre_hex_grid <- st_transform(pre_hex_grid,crs=4326)
post_hex_grid <- st_transform(post_hex_grid,crs=4326)
delta_hex_grid <- st_transform(delta_hex_grid,crs=4326)
sites <- st_as_sf(sites.df,coords=c('Easting','Northing'),crs=6684) |> st_transform(crs=4326)
win  <- ne_countries(scale=10,returnclass='sf') |> st_combine() #download background map

# Make Figures ----

# Figure 1 (Sample Distribution Map) ----
fig1 <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour='darkgrey') +
	geom_sf(data=st_geometry(sites),size=0.8) +
	geom_sf(data=st_geometry(hexgrid_plot),fill=NA,colour='lightgrey')+
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8))

ggsave(here('figures','figure1.pdf'),plot=fig1,width=4,height=4)
# Figure 2 (Posterior before) ----#

fig2.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=pre_hex_grid,aes(fill=p_mean),color='darkgrey') +
	geom_sf(data=subset(pre_hex_grid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.growth) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='Posterior Mean',x=132,y=45) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())


fig2.lo <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=pre_hex_grid,aes(fill=p_lo90),color='darkgrey') +
	geom_sf(data=subset(pre_hex_grid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.growth) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(lower)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 

fig2.hi <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=pre_hex_grid,aes(fill=p_hi90),color='darkgrey') +
	geom_sf(data=subset(pre_hex_grid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.growth) +
	# 	geom_sf(data=st_geometry(sites),size=0.8) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(upper)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 


fig2 <- grid.arrange(fig2.main,fig2.lo,fig2.hi,layout_matrix = rbind(c(1,2),c(1,3)),widths=c(2,1),heights=c(1,1),padding=unit(0.2,'cm'))
ggsave(here('figures','figure2.pdf'),plot=fig2,width=7,height=5)

# Figure 3 (Posterior after) ----

fig3.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=post_hex_grid,aes(fill=p_mean),color='darkgrey') +
	geom_sf(data=subset(post_hex_grid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.growth) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='Posterior Mean',x=132,y=45) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())


fig3.lo <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=post_hex_grid,aes(fill=p_lo90),color='darkgrey') +
	geom_sf(data=subset(post_hex_grid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.growth) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(lower)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 

fig3.hi <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=post_hex_grid,aes(fill=p_hi90),color='darkgrey') +
	geom_sf(data=subset(post_hex_grid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Annual Growth Rate (%)',limits=minmax.growth) +
	# 	geom_sf(data=st_geometry(sites),size=0.8) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(upper)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 


fig3 <- grid.arrange(fig3.main,fig3.lo,fig3.hi,layout_matrix = rbind(c(1,2),c(1,3)),widths=c(2,1),heights=c(1,1),padding=unit(0.2,'cm'))
ggsave(here('figures','figure3.pdf'),plot=fig3,width=7,height=5)

# Figure 4 (Probability Decrease) ----
fig4.left <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=pre_hex_grid,aes(fill=prob_decrease),color='darkgrey') +
	geom_sf(data=subset(pre_hex_grid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='darkgreen',mid='white',high='darkorange',midpoint=0.5,name='Probability Negative Growth',limits=c(0,1)) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='Before Eruption',x=132,y=45) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(1, 1, 1, 1),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())

fig4.right <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=post_hex_grid,aes(fill=prob_decrease),color='darkgrey') +
	geom_sf(data=subset(post_hex_grid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='darkgreen',mid='white',high='darkorange',midpoint=0.5,name='Probability Negative Growth',limits=c(0,1)) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='After Eruption',x=132,y=45) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='none',plot.margin = margin(1, 1, 1, 1),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())

fig4 <- grid.arrange(fig4.left,fig4.right,ncol=2,padding=unit(0.2,'cm'))
ggsave(here('figures','figure4.pdf'),plot=fig4,width=7,height=5)



# Figure 5 (Growth Rate Difference) ----

fig5.main <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=delta_hex_grid,aes(fill=p_mean),color='darkgrey') +
	geom_sf(data=subset(delta_hex_grid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Delta Growth Rate (%)',limits=minmax.delta) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='Posterior Mean',x=132,y=45) +
	xlab('Longitude') +
	ylab('Latitude') +
	theme(legend.position='inside',legend.position.inside=c(0.3,0.6),plot.margin = margin(0, 0, 0, 0),legend.text=element_text(size=7),legend.title=element_text(size=8),legend.key.size=unit(0.2,'in'),legend.background=element_blank())


fig5.lo <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=delta_hex_grid,aes(fill=p_lo90),color='darkgrey') +
	geom_sf(data=subset(delta_hex_grid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Delta Growth Rate (%)',limits=minmax.delta) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(lower)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 

fig5.hi <- ggplot(st_geometry(win)) +
	geom_sf(fill='darkgrey',colour=NA) +
	geom_sf(data=delta_hex_grid,aes(fill=p_hi90),color='darkgrey') +
	geom_sf(data=subset(delta_hex_grid,obs==FALSE),fill='darkgrey',color='darkgrey') +
	scale_fill_gradient2(low='blue',mid='white',high='red',midpoint=0,name='Delta Growth Rate (%)',limits=minmax.delta) +
	# 	geom_sf(data=st_geometry(sites),size=0.8) +
	coord_sf(xlim=c(129.5,145),ylim=c(31,45.8)) +
	annotate('text',label='90% HPD \n(upper)',x=133,y=45,size=3) +
	theme(legend.position='none', axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank() , axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank(),plot.margin = margin(0, 0.5, 0, 0)) 


fig5 <- grid.arrange(fig5.main,fig5.lo,fig5.hi,layout_matrix = rbind(c(1,2),c(1,3)),widths=c(2,1),heights=c(1,1),padding=unit(0.2,'cm'))
ggsave(here('figures','figure5.pdf'),plot=fig5,width=7,height=5)

