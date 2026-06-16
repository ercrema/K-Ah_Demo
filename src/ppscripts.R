# WAIC processing ----
deltaWAIC  <- function(x){return(x - min(x))}
WAICweights <- function(x){return(exp(-0.5*deltaWAIC(x))/sum(exp(-0.5*deltaWAIC(x))))}	

# Summarise posteriors for model A ICAR
post_summary_a <- function(r1,r2,n,p=0.9)
{
	res <- data.frame(hexid=1:n)
	res$r1.p_mean <- apply(r1,2,mean)*100
	res$r1.p_lo <- apply(r1,2,function(x,p){HPDinterval(mcmc(x),prob=p)[1]},p=p) * 100
	res$r1.p_hi <- apply(r1,2,function(x,p){HPDinterval(mcmc(x),prob=p)[2]},p=p) * 100
	res$r1.p_decrease <- apply(r1,2,function(x){sum(x<0)/length(x)})
	res$r2.p_mean <- apply(r2,2,mean)*100
	res$r2.p_lo <- apply(r2,2,function(x,p){HPDinterval(mcmc(x),prob=p)[1]},p=p) * 100
	res$r2.p_hi <- apply(r2,2,function(x,p){HPDinterval(mcmc(x),prob=p)[2]},p=p) * 100
	res$r2.p_decrease <- apply(r2,2,function(x){sum(x<0)/length(x)})
	res$dr2r1.p_mean <- apply(r2-r1,2,mean)*100
	res$dr2r1.p_lo <- apply(r2-r1,2,function(x,p){HPDinterval(mcmc(x),prob=p)[1]},p=p) * 100
	res$dr2r1.p_hi <- apply(r2-r1,2,function(x,p){HPDinterval(mcmc(x),prob=p)[2]},p=p) * 100
	res$dr2r1.p_decrease <- apply(r2-r1,2,function(x){sum(x<0)/length(x)})
	minmaxr1r2.lohi <- range(c(res$r1.p_lo,res$r1.p_hi,res$r2.p_lo,res$r2.p_hi)) |> round(1)
	minmaxr1r2.mean <- range(c(res$r1.p_mean,res$r2.p_mean))
	return(list(post=res,minmax.lohi=minmaxr1r2.lohi,minmax.mean=minmaxr1r2.mean))
}


	
# Summarise posteriors for model B ICAR
post_summary_b <- function(eta,n,p=0.9)
{
	res <- data.frame(hexid=1:n)
	res$p_mean <- apply(eta,2,mean) 
	res$p_lo90 <- apply(eta,2,function(x,p){HPDinterval(mcmc(x),prob=p)[1]},p=p) 
	res$p_hi90 <- apply(eta,2,function(x,p){HPDinterval(mcmc(x),prob=p)[2]},p=p) 
	res$p_decrease <- apply(eta,2,function(x){sum(x<0)/length(x)})
	minmax.eta.lohi  <- range(c(res$p_lo90,res$p_hi90))
	minmax.eta.mean  <- range(c(res$p_mean))
	return(list(post=res,minmax.lohi=minmax.eta.lohi,minmax.mean=minmax.eta.mean))
}

