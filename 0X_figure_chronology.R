library(here)
library(plotrix)
pdf(file=here('figures_and_tables','figure_chronology.pdf'),width=4.75,height=3.5)
par(mar=c(3,1,1,1))
plot(NULL,xlim=c(11000,5500),ylim=c(0,7),xlab='',ylab='',axes=F)
axis(1, at=seq(11000, 5000, -1000), labels=paste(seq(11, 5, -1), 'k', sep=''), cex.axis=0.87, mgp=c(3, 0.5, 0))
axis(1,at=seq(11000,5000,-500),labels=F,tck=-0.01)
axis(1,at=seq(11000,5000,-100),labels=F,tck=-0.005)
mtext('cal BP',side=1,line=1,at=5500,cex=0.7)
rect(xleft=11600,xright=9300,ybottom=-1,ytop=11,col=adjustcolor('lightblue',0.4),border=NA)
#rect(xleft=10800,xright=10500,ybottom=-1,ytop=11,col=adjustcolor('lightblue',1),border=NA)
rect(xleft=9200,xright=7800,ybottom=-1,ytop=11,col=adjustcolor('orange',0.4),border=NA)
rect(xleft=8000,xright=6500,ybottom=0,ytop=7,col=NULL,lwd=1.5,lty=3)
abline(v=7250,lty=2,col='indianred',lwd=1.5)
text(7200,6.2,'K-Ah eruption',cex=0.6,col='indianred',pos=1,srt=90)

text(6850,2,'Window of analysis \n (current paper)',cex=0.6,pos=1,srt=90)


lines(c(11600,9300),c(1,1),lwd=2) #Kaigara-Entomon
lines(c(10170,8950),c(1.5,1.5),lwd=2) #Oshigatamon
lines(c(9180,8600),c(2,2),lwd=2) #Hiragakoi
lines(c(9060,8500),c(2.5,2.5),lwd=2) #Senokan A
lines(c(8400,7790),c(3,3),lwd=2) #Senokan B
lines(c(7720,7410),c(3.5,3.5),lwd=2) #Kubama
lines(c(7600,7300),c(4,4),lwd=2) #Todoroki A
lines(c(7350,7050),c(4.5,4.5),lwd=2) #Nishinosono
lines(c(7100,6900),c(5,5),lwd=2) #Todoroki B1
lines(c(6900,6500),c(5.5,5.5),lwd=2) #Todoroki B2
lines(c(6700,6300),c(6,6),lwd=2) #Todoroki B3
lines(c(6300,5500),c(6.5,6.5),lwd=2) #Sobata

text(11200,6.5,'Villages \nrise and demise',cex=0.6,pos=4)
text(9200,6.5,'Aggregation sites \nrise and demise',cex=0.6,pos=4)
# text(10700,5,'Large circular settlements (peak)',cex=0.6,pos=1,srt=90)
off=-0.2
text(x=median(c(11600,9300)),y=1+off,'Kaigara-Entomon',cex=0.5,pos=4)
text(x=median(c(10170,8950)),y=1.5+off,'Oshigatamon',cex=0.5,pos=4)
text(x=median(c(9180,8600)),y=2.0+off,'Hiragakoi',cex=0.5,pos=4)
text(x=median(c(9060,8500)),y=2.5+off,'Senokan A',cex=0.5,pos=4)
text(x=median(c(8400,7790)),y=3.0+off,'Senokan B',cex=0.5,pos=4)
text(x=median(c(7720,7410)),y=3.5+off,'Kubama',cex=0.5,pos=4)
text(x=median(c(7600,7300)),y=4.0+off,'Todoroki A',cex=0.5,pos=4)
text(x=median(c(7350,7050)),y=4.5+off,'Nishinosono',cex=0.5,pos=4)
text(x=median(c(7100,6900)),y=5.0+off,'Todoroki B1',cex=0.5,pos=4)
text(x=median(c(6900,6500)),y=5.5+off,'Todoroki B2',cex=0.5,pos=4)
text(x=median(c(6700,6300)),y=6.0+off,'Todoroki B3',cex=0.5,pos=4)
text(x=median(c(6300,5500)),y=6.5+off,'Sobata',cex=0.5,pos=4)
dev.off()
