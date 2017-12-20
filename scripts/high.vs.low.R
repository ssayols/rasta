###################################
##
## plot the counts of the highly vs. lowly duplicated dataset 
##
###################################

pdf("high.vs.low2.pdf")
oldpar <- par(no.readonly = TRUE)

#layoutRatioWidth = c(0.75,0.25); layoutRatioHeight = c(0.25, 0.75)
#layout(matrix(c(2, 1, 0, 3), nrow = 2), widths = layoutRatioWidth, heights = layoutRatioHeight, respect = FALSE)
layoutRatioWidth = c(0.375,0.125,0.375,0.125); layoutRatioHeight = c(0.125, 0.375,0.125, 0.375)
layout(matrix(c(2,1,8,7,0,3,0,9,5,4,11,10,0,6,0,12),nrow=4),widths=layoutRatioWidth,heights=layoutRatioHeight,respect=FALSE)

f <- list.files("../results/rasta.good/",pattern="*.csv")
sapply(f,function(f) {
	x <- merge(read.csv(paste0("../results/rasta.good/",f),head=T),
			   read.csv(paste0("../results/rasta.bad/",f),head=T),
			   by="X",all=T)
	rownames(x) <- x[,1]; x <- as.matrix(x[,-1])

	fig <- function(x,y,xlab,ylab) {
		par(mar=c(4,4,1,1))
		smoothScatter(x,y,xlab=xlab,ylab=ylab,main=f,xlim=c(0,20),ylim=c(0,20))
		abline(a=0,b=1,lty=2,col="red")
		text
		par(mar=c(1,4,1,1))
		z <- density(x,na.rm=T)
		plot(z$x,z$y,xlab="",ylab="Density",xaxt="n",yaxt="n",type="l",main="")
		par(mar=c(4,1,1,1))
		z <- density(y,na.rm=T)
		plot(z$y,z$x,ylab="",xlab="Density",xaxt="n",yaxt="n",type="l",main="")
		invisible(0)
	}
	fig(log2(x[,1]),log2(x[,4]),xlab="log2(raw counts) low dupl",ylab="log2(raw counts) high dupl")
	fig(log2(x[,3]),log2(x[,6]),xlab="log2(rasta counts) low dupl",ylab="log2(rasta counts) high dupl")
	fig(log2(x[,1]),log2(x[,5]),xlab="log2(raw counts) low dupl",ylab="log2(censored counts) high dupl")
	fig(log2(x[,1]),log2(x[,6]),xlab="log2(raw counts) low dupl",ylab="log2(rasta counts) high dupl")
})

dev.off()
