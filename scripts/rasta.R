######################################
##
## Original RASTA implementation from http://www.stat.purdue.edu/~doerge/software/rasta.html
##
######################################
#Baumann and Doerge
#Bioinformatics June 4, 2013
#This code is for the analysis of RNA-seq data using the RASTA
#method to account for amplification bias

library(VGAM)
library(edgeR)
library(lattice)

#sample_data.txt has the following columns:
#  1. Read name (text)
#  2. Sample 1 read counts (numeric)
#  3. Sample 2 read counts (numeric)
#  4. Gene name (text)
#  5. Read start site (numeric)
#  6. Read end site (numeric)

sample_data <- read.table("sample_data.txt", header=FALSE)

#########################################################
## Compare the the digital read count to conventional 
## read count without amplification bias control

#Aggregate read counts by gene
samp2.reads.nocorrection <- aggregate(sample_data[,3], by=list(sample_data[,4]), sum)
samp1.reads.nocorrection <- aggregate(sample_data[,2], by=list(sample_data[,4]), sum)
data.nocorrection <- merge(samp2.reads.nocorrection, samp1.reads.nocorrection, by="Group.1")
#plot(log10(data.nocorrection[,3]), log10(data.nocorrection[,2]), col="blue")

#######################################################################
## Analyze the data with censoring at 1 (amplification bias correction)

censoring.data <- cbind(sample_data, rep(1, dim(sample_data)[1]))
cen.reads <- aggregate(censoring.data[,7], by=list(sample_data[,4]), sum)
data.censoring <- merge(data.nocorrection, cen.reads, by="Group.1")
#plot(log10(data.censoring[,3]), log10(data.censoring[,2]), col="black")
#points(log10(data.censoring[,3]), log10(data.censoring[,4]), col="red")


#######################################################################
## Analyze the data with RASTA (amplification bias correction)
# The read count list needs to be reformed for this analysis
 
data.temp <- split(sample_data, sample_data[,4])
data.reads <- lapply(data.temp, function(x) x[,3])
date()
data.counts <- matrix(nrow=length(data.reads), ncol=2)
for(i in 1:length(data.reads)){
	test1 <- data.reads[[i]]
	if(length(test1) > 1) {
	temp1 <- try(hc <- hclust(dist(unique(test1),method="canberra"),method="complete"),silent=TRUE)
	temp2 <- try(ct <- cutree(hc,k=2),silent=TRUE)
	temp3 <- try(cmeans <- c(mean(test1[which(ct==1)]), mean(test1[which(ct==2)])),silent=TRUE)
	temp4 <- try(min.group <- which(cmeans==min(cmeans)),silent=TRUE)
	temp5 <- try(fit <- vglm(test1[which(ct==min.group)] ~ 1, pospoisson),silent=TRUE)
	temp6 <- try(pospois.rate <- Coef(fit),silent=TRUE)
				
	if(class(temp1)=="try-error"){pospois.rate <- Coef(vglm(test1~1,pospoisson))}
	test2 <- NULL	
	for(j in 1:length(test1)){
			test2[j] <- min(test1[j], max(1,qpois(.95,pospois.rate)))
		}
	data.counts[i,1] <- as.character(names(data.reads)[[i]])
	data.counts[i,2] <- as.numeric(sum(test2))
	} else {
	data.counts[i,1] <- as.character(names(data.reads)[[i]])
	data.counts[i,2] <- as.numeric(test1)
	}
}
date()



